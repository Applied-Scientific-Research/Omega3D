/*
 * Surfaces.h - Specialized class for trimesh surfaces in 3D
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega3D.h"
#include "MathHelper.h"
#include "VectorHelper.h"
#include "ElementBase.h"

#ifdef USE_GL
#include "GlState.h"
#include "RenderParams.h"
#include "OglHelper.h"
#include "ShaderHelper.h"
#include "glad.h"
#endif

#include <iostream>
#include <vector>
#include <array>
#include <algorithm> // for max_element
#include <optional>
#include <cassert>


// useful structure for panel strengths and BCs
// Strength<S> ps;
//   ps[0] means that the tangential-1 (vortex) strength is present
//   ps[1] means that the tangential-2 (vortex) strength is present
//   ps[2] means that the normal (source) strength is present
template <class S> using Strength = std::array<std::optional<Vector<S>>,3>;

// useful structure for basis vectors
// Basis<S> b;
//   b[0] is the pair of arrays of the normalized tangential-1 vectors, b[0][0] for x, b[0][1] for y
//   b[2] is the same for the normal vectors (+ is into fluid), b[2][0] for x, b[2][1] for y, ...
template <class S> using Basis = std::array<std::array<Vector<S>,Dimensions>,Dimensions>;


// 2-D elements
template <class S>
class Surfaces: public ElementBase<S> {
public:
  // constructor - accepts vector of (x,y,z) triples, vector of indicies, vector of BCs
  //               and makes one closed body
  //               last parameter (_val) is either fixed strength or boundary condition
  Surfaces(const std::vector<S>&   _x,
           const std::vector<Int>& _idx,
           const std::vector<S>&   _val,
           const elem_t _e,
           const move_t _m,
           std::shared_ptr<Body> _bp)
    : ElementBase<S>(0, _e, _m, _bp),
      np(0),
      source_str_is_unknown(true),
      vol(-1.0),
      //solved_omega(0.0),
      max_strength(-1.0) {

    // make sure input arrays are correctly-sized
    assert(_idx.size() % Dimensions == 0 && "Index array is not an even multiple of dimensions");
    const size_t nsurfs = _idx.size() / Dimensions;

    // always initialize the ps panel strength optionals
    if (this->E != inert) {
      //for (size_t d=0; d<num_unknowns_per_panel(); ++d) {
      for (size_t d=0; d<3; ++d) {
        Vector<S> new_s;
        ps[d] = std::move(new_s);
       }
    }

    // if no surfs, quit out now
    if (nsurfs == 0) {
      // but still initialize ux before we go (in case first bfeature is not enabled)
      if (_bp) this->ux = this->x;
      return;
    }

    assert(_val.size() % nsurfs == 0 && "Value array is not an even multiple of panel count");
    assert(_x.size() % Dimensions == 0 && "Position array is not an even multiple of dimensions");
    const size_t nnodes = _x.size() / Dimensions;

    std::cout << "  new collection with " << nsurfs << " panels and " << nnodes << " nodes" << std::endl;

    // pull out the node locations, they go in the base class
    for (size_t d=0; d<Dimensions; ++d) {
      this->x[d].resize(nnodes);
      for (size_t i=0; i<nnodes; ++i) {
        this->x[d][i] = _x[Dimensions*i+d];
      }
    }

    // save them as untransformed if we are given a Body pointer
    if (_bp) {
      this->ux = this->x;
    }

    // copy over the node indices (with a possible type change)
    bool idx_are_all_good = true;
    idx.resize(_idx.size());
    for (size_t i=0; i<3*nsurfs; ++i) {
      // make sure it exists in the nodes array
      if (_idx[i] >= nnodes) idx_are_all_good = false;
      idx[i] = _idx[i];
    }
    assert(idx_are_all_good && "Some indicies are bad");

    // compute all basis vectors and panel areas
    compute_bases(nsurfs);

    // are strengths/values on a per-node or per-panel basis? - per panel now

    // now, depending on the element type, put the value somewhere
    if (this->E == active) {
      // value is a fixed strength for the panel
      assert(_val.size() == 3*nsurfs && "Value array is not an even multiple of panel count");
      // value is a fixed strength for the panel: x1 and x2 vortex sheet strengths and source sheet str
      for (size_t d=0; d<3; ++d) {
        if (ps[d]) {
          ps[d]->resize(nsurfs);
          for (size_t i=0; i<nsurfs; ++i) {
            (*ps[d])[i] = _val[3*i+d];
          }
        }
      }

      // we still need total strengths (ts)
      vortex_sheet_to_panel_strength(nsurfs);

    } else if (this->E == reactive) {
      // value is a boundary condition
      // always have some value
      for (size_t d=0; d<3; ++d) {
        Vector<S> new_bc(nsurfs);
        bc[d] = std::move(new_bc);
        std::fill(bc[d]->begin(), bc[d]->end(), 0.0);
      }
      // but only set some
      const size_t nper = _val.size()/nsurfs;
      //assert(nper>0 and nper<4 && "Number of boundary conditions is not 1..3");
      assert(nper==3 && "Number of boundary condition values per panel is not 3");
      assert(nper*nsurfs == _val.size() && "Number of bcs is not a factor of nsurfs");
      for (size_t d=0; d<nper; ++d) {
        for (size_t i=0; i<nsurfs; ++i) {
          (*bc[d])[i] = _val[nper*i+d];
        }
      }

      // convert the bc from velocities to basis-vector quantities
      bcs_to_bcs(0,nsurfs);

      // make space for the unknown sheet strengths
      for (size_t d=0; d<3; ++d) {
        if (ps[d]) {
          ps[d]->resize(nsurfs);
          std::fill(ps[d]->begin(), ps[d]->end(), 0.0);
        }
      }

      // we still need total strengths (ts)
      vortex_sheet_to_panel_strength(nsurfs);

    } else if (this->E == inert) {
      // value is ignored (probably zero)
    }

    // velocity is per node, in the base class - just resize it here
    for (size_t d=0; d<Dimensions; ++d) {
      this->u[d].resize(nnodes);
    }

    // but panel velocity is per panel
    for (size_t d=0; d<Dimensions; ++d) {
      pu[d].resize(nsurfs);
    }

    // debug print
    if (false) {
      std::cout << "Nodes" << std::endl;
      for (size_t i=0; i<nnodes; ++i) {
        std::cout << "  " << i << " " << this->x[0][i] << " " << this->x[1][i] << " " << this->x[2][i] << std::endl;
      }
      std::cout << "Triangles" << std::endl;
      for (size_t i=0; i<nsurfs; ++i) {
        std::cout << "  " << i << " " << idx[3*i] << " " << idx[3*i+1] << " " << idx[3*i+2] << std::endl;
      }
    }

    // need to reset the base class n and the local np
    this->n = nnodes;
    np = nsurfs;

    // find geometric center
    if (this->M == bodybound) {
      set_geom_center();
    }
  }

  size_t                         get_npanels()     const { return np; }
  const S                        get_vol()         const { return vol; }
  const std::array<S,Dimensions> get_geom_center() const { return tc; }

  // panel geometry
  const std::vector<Int>&                 get_idx()  const { return idx; }
  const std::array<Vector<S>,Dimensions>& get_x1()   const { return b[0]; }
  const std::array<Vector<S>,Dimensions>& get_x2()   const { return b[1]; }
  const std::array<Vector<S>,Dimensions>& get_norm() const { return b[2]; }
  const Vector<S>&                        get_area() const { return area; }

  // override the ElementBase versions and send the panel-center vels and strengths
  const std::array<Vector<S>,Dimensions>& get_vel() const { return pu; }
  std::array<Vector<S>,Dimensions>&       get_vel()       { return pu; }

  // fixed or unknown surface strengths, or those due to rotation
  const std::array<Vector<S>,Dimensions>& get_str() const { return ts; }
  std::array<Vector<S>,Dimensions>&       get_str()       { return ts; }
  const Vector<S>&                  get_vort1_str() const { return *ps[0]; }
  Vector<S>&                        get_vort1_str()       { return *ps[0]; }
  const Vector<S>&                  get_vort2_str() const { return *ps[1]; }
  Vector<S>&                        get_vort2_str()       { return *ps[1]; }
  const bool                         have_src_str() const { return (bool)ps[2]; }
  const Vector<S>&                    get_src_str() const { return *ps[2]; }
  Vector<S>&                          get_src_str()       { return *ps[2]; }

  // and (reactive only) boundary conditions
  const Vector<S>&                  get_tang1_bcs() const { return *bc[0]; }
  const Vector<S>&                  get_tang2_bcs() const { return *bc[1]; }
  const Vector<S>&                   get_norm_bcs() const { return *bc[2]; }

  // find out the next row index in the BEM after this collection
  void set_first_row(const Int _i) { istart = _i; }
  const bool src_is_unknown() const { return source_str_is_unknown; }
  const Int num_unknowns_per_panel() const { return (source_str_is_unknown ? 3 : 2); }
  const Int get_first_row() const { return istart; }
  const Int get_num_rows()  const { return (get_npanels()*num_unknowns_per_panel() + (is_augmented() ? 3 : 0)); }
  const Int get_next_row()  const { return istart+get_num_rows(); }

  // assign the new strengths from BEM - do not let base class do this
  void set_str(const size_t ioffset, const size_t icnt, Vector<S> _in) {

    assert(ps[0] && "Strength array does not exist");

    // pop off the "unknown" rotation rate and save it
    if (is_augmented()) {
      assert(false && "Augmentation not supported!");
      //solved_omega = _in.back();
      //std::cout << "    solved rotation rate is " << solved_omega << std::endl;
      //omega_error = solved_omega - this->B->get_rotvel();
      //std::cout << "    error in rotation rate is " << omega_error << std::endl;
      //_in.pop_back();
    }

    assert(_in.size() == (*ps[0]).size()*num_unknowns_per_panel() && "Set strength array size does not match");
    //assert(ioffset == 0 && "Offset is not zero");

    // copy the BEM-solved strengths into the panel-strength data structures
    if (source_str_is_unknown and have_src_str()) {
      assert(_in.size() == get_npanels()*3 && "Set strength array size does not match");
      // vortex and source terms
      for (size_t i=0; i<get_npanels(); ++i) {
        (*ps[0])[i] = _in[3*i+0];
        (*ps[1])[i] = _in[3*i+1];
        (*ps[2])[i] = _in[3*i+2];
        //std::cout << "elem " << i << " with area " << area[i] << " has vs " << vs[0][i] << " " << vs[1][i] << std::endl;
      }
    } else {
      assert(_in.size() == get_npanels()*2 && "Set strength array size does not match");
      // vortex only
      for (size_t i=0; i<get_npanels(); ++i) {
        (*ps[0])[i] = _in[2*i+0];
        (*ps[1])[i] = _in[2*i+1];
        //std::cout << "elem " << i << " with area " << area[i] << " has vs " << vs[0][i] << " " << vs[1][i] << std::endl;
      }
    }

    // now recompute the total panel strengths (ts)
    vortex_sheet_to_panel_strength(get_npanels());
  }

  // convert vortex sheet strength to absolute panel strength
  // HACK - throw away the source strength
  void vortex_sheet_to_panel_strength(const size_t num) {

    assert(ps[0]->size() == num && "Input array sizes do not match");
    assert(ps[0]->size() == b[0][0].size() && "Input array sizes do not match");

    // make sure the arrays are sized properly
    for (size_t d=0; d<3; ++d) {
      ts[d].resize(num);
    }

    // convenience references
    std::array<Vector<S>,3>& x1 = b[0];
    std::array<Vector<S>,3>& x2 = b[1];

    // now copy the values over (do them all)
    for (size_t i=0; i<num; ++i) {
      for (size_t d=0; d<Dimensions; ++d) {
        ts[d][i] = ((*ps[0])[i]*x1[d][i] + (*ps[1])[i]*x2[d][i]) * area[i];
      }
      //std::cout << "elem " << i << " has" << std::endl;
      //std::cout << "  x1 " << x1[0][i] << " " << x1[1][i] << " " << x1[2][i] << std::endl;
      //std::cout << "  x2 " << x2[0][i] << " " << x2[1][i] << " " << x2[2][i] << std::endl;
      //std::cout << "  ps " << (*ps[0])[i] << " " << (*ps[1])[i] << std::endl;
      //std::cout << "  ts " << ts[0][i] << " " << ts[1][i] << " " << ts[2][i] << std::endl;
      //std::cout << "elem " << i << " has s " << ts[0][i] << " " << ts[1][i] << " " << ts[2][i] << std::endl;
    }
  }

  // a little logic to see if we should augment the BEM equations for this object
  const bool is_augmented() const {
    bool augment = true;

    // don't augment the ground body, or the boundary to an internal flow
    if (this->B) {
      // is the body pointer ground?
      if (std::string("ground").compare(this->B->get_name()) == 0) {
        // now, does the object bound internal flow?
        if (vol < 0.0) augment = false;
      }
    } else {
      // nullptr for Body? no augment (old way of turning it off)
      augment = false;
    }
    // and only need to augment reactive surfaces (ones participating in BEM)
    if (this->E != reactive) augment = false;

    // force no augmentation at all
    augment = false;

    return augment;
  }

  const float get_max_bc_value() const {
    if (this->E == reactive) {
      S max_bc = 0.0;
      for (size_t i=0; i<get_npanels(); ++i) {
        const S this_bc = std::pow((*bc[0])[i],2) + std::pow((*bc[1])[i],2) + std::pow((*bc[2])[i],2);
        if (this_bc > max_bc) max_bc = this_bc;
      }
      max_bc = std::sqrt(max_bc);
      return (float)max_bc;
    } else {
      return (float)0.0;
    }
  }

  // add more nodes and panels to this collection
  void add_new(const std::vector<S>&   _x,
               const std::vector<Int>& _idx,
               const std::vector<S>&   _val) {

    // remember old sizes of nodes and element arrays
    const size_t nnold = this->n;
    const size_t neold = get_npanels();

    // make sure input arrays are correctly-sized
    assert(_idx.size() % Dimensions == 0 && "Index array is not an even multiple of dimensions");
    const size_t nsurfs = _idx.size() / Dimensions;
    // if no surfs, quit out now
    if (nsurfs == 0) return;

    assert(_val.size() % nsurfs == 0 && "Value array is not an even multiple of panel count");
    assert(_x.size() % Dimensions == 0 && "Position array is not an even multiple of dimensions");
    const size_t nnodes = _x.size() / Dimensions;

    std::cout << "  adding " << nsurfs << " new surface panels and " << nnodes << " new points to collection..." << std::endl;

    // DON'T call the method in the base class, because we do things differently here
    //ElementBase<S>::add_new(_in);

    // pull out the node locations, they are base class
    for (size_t d=0; d<Dimensions; ++d) {
      this->x[d].resize(nnold+nnodes);
      for (size_t i=0; i<nnodes; ++i) {
        this->x[d][nnold+i] = _x[Dimensions*i+d];
      }
    }

    // save them as untransformed if we have a Body pointer
    if (this->B) {
      for (size_t d=0; d<Dimensions; ++d) {
        (*this->ux)[d].resize(nnold+nnodes);
        for (size_t i=nnold; i<nnold+nnodes; ++i) {
          (*this->ux)[d][i] = this->x[d][i];
        }
      }
    }

    // copy over the node indices, taking care to offset into the new array
    bool idx_are_all_good = true;
    idx.resize(3*neold + _idx.size());
    for (size_t i=0; i<3*nsurfs; ++i) {
      // make sure it exists in the nodes array
      if (_idx[i] >= nnold+nnodes) idx_are_all_good = false;
      idx[3*neold+i] = nnold + _idx[i];
    }
    assert(idx_are_all_good && "Some indicies are bad");

    // compute all basis vectors and panel areas
    compute_bases(neold+nsurfs);

    // now, depending on the element type, put the value somewhere
    if (this->E == active) {
      // value is a fixed strength for the element
      assert(_val.size() == 3*nsurfs && "Value array is not an even multiple of panel count");
      // value is a fixed strength for the panel: x1 and x2 vortex sheet strengths and source sheet str
      for (size_t d=0; d<3; ++d) {
        if (ps[d]) {
          ps[d]->resize(neold+nsurfs);
          for (size_t i=0; i<nsurfs; ++i) {
            (*ps[d])[neold+i] = _val[3*i+d];
          }
        }
      }

      // and ensure that the raw strengths are resized and set
      vortex_sheet_to_panel_strength(neold+nsurfs);

    } else if (this->E == reactive) {
      // value is a boundary condition
      // always resize arrays
      for (size_t d=0; d<3; ++d) {
        bc[d]->resize(neold+nsurfs);
      }
      const size_t nper = _val.size()/nsurfs;
      //assert(nper>0 and nper<4 && "Number of boundary conditions is not 1..3");
      assert(nper==3 && "Number of boundary condition values is not 3");
      for (size_t d=0; d<nper; ++d) {
        for (size_t i=0; i<nsurfs; ++i) {
          (*bc[d])[neold+i] = _val[nper*i+d];
        }
      }

      // convert the bc from velocities to basis-vector quantities
      bcs_to_bcs(neold,nsurfs);

      // upsize vortex sheet and raw strength arrays, too
      for (size_t d=0; d<3; ++d) {
        if (ps[d]) {
          ps[d]->resize(neold+nsurfs);
          //std::fill(ps[d]->begin(), ps[d]->end(), 0.0);
        }
      }
      for (size_t d=0; d<3; ++d) {
        ts[d].resize(neold+nsurfs);
      }

      // and ensure that the raw strengths are resized and set
      vortex_sheet_to_panel_strength(neold+nsurfs);

    } else if (this->E == inert) {
      // value is ignored (probably zero)
    }

    // velocity is in the base class - just resize it here
    for (size_t d=0; d<Dimensions; ++d) {
      this->u[d].resize(nnold+nnodes);
    }

    // panel velocity is here
    for (size_t d=0; d<Dimensions; ++d) {
      pu[d].resize(neold+nsurfs);
    }

    // debug print
    if (false) {
      std::cout << "Nodes" << std::endl;
      for (size_t i=0; i<nnold+nnodes; ++i) {
        std::cout << "  " << i << " " << this->x[0][i] << " " << this->x[1][i] << " " << this->x[2][i] << std::endl;
      }
      std::cout << "Tris" << std::endl;
      for (size_t i=0; i<neold+nsurfs; ++i) {
        std::cout << "  " << i << " " << idx[3*i] << " " << idx[3*i+1] << " " << idx[3*i+2] << std::endl;
      }
    }

    // need to reset the base class n
    this->n += nnodes;
    np += nsurfs;

    // re-find geometric center
    if (this->M == bodybound) {
      set_geom_center();
    }
  }

  void add_body_motion(const S _factor, const double _time) {
    // no need to call base class now
    //ElementBase<S>::add_body_motion(_factor);

    // if no body pointer, or attached to ground, return
    if (not this->B) return;
    if (std::string("ground").compare(this->B->get_name()) == 0) return;

    // make sure we've calculated transformed center (we do this when we do volume)
    //assert(vol > 0.0 && "Have not calculated transformed center, or volume is negative");
    // and we trust that we've transformed utc to tc

    // do this for all nodes - what about panels?
    for (size_t i=0; i<get_npanels(); ++i) {

      // apply the translational velocity
      const std::array<double,Dimensions> thisvel = this->B->get_vel(_time);
      for (size_t d=0; d<Dimensions; ++d) {
        pu[d][i] += _factor * (float)thisvel[d];
      }

      // now compute the rotational velocity with respect to the geometric center
      const Vec rotvel = this->B->get_rotvel_vec(_time);
      // indices for this panel
      const Int id0 = idx[3*i];
      const Int id1 = idx[3*i+1];
      const Int id2 = idx[3*i+2];
      // panel center minus body center
      const S xc = (this->x[0][id0] + this->x[0][id1] + this->x[0][id2])/3.0 - tc[0];
      const S yc = (this->x[1][id0] + this->x[1][id1] + this->x[1][id2])/3.0 - tc[1];
      const S zc = (this->x[2][id0] + this->x[2][id1] + this->x[2][id2])/3.0 - tc[2];
      // add rotational velocity
      pu[0][i] += _factor * (rotvel[1]*zc - rotvel[2]*yc);
      pu[1][i] += _factor * (rotvel[2]*xc - rotvel[0]*zc);
      pu[2][i] += _factor * (rotvel[0]*yc - rotvel[1]*xc);
    }

    //std::cout << "in add_body_motion, pu is" << std::endl;
    //for (size_t i=0; i<100; ++i) {
    //  std::cout << i << " " << pu[0][i] << " " << pu[1][i] << " " << pu[2][i] << " " << std::endl;
    //}
  }
 
  void zero_strengths() {
    // call base class first
    ElementBase<S>::zero_strengths();

    // and reset any panel vortex or source strengths
    for (size_t d=0; d<3; ++d) {
      if (ps[d]) {
        std::fill(ps[d]->begin(), ps[d]->end(), 0.0);
      }
    }
  }

  // three ways to add source and vortex rotational strengths to the surface
  // first: as a multiple of the current, defined rotation rate
  void add_rot_strengths(const S _factor) {
    // if no parent Body, forget it
    if (not this->B) return;
    // get current rotation rate
    std::array<double,Dimensions> rotvel = this->B->get_rotvel_vec();
    for (size_t d=0; d<3; ++d) rotvel[d] *= _factor;
    // call parent
    add_rot_strengths_base(rotvel);
  }

  // second: assuming unit rotation rate (for finding the BEM influence matrix)
  void add_unit_rot_strengths(const int _d) {
    assert((_d>0 and _d<3) && "Invalid dimension");
    std::array<double,Dimensions> rotvel = {{0.0, 0.0, 0.0}};
    rotvel[_d] = 1.0;
    add_rot_strengths_base(rotvel);
  }

  // third: as a multiple of the rotation rate
  void add_solved_rot_strengths(const S _factor) {
    if (is_augmented()) {
      // use the augmented-BEM result for rotation rate
      assert(false && "Augmentation not yet supported");
      //add_rot_strengths_base(_factor * solved_omega);
    } else {
      // use the predefined rotationrate
      add_rot_strengths(_factor);
    }
  }

  // augment the strengths with a value equal to that which accounts for
  //   the solid-body rotation of the object
  // NOTE: this needs to provide both the vortex AND source strengths!
  // AND we don't have the time - assume bodies have been transformed
  void add_rot_strengths_base(const std::array<double,Dimensions> _rotvel) {

    // No - ALWAYS allow this to run
    // if no rotation, or no parent Body, or attached to ground, then no problem!
    //if (not this->B) return;
    //if (std::string("ground").compare(this->B->get_name()) == 0) return;

    // if no rotation, then we don't need to add anything
    //const auto rotvel = this->B->get_rotvel_vec();
    //if (std::abs(rotvel[0]) + std::abs(rotvel[1]) + std::abs(rotvel[2])
    //    > std::numeric_limits<double>::epsilon()) return;

    // make sure we've calculated transformed center (we do this when we do volume)
    //assert(vol > 0.0 && "Have not calculated transformed center, or volume is negative");
    // and we trust that we've transformed utc to tc

    // have we made source strength vector yet? or is it the right size?
    if (ps[2]) {
      ps[2]->resize(get_npanels());
    } else {
      // value is a fixed strength for the segment
      Vector<S> new_ss(get_npanels());
      *ps[2] = std::move(new_ss);
    }

    //std::cout << "Inside add_rot_strengths, sizes are: " << get_npanels() << " " << ps[2]->size() << std::endl;
    assert(ts[0].size() == get_npanels() && "Strength array is not the same as panel count");

    // get basis vectors
    std::array<Vector<S>,3>& t1   = b[0];
    std::array<Vector<S>,3>& t2   = b[1];
    std::array<Vector<S>,3>& norm = b[2];

    // still here? let's do it. use the untransformed coordinates
    for (size_t i=0; i<get_npanels(); i++) {
      // indices for this panel
      const Int id0 = idx[3*i];
      const Int id1 = idx[3*i+1];
      const Int id2 = idx[3*i+2];
      // vector from object geometric center to panel center
      const S dx = ((*this->ux)[0][id0] + (*this->ux)[0][id1] + (*this->ux)[0][id2])/3.0 - utc[0];
      const S dy = ((*this->ux)[1][id0] + (*this->ux)[1][id1] + (*this->ux)[1][id2])/3.0 - utc[1];
      const S dz = ((*this->ux)[2][id0] + (*this->ux)[2][id1] + (*this->ux)[2][id2])/3.0 - utc[2];
      // velocity of the panel center
      const S ui = _rotvel[1]*dz - _rotvel[2]*dy;
      const S vi = _rotvel[2]*dx - _rotvel[0]*dz;
      const S wi = _rotvel[0]*dy - _rotvel[1]*dx;

      // the vortex strengths - we ADD to the existing
      const S new_vort1 = -1.0 * (ui*t1[0][i] + vi*t1[1][i] + wi*t1[2][i]);
      (*ps[0])[i] += new_vort1;
      const S new_vort2 = -1.0 * (ui*t2[0][i] + vi*t2[1][i] + wi*t2[2][i]);
      (*ps[1])[i] += new_vort2;

      // the source strength
      const S new_src = -1.0 * (ui*norm[0][i] + vi*norm[1][i] + wi*norm[2][i]);
      (*ps[2])[i] += new_src;

      // debug print
      if (std::abs(_rotvel[0]) > 0.0 and false) {
        std::cout << "  panel " << i << " at " << dx << " " << dy << " " << dz << " adds to vortex str "
                  << new_vort1 << " " << new_vort2 << " and source str " << new_src << std::endl;
      }
    }
  }

  // calculate the geometric center of all geometry in this object
  void set_geom_center() {

    // we must have an attached body and a set of untransformed coordinates
    assert(this->B && "Body pointer has not been set");
    assert(this->ux && "Untransformed positions have not been set");

    std::cout << "  inside Surfaces::set_geom_center with " << get_npanels() << " panels" << std::endl;

    // iterate over panels, accumulating vol and CM
    double vsum = 0.0;
    double xsum = 0.0;
    double ysum = 0.0;
    double zsum = 0.0;
    for (size_t i=0; i<get_npanels(); i++) {
      const size_t jp0 = idx[3*i];
      const size_t jp1 = idx[3*i+1];
      const size_t jp2 = idx[3*i+2];
      // assume a triangle from 0,0 to two ends of each panel
      double thisvol = (*this->ux)[0][jp0] * (double)(*this->ux)[1][jp1] * (*this->ux)[2][jp2]
                     - (*this->ux)[0][jp0] * (double)(*this->ux)[1][jp2] * (*this->ux)[2][jp1]
                     - (*this->ux)[0][jp1] * (double)(*this->ux)[1][jp0] * (*this->ux)[2][jp2]
                     + (*this->ux)[0][jp1] * (double)(*this->ux)[1][jp2] * (*this->ux)[2][jp0]
                     + (*this->ux)[0][jp2] * (double)(*this->ux)[1][jp0] * (*this->ux)[2][jp1]
                     - (*this->ux)[0][jp2] * (double)(*this->ux)[1][jp1] * (*this->ux)[2][jp0];
      thisvol /= 6.0;
      // add this to the running sums
      //std::cout << "    and area " << thisvol << " and center " << xc << " " << yc << std::endl;
      vsum += thisvol;
      xsum += 0.25 * thisvol * ((*this->ux)[0][jp0] + (*this->ux)[0][jp1] + (*this->ux)[0][jp2]);
      ysum += 0.25 * thisvol * ((*this->ux)[1][jp0] + (*this->ux)[1][jp1] + (*this->ux)[1][jp2]);
      zsum += 0.25 * thisvol * ((*this->ux)[2][jp0] + (*this->ux)[2][jp1] + (*this->ux)[2][jp2]);
    }
    vol = (S)vsum;
    utc[0] = (S)xsum / vol;
    utc[1] = (S)ysum / vol;
    utc[2] = (S)zsum / vol;

    std::cout << "    geom center is " << utc[0] << " " << utc[1] << " " << utc[2] << " and vol is " << vol << std::endl;
  }

  // need to maintain the 3x3 set of basis vectors for each panel
  // this also calculates the triangle areas
  // always recalculate everything!
  void compute_bases(const Int nnew) {

    assert(3*nnew == idx.size() && "Array size mismatch");

    // resize any vectors
    for (size_t i=0; i<Dimensions; ++i) {
      for (size_t j=0; j<Dimensions; ++j) {
        b[i][j].resize(nnew);
      }
    }
    area.resize(nnew);

    // we'll reuse these vectors
    std::array<S,Dimensions> x1, x2, norm;

    // update everything
    for (size_t i=0; i<nnew; ++i) {
      const size_t id0 = idx[3*i];
      const size_t id1 = idx[3*i+1];
      const size_t id2 = idx[3*i+2];
      //std::cout << "elem near " << this->x[0][id0] << " " << this->x[1][id0] << " " << this->x[2][id0] << std::endl;

      // x1 vector is along direction from node 0 to node 1
      for (size_t j=0; j<3; ++j) x1[j] = this->x[j][id1] - this->x[j][id0];
      const S base = length(x1);
      for (size_t j=0; j<3; ++j) x1[j] *= (1.0/base);
      //std::cout << "  has x1 " << x1[0] << " " << x1[1] << " " << x1[2] << std::endl;

      // x2 vector is perpendicular to x1, and points toward node 2
      for (size_t j=0; j<3; ++j) x2[j] = this->x[j][id2] - this->x[j][id0];
      const S dp = dot_product<S>(x2, x1);
      for (size_t j=0; j<3; ++j) x2[j] -= dp*x1[j];
      const S height = length(x2);
      for (size_t j=0; j<3; ++j) x2[j] *= (1.0/height);
      //std::cout << "  has x2 " << x2[0] << " " << x2[1] << " " << x2[2] << std::endl;

      // now we have the area
      area[i] = 0.5 * base * height;

      // x3 is the normal vector, pointing into the fluid
      cross_product(x1, x2, norm);
      //std::cout << "  norm " << norm[0] << " " << norm[1] << " " << norm[2] << std::endl;

      // and assign
      for (size_t j=0; j<Dimensions; ++j) b[0][j][i] = x1[j];
      for (size_t j=0; j<Dimensions; ++j) b[1][j][i] = x2[j];
      for (size_t j=0; j<Dimensions; ++j) b[2][j][i] = norm[j];

      //std::cout << "elem near " << this->x[0][id0] << " " << this->x[1][id0] << " " << this->x[2][id0] << " has norm " << b[2][0][i] << " " << b[2][1][i] << " " << b[2][2][i] << std::endl;
    }
  }

  // need to convert the boundary conditions from global velocities into basis-frame vels
  void bcs_to_bcs(const Int nold, const Int nnew) {

    // get basis vectors - we just computed them
    std::array<Vector<S>,3>& t1 = b[0];
    std::array<Vector<S>,3>& t2 = b[1];
    std::array<Vector<S>,3>& nm = b[2];

    // update only the new panels
    for (size_t i=nold; i<nold+nnew; ++i) {
      const std::array<S,3> oldbc = {{ (*bc[0])[i], (*bc[1])[i], (*bc[2])[i] }};
      // normal velocity component is easy
      (*bc[2])[i] = nm[0][i]*oldbc[0] + nm[0][i]*oldbc[1] + nm[0][i]*oldbc[2];
      // but for tangential we need to recognize that t1 is the direction of the circulation,
      //   so the corresponding velocity is perpendicular to it!
      (*bc[0])[i] = +(t2[0][i]*oldbc[0] + t2[0][i]*oldbc[1] + t2[0][i]*oldbc[2]);
      (*bc[1])[i] = -(t1[0][i]*oldbc[0] + t1[0][i]*oldbc[1] + t1[0][i]*oldbc[2]);

      //std::cout << "elem near " << this->x[0][id0] << " " << this->x[1][id0] << " " << this->x[2][id0] << " has norm " << b[2][0][i] << " " << b[2][1][i] << " " << b[2][2][i] << std::endl;
    }
  }

  // when transforming a body-bound object to a new time, we must also transform the geometric center
  void transform(const double _time) {
    // must explicitly call the method in the base class
    ElementBase<S>::transform(_time);

    // and recalculate the basis vectors
    compute_bases(np);

    if (this->B and this->M == bodybound) {
      // do the transform with an affine matrix
      Eigen::Transform<double,3,Eigen::Affine> xform = this->B->get_transform_mat();
      const Eigen::Vector3d _pre = {utc[0], utc[1], utc[2]};
      const Eigen::Vector3d _post = xform * _pre;
      tc[0] = _post(0);
      tc[1] = _post(1);
      tc[2] = _post(2);

      // HACK - might need to compute_bases() from here - but only if rotation happened

    } else {
      // copy utc to tc
      tc[0] = utc[0];
      tc[1] = utc[1];
      tc[2] = utc[2];
    }
  }


  void zero_vels() {
    // zero the local, panel-center vels
    for (size_t d=0; d<Dimensions; ++d) {
      std::fill(pu[d].begin(), pu[d].end(), 0.0);
    }
    // then explicitly call the method in the base class to zero theirs
    ElementBase<S>::zero_vels();
  }

  void finalize_vels(const std::array<double,Dimensions>& _fs) {
    // finalize panel-center vels first
    const double factor = 0.25/M_PI;
    for (size_t d=0; d<Dimensions; ++d) {
      for (size_t i=0; i<get_npanels(); ++i) {
        pu[d][i] = _fs[d] + pu[d][i] * factor;
      }
    }
    // must explicitly call the method in the base class, too
    ElementBase<S>::finalize_vels(_fs);
  }

/*
  // up-size all arrays to the new size, filling with sane values
  void resize(const size_t _nnew) {
    const size_t currn = this->n;
    //std::cout << "  inside Surfaces::resize with " << currn << " " << _nnew << std::endl;

    // must explicitly call the method in the base class - this sets n
    ElementBase<S>::resize(_nnew);

    if (_nnew == currn) return;

    // radii here
    const size_t thisn = r.size();
    r.resize(_nnew);
    for (size_t i=thisn; i<_nnew; ++i) {
      r[i] = 1.0;
    }
  }

  //
  // 1st order Euler advection and stretch
  //
  void move(const double _time, const double _dt) {
    // must explicitly call the method in the base class
    ElementBase<S>::move(_time, _dt);

    // no specialization needed
    if (this->M == lagrangian and this->E != inert) {
      //std::cout << "  Stretching" << to_string() << " using 1st order" << std::endl;
      S thismax = 0.0;

      for (size_t i=0; i<this->n; ++i) {
        S this_s = (*this->s)[i];

        // compute stretch term
        std::array<S,2> wdu = {0.0};

        // add Cottet SFS

        // update strengths
        (*this->s)[i] = this_s + _dt * wdu[0];

        // check for max strength
        S thisstr = std::abs((*this->s)[i]);
        if (thisstr > thismax) thismax = thisstr;

      }
      if (max_strength < 0.0) {
        max_strength = thismax;
      } else {
        max_strength = 0.1*thismax + 0.9*max_strength;
      }
      //std::cout << "  New max_strength is " << max_strength << std::endl;
    } else {
      //std::cout << "  Not stretching" << to_string() << std::endl;
      max_strength = 1.0;
    }
  }
*/

  //
  // return a particle version of the panels (useful during Diffusion)
  // offset is in world units - NOT scaled
  //
  std::vector<S> represent_as_particles(const S _offset, const S _vdelta) {

    // how many panels?
    const size_t num_pts = get_npanels();

    // recompute the total strengths (ts)
    vortex_sheet_to_panel_strength(num_pts);

    // init the output vector (x, y, z, sx, sy, sz, r)
    std::vector<S> px(num_pts*7);

    // get basis vectors
    std::array<Vector<S>,3>&   x1 = b[0];
    std::array<Vector<S>,3>&   x2 = b[1];
    std::array<Vector<S>,3>& norm = b[2];
    Vector<S>&                bc1 = *bc[0];
    Vector<S>&                bc2 = *bc[1];

    for (size_t i=0; i<num_pts; i++) {
      const Int id0 = idx[3*i];
      const Int id1 = idx[3*i+1];
      const Int id2 = idx[3*i+2];
      const Int idx = 7*i;
      // start at center of panel
      for (size_t j=0; j<3; ++j) px[idx+j] = (1./3.) * (this->x[j][id0] + this->x[j][id1] + this->x[j][id2]);
      // push out a fixed distance
      // this assumes properly resolved, vdelta and dt
      for (size_t j=0; j<3; ++j) px[idx+j] += _offset * norm[j][i];
      // the panel strength is the solved strength plus the boundary condition
      for (size_t j=0; j<3; ++j) px[idx+3+j] = ts[j][i];
      // add on the (vortex) bc values here
      if (this->E == reactive) {
        for (size_t j=0; j<3; ++j) px[idx+3+j] += (bc1[i]*x1[j][i] + bc2[i]*x2[j][i]) * area[i];
      }
      // IGNORE SOURCE SHEET STRENGTHS
      // and the core size
      px[idx+6] = _vdelta;
      //std::cout << "  new part at " << px[idx+0] << " " << px[idx+1] << " " << px[idx+2];
      //std::cout << "     with str " << px[idx+3] << " " << px[idx+4] << " " << px[idx+5] << std::endl;
    }

    return px;
  }

  // find the new peak vortex sheet strength magnitude
  S get_max_str() {
    if (this->E != inert) {
      S max_ps = 0.0;
      for (size_t i=0; i<get_npanels(); ++i) {
        const S this_ps = std::pow((*ps[0])[i],2) + std::pow((*ps[1])[i],2);
        if (this_ps > max_ps) max_ps = this_ps;
      }
      max_ps = std::sqrt(max_ps);
      return (S)max_ps;
    } else {
      return (S)1.0;
    }
  }

  // smooth the peak strength magnitude
  void update_max_str() {
    S thismax = get_max_str();

    // and slowly update the current saved value
    if (max_strength < 0.0) {
      max_strength = thismax;
    } else {
      max_strength = 0.1*thismax + 0.9*max_strength;
    }
  }

  // add and return the total circulation of all elements
  //   specializing the one in ElementBase because we need
  //   to scale by panel area here
  std::array<S,Dimensions> get_total_circ(const double _time) {

    // here is the return vector
    std::array<S,Dimensions> circ;
    circ.fill(0.0);

    if (this->E != inert) {
      // make this easy - represent as particles
      std::vector<S> pts = represent_as_particles(0.0, 1.0);

      // now compute impulse of those
      for (size_t i=0; i<get_npanels(); ++i) {
        const size_t idx = 7*i;
        circ[0] += pts[idx+3+0];
        circ[1] += pts[idx+3+1];
        circ[2] += pts[idx+3+2];
      }
    }

    return circ;
  }

  // add and return the total circulation of all elements
  std::array<S,3> get_body_circ(const double _time) {
    std::array<S,3> circ;
    circ.fill(0.0);

    // do not call the parent
    if (this->B) {
      // we're attached to a body - great! what's the rotation rate?
      const Vec rotvel = this->B->get_rotvel_vec(_time);
      circ[0] = 2.0 * vol * rotvel[0];
      circ[1] = 2.0 * vol * rotvel[1];
      circ[2] = 2.0 * vol * rotvel[2];
    } else {
      // we are fixed, thus not rotating
    }

    return circ;
  }

  // add and return the total impulse of all elements
  std::array<S,Dimensions> get_total_impulse() {

    // here is the return vector
    std::array<S,Dimensions> imp;
    imp.fill(0.0);

    if (this->E != inert) {
      // make this easy - represent as particles
      std::vector<S> pts = represent_as_particles(0.0, 1.0);

      // now compute impulse of those
      for (size_t i=0; i<get_npanels(); ++i) {
        const size_t idx = 7*i;
        imp[0] += pts[idx+3+1] * pts[idx+2] - pts[idx+3+2] * pts[idx+1];
        imp[1] += pts[idx+3+2] * pts[idx+0] - pts[idx+3+0] * pts[idx+2];
        imp[2] += pts[idx+3+0] * pts[idx+1] - pts[idx+3+1] * pts[idx+0];
      }
    }

    return imp;
  }


#ifdef USE_GL
  //
  // OpenGL functions
  //

  // helper function to clean up initGL
  void prepare_opengl_buffer(GLuint _prog, GLuint _idx, const GLchar* _name) {
    glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[_idx]);
    const GLint position_attribute = glGetAttribLocation(_prog, _name);
    // Specify how the data for position can be accessed
    glVertexAttribPointer(position_attribute, 1, get_gl_type<S>, GL_FALSE, 0, 0);
    // Enable the attribute
    glEnableVertexAttribArray(position_attribute);
  }

  // this gets done once - load the shaders, set up the vao
  void initGL(std::vector<float>& _projmat,
              float*              _poscolor,
              float*              _negcolor,
              float*              _defcolor) {

    //std::cout << "inside Surfaces.initGL" << std::endl;
    std::cout << "inside Surfaces.initGL with E=" << this->E << " and M=" << this->M << std::endl;

    // generate the opengl state object with space for 7 vbos and 1 shader program
    mgl = std::make_shared<GlState>(7,1);

    // Allocate space, but don't upload the data from CPU to GPU yet
    for (size_t i=0; i<Dimensions; ++i) {
      glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i]);
      glBufferData(GL_ARRAY_BUFFER, 0, this->x[i].data(), GL_STATIC_DRAW);
    }

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mgl->vbo[Dimensions]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 0, idx.data(), GL_STATIC_DRAW);

    if (this->E != inert) {
      for (size_t i=0; i<Dimensions; ++i) {
        glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i+4]);
        glBufferData(GL_ARRAY_BUFFER, 0, ts[i].data(), GL_STATIC_DRAW);
      }
    }

    // Load and create the blob-drawing shader program
    mgl->spo[0] = create_draw_surface_tri_prog();

    // Now do the arrays
    prepare_opengl_buffer(mgl->spo[0], 0, "px");
    prepare_opengl_buffer(mgl->spo[0], 1, "py");
    prepare_opengl_buffer(mgl->spo[0], 2, "posz");
    prepare_opengl_buffer(mgl->spo[0], 3, "sx");
    prepare_opengl_buffer(mgl->spo[0], 4, "sy");
    prepare_opengl_buffer(mgl->spo[0], 5, "sz");

    // and for the compute shaders! (later)

    // Get the location of the attributes that enters in the vertex shader
    mgl->projmat_attribute = glGetUniformLocation(mgl->spo[0], "Projection");

    // upload the projection matrix
    glUniformMatrix4fv(mgl->projmat_attribute, 1, GL_FALSE, _projmat.data());

    // locate where the colors and color scales go
    mgl->pos_color_attribute = glGetUniformLocation(mgl->spo[0], "pos_color");
    mgl->neg_color_attribute = glGetUniformLocation(mgl->spo[0], "neg_color");
    mgl->def_color_attribute = glGetUniformLocation(mgl->spo[0], "def_color");
    mgl->str_scale_attribute = glGetUniformLocation(mgl->spo[0], "str_scale");

    // send the current values
    glUniform4fv(mgl->pos_color_attribute, 1, (const GLfloat *)_poscolor);
    glUniform4fv(mgl->neg_color_attribute, 1, (const GLfloat *)_negcolor);
    glUniform4fv(mgl->def_color_attribute, 1, (const GLfloat *)_defcolor);
    glUniform1f (mgl->str_scale_attribute, (const GLfloat)1.0);
    //std::cout << "init pos color as " << _poscolor[0] << " " << _poscolor[1] << " " << _poscolor[2] << " " << _poscolor[3] << std::endl;

    // and indicate the fragment color output
    glBindFragDataLocation(mgl->spo[0], 0, "frag_color");

    glBindVertexArray(0);
  }

  // this gets done every time we change the size of the index array
  void updateGL() {
    //std::cout << "inside Surfaces.updateGL" << std::endl;

    // has this been init'd yet?
    if (not mgl) return;
    if (glIsVertexArray(mgl->vao) == GL_FALSE) return;

    const size_t vlen = this->x[0].size()*sizeof(S);
    if (vlen > 0) {
      glBindVertexArray(mgl->vao);

      // Indicate and upload the data from CPU to GPU
      for (size_t i=0; i<Dimensions; ++i) {
        // the positions
        glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i]);
        glBufferData(GL_ARRAY_BUFFER, vlen, this->x[i].data(), GL_DYNAMIC_DRAW);
      }

      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mgl->vbo[Dimensions]);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(Int)*idx.size(), idx.data(), GL_DYNAMIC_DRAW);

      // here is where we split on element type: active/reactive vs. inert
      if (this->E == inert) {
        // just don't upload strengths

      } else { // this->E is active or reactive
        const size_t slen = ts[0].size()*sizeof(S);
        for (size_t i=0; i<Dimensions; ++i) {
          glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i+4]);
          glBufferData(GL_ARRAY_BUFFER, slen, ts[i].data(), GL_DYNAMIC_DRAW);
        }
      }

      glBindVertexArray(0);

      // must tell draw call how many elements are there - or, really, how many indices
      mgl->num_uploaded = idx.size();
    }
  }

  // OpenGL3 stuff to draw triangles, called once per frame
  void drawGL(std::vector<float>& _modelviewmat,
              std::vector<float>& _projmat,
              RenderParams&       _rparams,
              const float         _vdelta) {

    //std::cout << "inside Surfaces.drawGL" << std::endl;

    // has this been init'd yet?
    if (not mgl) {
      initGL(_projmat, _rparams.pos_circ_color,
                       _rparams.neg_circ_color,
                       _rparams.default_color);
      updateGL();
    }

    if (mgl->num_uploaded > 0) {
      glBindVertexArray(mgl->vao);

      // get blending ready
      glDisable(GL_DEPTH_TEST);
      glEnable(GL_BLEND);
      glBlendFunc(GL_ONE, GL_ONE);

      glLineWidth(2.0);

      // here is where we split on element type: active/reactive vs. inert
      //if (this->E == inert) {
      //} else { // this->E is active or reactive

      // draw as triangles
      glUseProgram(mgl->spo[0]);

      // multiply the two matrices to get the MVP matrix
      Eigen::Matrix<float,4,4> mvp;
      Eigen::Map<Eigen::Matrix<float,4,4>> pmat(_projmat.data());
      Eigen::Map<Eigen::Matrix<float,4,4>> mvmat(_modelviewmat.data());
      mvp = pmat * mvmat;

      // upload the current projection matrix
      glUniformMatrix4fv(mgl->projmat_attribute, 1, GL_FALSE, mvp.data());

      // upload the current color values
      glUniform4fv(mgl->pos_color_attribute, 1, (const GLfloat *)_rparams.pos_circ_color);
      glUniform4fv(mgl->neg_color_attribute, 1, (const GLfloat *)_rparams.neg_circ_color);
      glUniform4fv(mgl->def_color_attribute, 1, (const GLfloat *)_rparams.default_color);
      glUniform1f (mgl->str_scale_attribute, (const GLfloat)max_strength);

      // the one draw call here
      glDrawElements(GL_TRIANGLES, mgl->num_uploaded, get_gl_type<Int>, 0);

      // return state
      glEnable(GL_DEPTH_TEST);
      glDisable(GL_BLEND);
      glBindVertexArray(0);
    }
  }
#endif

  std::string to_string() const {
    std::string retstr = " " + std::to_string(get_npanels()) + ElementBase<S>::to_string() + " Panels";
    return retstr;
  }

protected:
  // ElementBase.h has x, s, u, ux on the *nodes*

  size_t np;				// number of panels

  // element-wise variables special to triangular panels
  std::vector<Int>                 idx;	// indexes into the x array
  Vector<S>                       area; // panel areas
  Basis<S>                           b; // transformed basis vectors: x1 is b[0], x2 is b[1], normal is b[2], normal x is b[2][0]
  std::array<Vector<S>,Dimensions>  pu; // velocities on panel centers - "u" is node vels in ElementBase

  // strengths and BCs
  Strength<S>                       ps; // panel-wise strengths per area (for "active" and "reactive")
                                        // vortex sheet strengths are 0 and 1, source is in 2
  Strength<S>                       bc; // boundary condition for the elements (only when "reactive")
  std::array<Vector<S>,Dimensions>  ts; // total element strengths (do not use s in ElementBase)
  bool           source_str_is_unknown; // should the BEM solve for source strengths?

  // parameters for the encompassing body
  Int                           istart; // index of first entry in RHS vector and A matrix
  S                                vol; // volume of the body - for augmented BEM solution
  std::array<S,Dimensions>         utc; // untransformed geometric center
  std::array<S,Dimensions>          tc; // transformed geometric center

  // augmented-BEM-related (see Omega2D for more)
  //std::array<S,Dimensions> solved_omega; // rotation rate returned from augmented row in BEM

private:
#ifdef USE_GL
  std::shared_ptr<GlState> mgl;		// for drawing only
#endif
  float max_strength;
};

