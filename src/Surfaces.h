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
//#include <memory>
//#include <optional>
//#include <random>
#include <cassert>


// 1-D elements
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
      vol(-1.0),
      max_strength(-1.0) {

    // make sure input arrays are correctly-sized
    assert(_x.size() % Dimensions == 0 && "Position array is not an even multiple of dimensions");
    assert(_idx.size() % Dimensions == 0 && "Index array is not an even multiple of dimensions");
    const size_t nnodes = _x.size() / Dimensions;
    const size_t nsurfs = _idx.size() / Dimensions;
    assert(_val.size() % nsurfs == 0 && "Value array is not an even multiple of panel count");

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
    compute_bases();

    // now, depending on the element type, put the value somewhere
    if (this->E == active) {
      // value is a fixed strength for the segment
      assert(_val.size() == 2*nsurfs);
      // value is a fixed strength for the panel: x1 and x2 vortex sheet strengths
      for (size_t d=0; d<2; ++d) {
        vs[d].resize(nsurfs);
        for (size_t i=0; i<nsurfs; ++i) {
          vs[d][i] = _val[2*i+d];
        }
      }

      // we still need general strengths
      // NOTE that this means that a vector in ElementBase is NOT sized to n, but to nsurfs!
      vortex_sheet_to_panel_strength();

    } else if (this->E == reactive) {
      // value is a boundary condition
      const size_t nper = _val.size()/nsurfs;
      assert(nper>0 and nper<4);
      bc.resize(nper);
      for (size_t d=0; d<nper; ++d) {
        bc[d].resize(nsurfs);
        for (size_t i=0; i<nsurfs; ++i) {
          bc[d][i] = _val[nper*i+d];
        }
      }

      // make space for panel-centric strengths
      for (size_t d=0; d<2; ++d) {
        vs[d].resize(nsurfs);
        std::fill(vs[d].begin(), vs[d].end(), 0.0);
      }

      // we still need general strengths
      vortex_sheet_to_panel_strength();

    } else if (this->E == inert) {
      // value is ignored (probably zero)
    }

    // velocity is per node, in the base class - just resize it here
    for (size_t d=0; d<Dimensions; ++d) {
      this->u[d].resize(nnodes);
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

    // need to reset the base class n
    this->n = nnodes;

    // find geometric center
    if (this->M == bodybound) {
      set_geom_center();
    }
  }

  size_t                         get_npanels()     const { return idx.size()/Dimensions; }
  const S                        get_vol()         const { return vol; }
  const std::array<S,Dimensions> get_geom_center() const { return tc; }

  // callers should never have to change this array
  const std::vector<Int>&        get_idx()         const { return idx; }
  const std::vector<Vector<S>>&  get_bcs()         const { return bc; }

  // override the ElementBase versions and send the panel-center vels
  const std::array<Vector<S>,Dimensions>& get_vel() const { return pu; }
  std::array<Vector<S>,Dimensions>&       get_vel()       { return pu; }

  // vortex and source strengths
  const std::array<Vector<S>,2>&  get_vort_str() const { return vs; }
  const bool                      have_src_str() const { return (bool)ss; }
  const Vector<S>&                get_src_str()  const { return *ss; }
  const std::array<Vector<S>,3>&  get_x1()       const { return b[0]; }
  const std::array<Vector<S>,3>&  get_x2()       const { return b[1]; }
  const std::array<Vector<S>,3>&  get_norm()     const { return b[2]; }
  const Vector<S>&                get_area()     const { return area; }

  // find out the next row index in the BEM after this collection
  void set_first_row(const Int _i) { istart = _i; }
  const Int get_first_row() const { return istart; }
  const Int get_num_rows()  const { return bc.size()*bc[0].size() + (is_augmented() ? 3 : 0); }
  const Int get_next_row()  const { return istart+get_num_rows(); }

  // assign the new strengths from BEM - do not let base class do this
  void set_str(const size_t ioffset, const size_t icnt, Vector<S> _in) {
    assert(vs.size() == 2 && "Vortex strength array not initialized");
    assert(_in.size() == vs[0].size()*2 && "Set strength array size does not match");
    assert(ioffset == 0 && "Offset is not zero");

    // copy over the strengths
    for (size_t i=0; i<get_npanels(); ++i) {
      vs[0][i] = _in[2*i+0];
      vs[1][i] = _in[2*i+1];
      //std::cout << "elem " << i << " with area " << area[i] << " has vs " << vs[0][i] << " " << vs[1][i] << std::endl;
    }

    // now recompute the absolute panel strengths
    vortex_sheet_to_panel_strength();
  }

  // convert vortex sheet strength to absolute panel strength
  void vortex_sheet_to_panel_strength() {
    // make sure (*s) even exists
    if (not this->s) {
      std::array<Vector<S>,3> new_s;
      this->s = std::move(new_s);
    }

    // make sure the arrays are sized properly
    for (size_t d=0; d<3; ++d) {
      (*this->s)[d].resize(get_npanels());
    }

    // convenience references
    std::array<Vector<S>,3>& x1 = b[0];
    std::array<Vector<S>,3>& x2 = b[1];

    // now copy the values over (do them all)
    for (size_t i=0; i<get_npanels(); ++i) {
      for (size_t d=0; d<Dimensions; ++d) {
        (*this->s)[d][i] = (vs[0][i]*x1[d][i] + vs[1][i]*x2[d][i]) * area[i];
      }
      //std::cout << "elem " << i << " has" << std::endl;
      //std::cout << "  x1 " << x1[0][i] << " " << x1[1][i] << " " << x1[2][i] << std::endl;
      //std::cout << "  x2 " << x2[0][i] << " " << x2[1][i] << " " << x2[2][i] << std::endl;
      //std::cout << "  vs " << vs[0][i] << " " << vs[1][i] << std::endl;
      //std::cout << "   s " << (*this->s)[0][i] << " " << (*this->s)[1][i] << " " << (*this->s)[2][i] << std::endl;
      //std::cout << "elem " << i << " has s " << (*this->s)[0][i] << " " << (*this->s)[1][i] << " " << (*this->s)[2][i] << std::endl;
    }
  }

  // a little logic to see if we should augment the BEM equations for this object
  const bool is_augmented() const {
    bool augment = true;
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
    //if (FORCE_NO_AUGMENTATION) augment = false;
    //augment = false;
    return augment;
  }

  // add more nodes and panels to this collection
  void add_new(const std::vector<S>&   _x,
               const std::vector<Int>& _idx,
               const std::vector<S>&   _val) {

    // remember old sizes of nodes and element arrays
    const size_t nnold = this->n;
    const size_t neold = get_npanels();

    // make sure input arrays are correctly-sized
    assert(_x.size() % Dimensions == 0 && "Position array is not an even multiple of dimensions");
    assert(_idx.size() % Dimensions == 0 && "Index array is not an even multiple of dimensions");
    const size_t nnodes = _x.size() / Dimensions;
    const size_t nsurfs = _idx.size() / Dimensions;

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
    idx.resize(2*neold + _idx.size());
    for (size_t i=0; i<2*nsurfs; ++i) {
      // make sure it exists in the nodes array
      if (_idx[i] >= nnold+nnodes) idx_are_all_good = false;
      idx[2*neold+i] = nnold + _idx[i];
    }
    assert(idx_are_all_good && "Some indicies are bad");

    // compute all basis vectors and panel areas
    compute_bases();

    // now, depending on the element type, put the value somewhere
    if (this->E == active) {
      // value is a fixed strength for the element
      assert(_val.size() == 2*nsurfs);
      // value is a fixed strength for the element
      for (size_t d=0; d<2; ++d) {
        vs[d].resize(neold+nsurfs);
        for (size_t i=0; i<nsurfs; ++i) {
          vs[d][neold+i] = _val[2*i+d];
        }
      }

      // and ensure that the raw strengths are resized and set
      vortex_sheet_to_panel_strength();

    } else if (this->E == reactive) {
      // value is a boundary condition
      // make sure we have the same number of components in the new array as in the old
      assert(bc.size() == _val.size()/nsurfs);
      // copy them into place
      for (size_t d=0; d<3; ++d) {
        bc[d].resize(neold+nsurfs);
        for (size_t i=0; i<nsurfs; ++i) {
          bc[d][neold+i] = _val[3*i+d];
        }
      }
      // upsize vortex sheet and raw strength arrays, too
      for (size_t d=0; d<2; ++d) {
        vs[d].resize(neold+nsurfs);
      }
      for (size_t d=0; d<3; ++d) {
        (*this->s)[d].resize(neold+nsurfs);
      }

    } else if (this->E == inert) {
      // value is ignored (probably zero)
    }

    // velocity is in the base class - just resize it here
    for (size_t d=0; d<Dimensions; ++d) {
      this->u[d].resize(nnold+nnodes);
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

    // make sure we've calculated transformed center (we do this when we do vol)
    assert(vol > 0.0 && "Have not calculated transformed center");
    // and we trust that we've transformed utc to tc

    // apply a factor times the body motion
    for (size_t i=0; i<this->get_n(); ++i) {

      // apply the translational velocity
      std::array<double,Dimensions> thisvel = this->B->get_vel(_time);
      for (size_t d=0; d<Dimensions; ++d) {
        this->u[d][i] += _factor * (float)thisvel[d];
      }

      // now compute the rotational velocity with respect to the geometric center - NEEDS WORK
      double thisrotvel = this->B->get_rotvel(_time);
      // center of this panel
      Int id0 = idx[2*i];
      Int id1 = idx[2*i+1];
      // panel center
      const S xc = 0.5 * (this->x[0][id1] + this->x[0][id0]);
      const S yc = 0.5 * (this->x[1][id1] + this->x[1][id0]);
      // add rotational velocity
      this->u[0][i] -= _factor * (float)thisrotvel * (yc - tc[1]);
      this->u[1][i] += _factor * (float)thisrotvel * (xc - tc[0]);
    }
  }
 
  void zero_strengths() {
    // call base class first
    ElementBase<S>::zero_strengths();

    // and reset the source strengths here
    if (ss) {
      std::fill(ss->begin(), ss->end(), 0.0);
    }
  }

  // augment the strengths with a value equal to that which accounts for
  //   the solid-body rotation of the object
  // NOTE: this needs to provide both the vortex AND source strengths!
  // AND we don't have the time - assume bodies have been transformed
  void add_rot_strengths(const S _constfac, const S _rotfactor) {

    // if no rotation, strengths, or no parent Body, or attached to ground, then no problem!
    if (not this->B) return;
    if (not this->s) return;
    if (std::string("ground").compare(this->B->get_name()) == 0) return;

    const S rotvel = (S)this->B->get_rotvel();
    //if (std::abs(rotvel) < std::numeric_limits<float>::epsilon()) return;

    // make sure we've calculated transformed center (we do this when we do vol)
    assert(vol > 0.0 && "Have not calculated transformed center");
    // and we trust that we've transformed utc to tc

    // have we made ss yet? or is it the right size?
    if (ss) {
      ss->resize(this->s->size());
    } else {
      // value is a fixed strength for the segment
      Vector<S> new_ss(this->s->size());
      ss = std::move(new_ss);
    }

    //std::cout << "Inside add_rot_strengths, sizes are: " << get_npanels() << " " << this->s->size() << " " << ss->size() << std::endl;
    assert(this->s->size() == get_npanels() && "Strength array is not the same as panel count");

    // what is the actual factor that we will add?
    const S factor = _constfac + rotvel*_rotfactor;

    // still here? let's do it. use the untransformed coordinates
    for (size_t i=0; i<get_npanels(); i++) {
      const size_t j   = idx[2*i];
      const size_t jp1 = idx[2*i+1];
      // vector from object geometric center to panel center
      const S dx = 0.5 * ((*this->ux)[0][j] + (*this->ux)[0][jp1]) - utc[0];
      const S dy = 0.5 * ((*this->ux)[1][j] + (*this->ux)[1][jp1]) - utc[1];
      // velocity of the panel center
      const S ui = -factor * dy;
      const S vi =  factor * dx;

      // panel tangential vector, fluid to the left, body to the right
      S panelx = (*this->ux)[0][jp1] - (*this->ux)[0][j];
      S panely = (*this->ux)[1][jp1] - (*this->ux)[1][j];
      const S panell = 1.0 / std::sqrt(panelx*panelx + panely*panely);
      panelx *= panell;
      panely *= panell;

      // the vortex strength - we ADD to the existing
      (*this->s)[i] += -1.0 * (ui*panelx + vi*panely);

      // the source strength
      (*ss)[i] += -1.0 * (ui*panely - vi*panelx);

      // debug print
      if (_rotfactor > 0.0 and false) {
        std::cout << "  panel " << i << " at " << dx << " " << dy
                  << " adds to vortex str " << (-1.0 * (ui*panelx + vi*panely))
                  << " and source str " << (-1.0 * (ui*panely - vi*panelx)) << std::endl;
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
    utc[0] = (S)xsum;
    utc[1] = (S)ysum;
    utc[2] = (S)zsum;

    std::cout << "    geom center is " << utc[0] << " " << utc[1] << " " << utc[2] << " and vol is " << vol << std::endl;
  }

  // need to maintain the 3x3 set of basis vectors for each panel
  // this also calculates the triangle areas
  void compute_bases() {

    // how many panels do we have now?
    const Int nnew = get_npanels();

    // how big is my set of basis vectors?
    const Int norig = b[0][0].size();

    // resize any vectors
    for (size_t i=0; i<3; ++i) {
      for (size_t j=0; j<3; ++j) {
        b[i][j].resize(nnew);
      }
    }
    area.resize(nnew);

    // we'll reuse these vectors
    std::array<S,3> x1, x2, norm;

    // update what we need
    for (size_t i=norig; i<nnew; ++i) {
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
      for (size_t j=0; j<3; ++j) b[0][j][i] = x1[j];
      for (size_t j=0; j<3; ++j) b[1][j][i] = x2[j];
      for (size_t j=0; j<3; ++j) b[2][j][i] = norm[j];

      //std::cout << "elem near " << this->x[0][id0] << " " << this->x[1][id0] << " " << this->x[2][id0] << " has norm " << b[2][0][i] << " " << b[2][1][i] << " " << b[2][2][i] << std::endl;
    }
  }

  // when transforming a body-bound object to a new time, we must also transform the geometric center
  void transform(const double _time) {
    // must explicitly call the method in the base class
    ElementBase<S>::transform(_time);

    if (this->B) {
      // prepare for the transform
      std::array<double,Dimensions> thispos = this->B->get_pos();
      const double theta = this->B->get_orient();
      const S st = std::sin(theta);
      const S ct = std::cos(theta);

      // transform the utc to tc here
      tc[0] = (S)thispos[0] + utc[0]*ct - utc[1]*st;
      tc[1] = (S)thispos[1] + utc[0]*st + utc[1]*ct;
      tc[2] = (S)thispos[2];

    } else {
      // transform the utc to tc here
      tc[0] = utc[0];
      tc[1] = utc[1];
      tc[2] = utc[2];
    }
  }


  void zero_vels() {
    // must explicitly call the method in the base class to zero the node vels
    ElementBase<S>::zero_vels();

    // but also zero the panel-center vels
    for (size_t d=0; d<Dimensions; ++d) {
      pu[d].resize(get_npanels());
      std::fill(pu[d].begin(), pu[d].end(), 0.0);
    }
  }

  void finalize_vels(const std::array<double,Dimensions>& _fs) {
    // must explicitly call the method in the base class, too
    ElementBase<S>::finalize_vels(_fs);

    // but also zero the panel-center vels
    const S factor = 0.25/M_PI;
    for (size_t d=0; d<Dimensions; ++d) {
      for (size_t i=0; i<pu[d].size(); ++i) {
        pu[d][i] = _fs[d] + pu[d][i] * factor;
      }
    }
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
  // offset is scaled by vdelta
  //
  std::vector<S> represent_as_particles(const S _offset, const S _vdelta) {

    // how many panels?
    const size_t num_pts = get_npanels();

    // init the output vector (x, y, z, sx, sy, sz, r)
    std::vector<S> px(num_pts*7);

    // get basis vectors
    std::array<Vector<S>,3>& norm = b[2];

    // how far above the surface
    const S dn = _offset * _vdelta;

    for (size_t i=0; i<num_pts; i++) {
      const Int id0 = idx[3*i];
      const Int id1 = idx[3*i+1];
      const Int id2 = idx[3*i+2];
      const Int idx = 7*i;
      // start at center of panel
      for (size_t j=0; j<3; ++j) px[idx+j] = (1./3.) * (this->x[j][id0] + this->x[j][id1] + this->x[j][id2]);
      // push out a fixed distance
      // this assumes properly resolved, vdelta and dt
      for (size_t j=0; j<3; ++j) px[idx+j] += dn * norm[j][i];
      // complete the element with a strength
      for (size_t j=0; j<3; ++j) px[idx+3+j] = (*this->s)[j][i];
      // and the core size
      px[idx+6] = _vdelta;
      //std::cout << "  new part at " << px[idx+0] << " " << px[idx+1] << " " << px[idx+2];
      //std::cout << "     with str " << px[idx+3] << " " << px[idx+4] << " " << px[idx+5] << std::endl;
    }

    return px;
  }

  // find the new peak strength magnitude
  void update_max_str() {
    S thismax = ElementBase<S>::get_max_str();

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
  std::array<S,3> get_total_circ(const double _time) {
    std::array<S,3> circ;
    circ.fill(0.0);

    if (this->s) {
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
      // HACK - this assumes z-axis rotation only
      circ[2] = 2.0 * vol * (S)this->B->get_rotvel(_time);
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

    if (this->s) {
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

    if (this->s) {
      for (size_t i=0; i<Dimensions; ++i) {
        glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i+4]);
        glBufferData(GL_ARRAY_BUFFER, 0, (*this->s)[i].data(), GL_STATIC_DRAW);
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

      const size_t ilen = idx.size()*sizeof(Int);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mgl->vbo[Dimensions]);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, ilen, idx.data(), GL_DYNAMIC_DRAW);

      // here is where we split on element type: active/reactive vs. inert
      if (this->E == inert) {
        // just don't upload strengths

      } else { // this->E is active or reactive
        // the strengths
        if (this->s) {
          const size_t slen = (*this->s)[0].size()*sizeof(S);
          for (size_t i=0; i<Dimensions; ++i) {
            glBindBuffer(GL_ARRAY_BUFFER, mgl->vbo[i+4]);
            glBufferData(GL_ARRAY_BUFFER, slen, (*this->s)[i].data(), GL_DYNAMIC_DRAW);
          }
        }
      }

      glBindVertexArray(0);

      // must tell draw call how many elements are there - or, really, how many indices
      mgl->num_uploaded = idx.size();
    }
  }

  // OpenGL3 stuff to draw triangles, called once per frame
  void drawGL(std::vector<float>& _projmat,
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

      // draw as lines
      glUseProgram(mgl->spo[0]);

      // upload the current projection matrix
      glUniformMatrix4fv(mgl->projmat_attribute, 1, GL_FALSE, _projmat.data());

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

  // element-wise variables special to triangular panels
  std::vector<Int>		   idx;	// indexes into the x array
  std::array<Vector<S>,3>	pu;	// panel-center velocities (ElementBase stores *node* properties)

  std::vector<Vector<S>>	bc;	// boundary condition for the elements (normal) or (x1,x2) or (x1,x2,normal)
  std::array<Vector<S>,2>	vs;	// vortex sheet strengths of the elements (x1,x2)
  std::optional<Vector<S>> 	ss;	// source strengths which represent the vel inf of the rotating volume

  std::array<std::array<Vector<S>,3>,3> b;  // transformed basis vectors: x1 is b[0], x2 is b[1], normal is b[2], normal x is b[2][0]
  Vector<S>                     area;

  // parameters for the encompassing body
  Int istart;	// index of first entry in RHS vector and A matrix

  S vol;			// volume of the body - for augmented BEM solution
  std::array<S,Dimensions> utc;		// untransformed geometric center
  std::array<S,Dimensions>  tc;		// transformed geometric center

private:
#ifdef USE_GL
  std::shared_ptr<GlState> mgl;
#endif
  float max_strength;
};

