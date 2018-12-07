/*
 * Diffusion.h - a class for diffusion of strength from bodies to particles and among particles
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega3D.h"
#include "Core.h"
//#include "Body.h"
#include "Merge.h"
//#include "Panels.h"
//#include "Particles.h"
//#include "Reflect.h"
#include "Collection.h"
#include "CollectionHelper.h"
#include "VectorHelper.h"
#include "VRM.h"

#include <cstdlib>
//#include <iostream>
#include <vector>
#include <cassert>


//
// One step of convection of elements, allowing for boundary conditions
//
// templatized on 'S'torage, 'A'ccumulator, and 'I'ndex types
//
template <class S, class A, class I>
class Diffusion {
public:
  Diffusion()
    : vrm(),
      core_func(gaussian),
      is_inviscid(false),
      adaptive_radii(false),
      nom_sep_scaled(std::sqrt(8.0)),
      particle_overlap(1.5)
    {}

  void set_diffuse(const bool _do_diffuse) { is_inviscid = not _do_diffuse; }
  //void set_amr(const bool _do_amr) { adaptive_radii = _do_amr; }
  S get_nom_sep_scaled() const { return nom_sep_scaled; }
  S get_particle_overlap() const { return particle_overlap; }

  void step(const float,
            const float,
            const std::array<double,3>&,
            std::vector<Collection>&,
            std::vector<Collection>&);

private:
  // the VRM algorithm, template params are storage, compute, max moments
  // note that NNLS needs doubles for its compute type or else it will fail
  VRM<S,A,2> vrm;

  // other necessary variables
  CoreType core_func;

  // toggle inviscid
  bool is_inviscid;

  // toggle adaptive particle sizes
  bool adaptive_radii;

  // nominal separation normalized by h_nu
  S nom_sep_scaled;

  // particle core size is nominal separation times this
  S particle_overlap;
};


//
// take a diffusion step
//
template <class S, class A, class I>
void Diffusion<S,A,I>::step(const float                 _dt,
                            const float                 _re,
                            const std::array<double,3>& _fs,
                            std::vector<Collection>&    _vort,
                            std::vector<Collection>&    _bdry) {

  if (is_inviscid) return;

  std::cout << "  inside Diffusion::step with dt=" << _dt << std::endl;

/*
  // always re-run the BEM calculation before shedding
  if (_bdry.exists()) {

    // set up and performs the BEM
    _bdry.reset_vels();
    add_influence<S,A,I>(_vort, _bdry);
    _bdry.scale_and_add_freestream(_fs);
    _bdry.find_strengths();

    // generate particles just above the surface
    std::vector<S> newparts = _bdry.get_panels().diffuse_onto(0.0001*(S)_dt, _re, _vdelta);

    // add those particles to the main particle list
    _vort.add_new(newparts);
  }
*/

  //
  // diffuse strength among existing particles
  //

  // ensure that we have a current h_nu
  vrm.set_hnu(std::sqrt(_dt/_re));

  // ensure that it knows to allow or disallow adaptive radii
  //vrm.set_adaptive_radii(adaptive_radii);

  // loop over active vorticity
  for (auto &coll: _vort) {

    // but only check particles ("Points")
    if (std::holds_alternative<Points<float>>(coll)) {

      Points<float>& pts = std::get<Points<float>>(coll);
      std::cout << "    computing diffusion among " << pts.getn() << " particles" << std::endl;

      // none of these are passed as const, because both may be extended with new particles
      std::array<Vector<S>,Dimensions>& x = pts.get_pos();
      Vector<S>&                        r = pts.get_rad();
      std::array<Vector<S>,Dimensions>& s = pts.get_str();

      // and make vectors for the new values
      Vector<S> newr = r;
      Vector<S> dsx, dsy, dsz;
      dsx.resize(r.size());
      dsy.resize(r.size());
      dsz.resize(r.size());

      // finally call VRM
      vrm.diffuse_all(x[0], x[1], x[2], r, newr, s[0], s[1], s[2], dsx, dsy, dsz, core_func, particle_overlap);

      // apply the strength change to the particles
      //elem->increment_in_place();
      assert(dsx.size()==s[0].size());
      for (size_t i=0; i<s[0].size(); ++i) {
        s[0][i] += dsx[i];
      }
      assert(dsx.size()==s[1].size());
      for (size_t i=0; i<s[1].size(); ++i) {
        s[1][i] += dsy[i];
      }
      assert(dsx.size()==s[2].size());
      for (size_t i=0; i<s[2].size(); ++i) {
        s[2][i] += dsz[i];
      }

      // and update the strength
      //elem->update_max_str();

      // we probably have a different number of particles now, resize the u, ug, elong arrays
      pts.resize(r.size());
    }
  }

  // reflect interior particles to exterior because VRM does not account for panels
/*
  if (_bdry.exists()) {
    (void) reflect<S,I>(_bdry, _vort);
  }
*/

  //
  // diffuse strength from boundaries/bodies
  //

  // use those BEM strengths to shed now
  //if (_bdry.exists()) {
    // initialize the new particles vector
  //  std::vector<S> newparts = _bdry.get_panels().diffuse_onto((S)_dt, _re, _vdelta);

    // finally, add those particles to the main particle list
  //  _vort.add_new(newparts);
  //}

  //
  // merge any close particles to clean up potentially-dense areas
  //
  for (auto &coll: _vort) {

    // but only check particles ("Points")
    if (std::holds_alternative<Points<float>>(coll)) {

      Points<float>& pts = std::get<Points<float>>(coll);
      std::cout << "    merging among " << pts.getn() << " particles" << std::endl;

      // none of these are passed as const, because both may be extended with new particles
      std::array<Vector<S>,Dimensions>& x = pts.get_pos();
      Vector<S>&                        r = pts.get_rad();
      std::array<Vector<S>,Dimensions>& s = pts.get_str();

      // last two arguments are: relative distance, allow variable core radii
      (void)merge_close_particles<S>(x[0], x[1], x[2], r, s[0], s[1], s[2], 
                                     particle_overlap,
                                     0.3);

      // we probably have a different number of particles now, resize the u, ug, elong arrays
      pts.resize(r.size());
    }
  }

  //
  // clean up by removing the innermost layer - the one that will be represented by boundary strengths
  //
/*
  if (_bdry.exists()) {
    // may need to do this multiple times to clear out concave zones!

    // new way
    (void) clear_inner_layer<S,I>(_bdry, _vort, _vdelta/particle_overlap);

    // old way
    //for (auto & elem: _vort.get_collections()) {
      //std::vector<S> partmod = _bdry.get_panels().remove_inner_layer(_nomsep, elem->get_x());
      //_p.modify_particles(partmod);
    //}
  }
*/

  //if (n>0) std::cout << "  part 0 with str " << x[2] << " is at " << x[0] << " " << x[1] << std::endl;
}

