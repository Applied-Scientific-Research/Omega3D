/*
 * Simulation.cpp - a class to control a 3d vortex particle sim
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "Simulation.h"
#include "Points.h"
#include "VtkXmlHelper.h"
//#include "Influence.h"
//#include "CollectionHelper.h"
//#include "Core.h"
#include "Split.h"

#include <cassert>
#include <cmath>
#include <variant>

// constructor
Simulation::Simulation()
  : re(100.0),
    dt(0.01),
    fs{0.0,0.0,0.0},
    vort(),
    bdry(),
    fldpt(),
    diff(),
    conv(),
    description(),
    time(0.0),
    output_dt(0.0),
    end_time(0.0),
    use_end_time(false),
    max_steps(0),
    use_max_steps(false),
    sim_is_initialized(false),
    step_has_started(false),
    step_is_finished(false)
  {}

// addresses for use in imgui
float* Simulation::addr_re() { return &re; }
float* Simulation::addr_dt() { return &dt; }
float* Simulation::addr_fs() { return fs; }

// getters
float Simulation::get_hnu() { return std::sqrt(dt/re); }
float Simulation::get_ips() { return diff.get_nom_sep_scaled() * get_hnu(); }
float Simulation::get_vdelta() { return diff.get_particle_overlap() * get_ips(); }
float Simulation::get_time() { return (float)time; }
float Simulation::get_end_time() { return (float)end_time; }
bool Simulation::using_end_time() { return use_end_time; }
size_t Simulation::get_max_steps() { return max_steps; }
bool Simulation::using_max_steps() { return use_max_steps; }
float Simulation::get_output_dt() { return (float)output_dt; }
std::string Simulation::get_description() { return description; }

// setters
void Simulation::set_description(const std::string desc) { description = desc; }
void Simulation::set_end_time(const double net) { end_time = net; use_end_time = true; }
void Simulation::set_max_steps(const size_t nms) { max_steps = nms; use_max_steps = true; }
void Simulation::set_output_dt(const double nodt) { output_dt = nodt; }

// status
//size_t Simulation::get_npanels() { return bdry.get_npanels(); }
size_t Simulation::get_nparts() {
  size_t n = 0;
  for (auto &coll: vort) {
    std::visit([&n](auto& elem) { n += elem.get_n(); }, coll);
  }
  return n;
}

size_t Simulation::get_nfldpts() {
  size_t n = 0;
  for (auto &coll : fldpt) {
    std::visit([&n](auto& elem) { n += elem.get_n(); }, coll);
  }
  return n;
}

// like a setter
void Simulation::set_re_for_ips(const float _ips) {
  re = std::pow(diff.get_nom_sep_scaled(), 2) * dt / pow(_ips, 2);
  diff.set_diffuse(false);
}

void Simulation::set_diffuse(const bool _do_diffuse) {
  diff.set_diffuse(_do_diffuse);
}

const bool Simulation::get_diffuse() {
  return diff.get_diffuse();
}

//void Simulation::set_amr(const bool _do_amr) {
//  diff.set_amr(_do_amr);
//  diff.set_diffuse(true);
//}

#ifdef USE_GL
void Simulation::initGL(std::vector<float>& _projmat,
                        float*              _poscolor,
                        float*              _negcolor,
                        float*              _defcolor) {
  for (auto &coll : vort) {
    std::visit([=, &_projmat](auto& elem) { elem.initGL(_projmat, _poscolor, _negcolor, _defcolor); }, coll);
  }
  //for (auto &coll : bdry) {
  //  std::visit([=, &_projmat](auto& elem) { elem.initGL(_projmat, _poscolor, _negcolor, _defcolor); }, coll);
  //}
  for (auto &coll : fldpt) {
    std::visit([=, &_projmat](auto& elem) { elem.initGL(_projmat, _poscolor, _negcolor, _defcolor); }, coll);
  }
}

void Simulation::updateGL() {
  for (auto &coll : vort) {
    std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  }
  //for (auto &coll : bdry) {
  //  std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  //}
  for (auto &coll : fldpt) {
    std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  }
}

void Simulation::drawGL(std::vector<float>& _projmat,
                        RenderParams&       _rparams) {

  if (step_is_finished) {
    _rparams.tracer_size = get_ips() * _rparams.tracer_scale;
    for (auto &coll : vort) {
      std::visit([&](auto& elem) { elem.drawGL(_projmat, _rparams); }, coll);
    }
    //for (auto &coll : bdry) {
    //  std::visit([&](auto& elem) { elem.drawGL(_projmat, _rparams); }, coll);
    //}
    for (auto &coll : fldpt) {
      std::visit([&](auto& elem) { elem.drawGL(_projmat, _rparams); }, coll);
    }
  }
}
#endif

bool Simulation::is_initialized() { return sim_is_initialized; }

void Simulation::set_initialized() { sim_is_initialized = true; }

void Simulation::reset() {

  // must wait for step() to complete, if it's still working
  if (stepfuture.valid()) {
    stepfuture.wait();
    stepfuture.get();
  }

  // now reset everything else
  time = 0.0;
  vort.clear();
  bdry.clear();
  fldpt.clear();
  sim_is_initialized = false;
  step_has_started = false;
  step_is_finished = false;
}

//void Simulation::clear_bodies() {
//  bodies.clear();
//}

// Write a set of vtu files for the particles and panels
void Simulation::write_vtk() {
  static size_t frameno = 0;

  size_t idx = 0;
  for (auto &coll : vort) {
    // eventually all collections will support vtk output
    //std::visit([=](auto& elem) { elem.write_vtk(); }, coll);
    // only proceed if the collection is Points
    if (std::holds_alternative<Points<float>>(coll)) {
      Points<float>& pts = std::get<Points<float>>(coll);
      write_vtu_points<float>(pts, idx++, frameno);
    }
  }

  idx = 0;
  for (auto &coll : fldpt) {
    // eventually all collections will support vtk output
    //std::visit([=](auto& elem) { elem.write_vtk(); }, coll);
    // only proceed if the collection is Points
    if (std::holds_alternative<Points<float>>(coll)) {
      Points<float>& pts = std::get<Points<float>>(coll);
      write_vtu_points<float>(pts, idx++, frameno);
    }
  }

  frameno++;
}

//
// Check all aspects of the simulation for conditions that should stop the run
//
//std::string Simulation::check_simulation(const size_t _nff, const size_t _nbf) {
std::string Simulation::check_simulation(const size_t _nff) {
  std::string retstr;

  // Check for no bodies and no particles
  if (_nff == 0) {
    //retstr.append("No flow features and no bodies - try adding one or both.\n");
    retstr.append("No flow features. Add one, reset, and run.\n");
  }

  // Check for a body and no particles and no freestream
  // retstr.append("No flow features and zero freestream speed - try adding one or both.\n");

  // Check for excessive elongation
  float max_elong = 0.0;
  for (auto &coll: vort) {
    // only check particles ("Points")
    if (std::holds_alternative<Points<float>>(coll)) {
      Points<float>& pts = std::get<Points<float>>(coll);
      max_elong = std::max(max_elong, pts.get_max_elong());
    }
    //std::visit([&](auto& elem) { max_elong = std::max(max_elong, elem.get_max_elong()); }, coll);
  }
  if (max_elong > 1.5) retstr.append("Elongation threshold exceeded! Reset and reduce the time step size.\n");

  // Check for vorticity-based Courant number - too high and we're probably messing up
  // later

  return retstr;
}

//
// query and get() the future if possible
//
bool Simulation::test_for_new_results() {

  if (not step_has_started) {
    // if we haven't made an async call yet
    return true;

  } else if (is_future_ready(stepfuture)) {
    // if we did, and it's ready to give us the results
    stepfuture.get();

#ifdef USE_GL
    // tell flow objects to update their values to the GPU
    updateGL();
#endif

    // set flag indicating that at least one step has been solved
    step_is_finished = true;
    step_has_started = false;

    return true;
  }

  // async call is not finished, do not try calling it again
  return false;
}

//
// call this from a real-time GUI - will only run a new step if the last one is done
//
void Simulation::async_step() {
  step_has_started = true;
  stepfuture = std::async(std::launch::async, [this](){step();});
}

//
// here's the vortex method: convection and diffusion with operator splitting
//
void Simulation::step() {
  std::cout << std::endl << "Taking step at t=" << time << " with n=" << get_nparts() << std::endl;

  // we wind up using this a lot
  std::array<double,3> thisfs = {fs[0], fs[1], fs[2]};

  // for simplicity's sake, just run one full diffusion step here
  //diff.step(dt, re, get_vdelta(), thisfs, vort, bdry);
  diff.step(dt, re, thisfs, vort, bdry);

  // operator splitting requires one half-step diffuse (use coefficients from previous step, if available)
  //diff.step(0.5*dt, get_vdelta(), get_ips(), thisfs, vort, bdry);


  // advect with no diffusion (must update BEM strengths)
  //conv.advect_1st(dt, thisfs, vort, bdry, fldpt);
  conv.advect_2nd(dt, thisfs, vort, bdry, fldpt);

  // operator splitting requires another half-step diffuse (must compute new coefficients)
  //diff.step(0.5*dt, get_vdelta(), get_ips(), thisfs, vort, bdry);

  // step complete, now split any elongated particles
  for (auto &coll: vort) {
  
    // but only check particles ("Points")
    if (std::holds_alternative<Points<float>>(coll)) {

      Points<float>& pts = std::get<Points<float>>(coll);
      //std::cout << "    check split for " << pts.get_n() << " particles" << std::endl;
      std::cout << std::endl;

      // none of these are passed as const, because both may be extended with new particles
      std::array<Vector<float>,Dimensions>& x = pts.get_pos();
      Vector<float>&                        r = pts.get_rad();
      Vector<float>&                    elong = pts.get_elong();
      std::array<Vector<float>,Dimensions>& s = pts.get_str();

      // last two arguments are: relative distance, allow variable core radii
      (void)split_elongated<float>(x[0], x[1], x[2], r, elong, s[0], s[1], s[2],
                               diff.get_core_func(),
                               diff.get_particle_overlap(),
                               1.2);

      // we probably have a different number of particles now, resize the u, ug, elong arrays
      pts.resize(r.size());
    }
  }

  // update strength for coloring purposes (eventually should be taken care of automatically)
  //vort.update_max_str();

  // update dt and return
  time += (double)dt;
}

// set up the particles
// TODO - accept elem_t and move_t from the caller!
void Simulation::add_particles(std::vector<float> _invec) {

  if (_invec.size() == 0) return;

  // make sure we're getting full particles
  assert(_invec.size() % 7 == 0);

  // add the vdelta to each particle and pass it on
  const float thisvd = get_vdelta();
  for (size_t i=6; i<_invec.size(); i+=7) {
    _invec[i] = thisvd;
  }

  // if no collections exist
  if (vort.size() == 0) {
    // make a new collection
    vort.push_back(Points<float>(_invec, active, lagrangian));      // vortons

    // some examples of other collections
    //vort.push_back(Points<float>(5000, active, lagrangian));      // vortons
    //fldpt.push_back(Points<float>(2000, inert, lagrangian));      // tracer particles
    //fldpt.push_back(Points<float>(100, inert, fixed));            // static field points
    //bdry.push_back(Panels<float>(500, reactive, bodybound));    // panels

  } else {
    // THIS MUST USE A VISITOR

    // helper struct for dispatching through a variant
    //struct AddElemsVisitor {
    //  void operator()(Points<float> const& coll) { coll.add_new(what); }
    //  void operator()(Panels<float> const& coll) { coll.add_new(src); }
    //} visitor;

    // HACK - add all particles to first collection
    std::visit([&](auto& elem) { elem.add_new(_invec); }, vort.back());

    // TODO - need to find a way to not always make a new collection!
    //        like, test a collection for matching elem_t and move_t and Points/Panels/etc
    //        and add to the most appropriate
    //for (auto &coll : vort) {
      //std::visit([&](auto& elem) { elem.add_new(_invec); }, coll);
    //}
  }
}

// add some tracer particles to new arch
void Simulation::add_fldpts(std::vector<float> _xyz, const bool _moves) {

  if (_xyz.size() == 0) return;

  // make sure we're getting full points
  assert(_xyz.size() % 3 == 0);

  const move_t move_type = _moves ? lagrangian : fixed;

  // add to new archtecture

  // if no collections exist
  if (fldpt.size() == 0) {
    // make a new collection
    fldpt.push_back(Points<float>(_xyz, inert, move_type));

  } else {
    // THIS MUST USE A VISITOR
    // HACK - add all particles to first collection
    //std::visit([&](auto& elem) { elem.add_new(_xyz); }, fldpt.back());
    auto& coll = fldpt.back();
    // eventually we will want to check every collection for matching element and move types
    // only proceed if the collection is Points
    if (std::holds_alternative<Points<float>>(coll)) {
      Points<float>& pts = std::get<Points<float>>(coll);
      pts.add_new(_xyz);
    }
  }
}

// add geometry
/*
void Simulation::add_boundary(bdryType _type, std::vector<float> _params) {
  if (_type == circle) {
    assert(_params.size() == 3);
    bdry.add(Circle<float>(_params[0], _params[1], _params[2]));
  }
}
*/

