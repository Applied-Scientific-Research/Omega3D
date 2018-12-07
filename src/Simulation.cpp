/*
 * Simulation.cpp - a class to control a 3D vortex particle sim
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "Simulation.h"
#include "Points.h"
#include "Influence.h"
#include "CollectionHelper.h"

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
    time(0.0),
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

// status
//size_t Simulation::get_npanels() { return bdry.get_npanels(); }
size_t Simulation::get_nparts() {
  size_t n = 0;
  for (auto &targ: vort) {
    //n += std::visit([=](auto& elem) { return elem.getn(); }, targ);
    std::visit([&n](auto& elem) { n += elem.getn(); }, targ);
  }
  return n;
}

size_t Simulation::get_n(std::vector<Collection>& _collect) {
  size_t n = 0;
  for (auto &targ: _collect) {
    //n += std::visit([=](auto& elem) { return elem.getn(); }, targ);
    std::visit([&n](auto& elem) { n += elem.getn(); }, targ);
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

//void Simulation::set_amr(const bool _do_amr) {
//  diff.set_amr(_do_amr);
//  diff.set_diffuse(true);
//}

#ifdef USE_GL
void Simulation::initGL(std::vector<float>& _projmat,
                        float*              _poscolor,
                        float*              _negcolor) {
  // never gets called
  //bdry.initGL(_projmat, _poscolor, _negcolor);
  //vort.initGL(_projmat, _poscolor, _negcolor);
}
void Simulation::updateGL() {
  //bdry.updateGL();
  //vort.updateGL();
  for (auto &coll : vort) {
    std::visit([=](auto& elem) { elem.updateGL(); }, coll);
  }
}
void Simulation::drawGL(std::vector<float>& _projmat,
                        float*              _poscolor,
                        float*              _negcolor) {
  if (step_is_finished) {
    //vort.drawGL(_projmat, _poscolor, _negcolor);
    for (auto &coll : vort) {
      std::visit([=, &_projmat](auto& elem) { elem.drawGL(_projmat, _poscolor, _negcolor); }, coll);
    }
    //bdry.drawGL(_projmat, _poscolor, _negcolor);
  }
}
#endif

//
// main must indicate that panels should be made
//   because initGL and updateGL need to send data soon
//
/*
void Simulation::init_bcs() {
  // are panels even made? do this first
  bdry.make_panels(get_ips());
}
*/

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
  //bdry.clear();
  sim_is_initialized = false;
  step_has_started = false;
  step_is_finished = false;
}

//
// Check all aspects of the simulation for conditions that should stop the run
//
std::string Simulation::check_simulation() {
  std::string retstr;

  // Check for no bodies and no particles
  // retstr.append("No flow features and no bodies - try adding one or both.\n");

  // Check for a body and no particles and no freestream
  // retstr.append("No flow features and zero freestream speed - try adding one or both.\n");

  // Check for excessive elongation
  float max_elong = 0.0;
  for (auto &coll: vort) {
    //std::visit([&](auto& elem) { max_elong = std::max(max_elong, elem.get_max_elong(); }, coll);
  }
  if (max_elong > 2.0) retstr.append("Elongation threshold exceeded! Reset and reduce the time step size.\n");

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
  std::cout << std::endl << "Taking step at t=" << time << " with n=" << get_n(vort) << std::endl;

  // we wind up using this a lot
  //std::array<double,3> thisfs = reinterpret_cast<std::array<float,3>&>(fs);
  std::array<double,3> thisfs = {fs[0], fs[1], fs[2]};


  // are panels even made? do this first
  //bdry.make_panels(get_ips());

  // for simplicity's sake, just run one full diffusion step here
  //diff.step(dt, re, get_vdelta(), thisfs, vort, bdry);

  // operator splitting requires one half-step diffuse (use coefficients from previous step, if available)
  //diff.step(0.5*dt, get_vdelta(), get_ips(), thisfs, vort, bdry);


  // advect with no diffusion (must update BEM strengths)
  //conv.advect_1st(dt, thisfs, vort, bdry, fldpt);
  conv.advect_2nd(dt, thisfs, vort, bdry, fldpt);

  // operator splitting requires another half-step diffuse (must compute new coefficients)
  //diff.step(0.5*dt, get_vdelta(), get_ips(), thisfs, vort, bdry);

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

// add geometry
/*
void Simulation::add_boundary(bdryType _type, std::vector<float> _params) {
  if (_type == circle) {
    assert(_params.size() == 3);
    bdry.add(Circle<float>(_params[0], _params[1], _params[2]));
  }
}
*/

