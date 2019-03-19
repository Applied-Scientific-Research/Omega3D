/*
 * main_batch.cpp - Driver code for Omega3D + Vc vortex particle method
 *                  and boundary element method solver, batch version
 *
 * (c)2017-9 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "FlowFeature.h"
#include "Simulation.h"
#include "RenderParams.h"

#include <iostream>
#include <vector>


// execution starts here

int main(int argc, char const *argv[]) {
  std::cout << std::endl << "Omega3D Batch" << std::endl;

  // Set up vortex particle simulation
  Simulation sim;
  std::vector< std::unique_ptr<FlowFeature> > ffeatures;
  //std::vector< std::unique_ptr<BoundaryFeature> > bfeatures;

  size_t nsteps = 0;
  RenderParams rparams;

  sim.set_diffuse(false);
  sim.set_re_for_ips(0.1);
  *(sim.addr_dt()) = 0.1;

  // a string to hold any error messages
  std::string sim_err_msg;

  // for starters, generate some vortons, particles, and field points
  //ffeatures.emplace_back(std::make_unique<SingularRing>(0.0, 0.0, 0.1, 0.0, 0.0, -1.0, 1.0, 1.0));
  //ffeatures.emplace_back(std::make_unique<SingularRing>(0.0, 0.0, -0.1, 0.0, 0.0, 1.0, 1.0, 1.0));
  ffeatures.emplace_back(std::make_unique<BlockOfRandom>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1000.0, 30000));
  //ffeatures.emplace_back(std::make_unique<BlockOfRandom>(10000, active, lagrangian));
  //ffeatures.emplace_back(std::make_unique<BlockOfRandom>(5000, inert, lagrangian));
  //ffeatures.emplace_back(std::make_unique<BlockOfRandom>(2000, inert, fixed));
 

  std::cout << std::endl << "Initializing simulation" << std::endl;

  // initialize particle distributions
  for (auto const& ff: ffeatures) {
    sim.add_particles( ff->init_particles(sim.get_ips()) );
  }

  // initialize solid objects
  //for (auto const& bf : bfeatures) {
  //  sim.add_boundary( bf->get_body(), bf->init_elements(sim.get_ips()) );
  //}

  // initialize measurement features
  //for (auto const& mf: mfeatures) {
  //  sim.add_fldpts( mf->init_particles(0.1*sim.get_ips()), mf->moves() );
  //}

  sim.set_initialized();

  //
  // Main loop
  //

  while (true) {

    // check flow for blow-up or errors
    sim_err_msg = sim.check_simulation(ffeatures.size());

    if (sim_err_msg.empty()) {
      // the last simulation step was fine, OK to continue

      // generate new particles from emitters
      //for (auto const& ff : ffeatures) {
      //  sim.add_particles( ff->step_particles(sim.get_ips()) );
      //}

      // begin a new dynamic step: convection and diffusion
      sim.step();

    } else {
      // the last step had some difficulty
      std::cout << std::endl << "ERROR: " << sim_err_msg;

      // stop the run
      break;
    }

    nsteps++;

    // for testing: always break after a few steps
    if (nsteps == 2) break;

    // check vs. stopping conditions
    if (sim.using_max_steps() and sim.get_max_steps() == nsteps) break;
    if (sim.using_end_time() and sim.get_end_time() <= sim.get_time()) break;

  } // end step

  sim.reset();
  std::cout << "Quitting" << std::endl;

  return 0;
}

