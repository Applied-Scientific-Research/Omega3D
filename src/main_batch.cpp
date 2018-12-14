/*
 * main_batch.cpp - Driver program for batch version of Omega3D - The Vorticity Flow Solver
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#include "Omega3D.h"
#include "Simulation.h"
#include "FlowFeature.h"

#include <iostream>
#include <vector>


// execution starts here

int main(int argc, char const *argv[]) {
  std::cout << std::endl << "Omega3D Batch Solver" << std::endl;

  // Set up vortex particle simulation
  Simulation sim;
  std::vector< std::unique_ptr<FlowFeature> > ffeatures;
  //std::vector< std::unique_ptr<BoundaryFeature> > bfeatures;

  size_t nsteps = 0;
  //static bool sim_is_running = false;
  //static bool begin_single_step = false;
  //const solution_t solver = direct_cpu;
  //std::array<double,Dimensions> fs = {0.0, 0.0, 0.0};
  //double time = 0.0;
  //const double dt = 0.01;
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
 
  // Main loop
  while (true) {

    // if particles are not yet created, make them
    if (not sim.is_initialized()) {
      std::cout << std::endl << "Initializing simulation" << std::endl;

      // initialize particle distributions
      for (auto const& ff: ffeatures) {
        sim.add_particles( ff->init_particles(sim.get_ips()) );
      }

      // initialize solid objects
      //for (auto const& bf : bfeatures) {
      //  sim.add_boundary( bf->get_type(), bf->get_params() );
      //}

      // initialize panels
      //sim.init_bcs();
      sim.set_initialized();
    }

    // check flow for blow-up or errors
    sim_err_msg = sim.check_simulation(ffeatures.size());

    if (sim_err_msg.empty()) {
      // no errors, OK to continue

      // generate new particles from emitters
      //for (auto const& ff : ffeatures) {
      //  sim.add_particles( ff->step_particles(sim.get_ips()) );
      //}

      // begin a dynamic step: convection and diffusion
      // no need for async call in the batch program
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
  }
  std::cout << std::endl << "Done" << std::endl;

  return 0;
}

