
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

  //static bool sim_is_running = false;
  //static bool begin_single_step = false;
  //const solution_t solver = direct_cpu;
  //std::array<double,Dimensions> fs = {0.0, 0.0, 0.0};
  //double time = 0.0;
  //const double dt = 0.01;

  // for starters, generate some vortons, particles, and field points
  ffeatures.emplace_back(std::make_unique<BlockOfRandom>(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.01, 10000));
  //ffeatures.emplace_back(std::make_unique<BlockOfRandom>(10000, active, lagrangian));
  //ffeatures.emplace_back(std::make_unique<BlockOfRandom>(5000, inert, lagrangian));
  //ffeatures.emplace_back(std::make_unique<BlockOfRandom>(2000, inert, fixed));
 
  // Main loop
  while (true) {

    // if particles are not yet created, make them
    if (not sim.is_initialized()) {
      std::cout << std::endl << "Initializing simulation" << std::endl;

      // initialize particle distributions
      //for (auto const& ff: ffeatures) {
      //  sim.add_particles( ff->init_particles(sim.get_ips()) );
      //}

      // initialize solid objects
      //for (auto const& bf : bfeatures) {
      //  sim.add_boundary( bf->get_type(), bf->get_params() );
      //}

      // initialize panels
      //sim.init_bcs();
      //sim.set_initialized();
    }

    // begin a dynamic step: convection and diffusion
    sim.step();

    // for testing: always break after one step
    break;
  }

  return 0;
}

