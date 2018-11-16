
#include "Omega3D.h"
#include "Collection.h"
#include "Points.h"
#include "Influence.h"

#include <iostream>
#include <vector>
#include <variant>


// execution starts here

int main(int argc, char const *argv[]) {

  //const solution_t solver = direct_cpu;
  std::array<double,Dimensions> fs = {0.0, 0.0, 0.0};
  double time = 0.0;
  const double dt = 0.01;
 
  std::vector<Collection> vort;		// the free vorticity
  std::vector<Collection> bdry;		// all boundaries
  std::vector<Collection> fldpt;	// tracers and field points
 
  vort.push_back(Points<float>(5000, active, lagrangian));	// vortons
  fldpt.push_back(Points<float>(2000, inert, lagrangian));	// tracer particles
  fldpt.push_back(Points<float>(100, inert, fixed));		// static field points
  //bdry.push_back(Panels<float>(500, reactive, bodybound));	// panels

  // need this for dispatching velocity influence calls, template param is accumulator type
  // should the solution_t be an argument to the constructor?
  InfluenceVisitor<float> visitor;


  // one-half diffusion step
  //tbd


  // this is one Euler convection step - how would we do one 2nd order step?
  // can we have an object for different forward integrators?

  if (bdry.size() > 0) std::cout << std::endl << "Solving for BEM RHS" << std::endl << std::endl;

  // loop over boundary collections
  for (auto &targ: bdry) {
    std::cout << "  Solving for velocities on" << to_string(targ) << std::endl;
    // zero velocities
    std::visit([=](auto& elem) { elem.zero_vels(); }, targ);
    // accumulate from vorticity
    for (auto &src: vort) {
      std::visit(visitor, src, targ);
    }
  }
 
  if (bdry.size() > 0) std::cout << std::endl << "Solving for BEM matrix" << std::endl << std::endl;

  // loop over boundary collections
  for (auto &targ: bdry) {
    std::cout << "  Solving for influence coefficients on" << to_string(targ) << std::endl;
    // assemble from all boundaries
    for (auto &src: bdry) {
      std::visit(visitor, src, targ);
    }
  }
 
  //std::cout << std::endl << "Solving BEM for strengths" << std::endl << std::endl;

  if (vort.size()+fldpt.size() > 0) std::cout << std::endl << "Solving for velocities" << std::endl << std::endl;

  // find the influence on every vorticity element
  for (auto &targ: vort) {
    std::cout << "  Solving for velocities on" << to_string(targ) << std::endl << std::flush;
    // zero velocities
    std::visit([=](auto& elem) { elem.zero_vels(); }, targ);
    // accumulate from vorticity
    for (auto &src: vort) {
      // how do I specify the solver?
      std::visit(visitor, src, targ);
    }
    // accumulate from boundaries
    for (auto &src: bdry) {
      std::visit(visitor, src, targ);
    }
    // add freestream and divide by 2pi
    std::visit([=](auto& elem) { elem.finalize_vels(fs); }, targ);
  }
 
  // find the influence on every field point/tracer element
  for (auto &targ: fldpt) {
    std::cout << "  Solving for velocities on" << to_string(targ) << std::endl;
    // zero velocities
    std::visit([=](auto& elem) { elem.zero_vels(); }, targ);
    // accumulate from vorticity
    for (auto &src: vort) {
      std::visit(visitor, src, targ);
    }
    // accumulate from boundaries
    for (auto &src: bdry) {
      std::visit(visitor, src, targ);
    }
    // add freestream and divide by 2pi
    std::visit([=](auto& elem) { elem.finalize_vels(fs); }, targ);
  }

  std::cout << std::endl << "Convection step" << std::endl << std::endl;

  // move every movable element
  for (auto &coll: vort) {
    std::visit([=](auto& elem) { elem.move(dt); }, coll);
  }
  for (auto &coll: bdry) {
    std::visit([=](auto& elem) { elem.move(dt); }, coll);
  }
  for (auto &coll: fldpt) {
    std::visit([=](auto& elem) { elem.move(dt); }, coll);
  }
 
  // one-half diffusion step
  //tbd

  // write output, restart files
  time += dt;

  std::cout << std::endl << "Done" << std::endl << std::endl;

  return 0;
}

