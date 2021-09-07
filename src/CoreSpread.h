/*
 * CoreSpread.h - the pure core-spreading method for diffusion in 3D
 *
 * (c)2020-1 Applied Scientific Research, Inc.
 *           Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Core.h"
#include "VectorHelper.h"

#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <array>
#include <vector>
#include <random>

//
// Class to hold CoreSpread parameters and temporaries
//
// templatized on storage type ST, compute type CT, and max number of moments to solve for
//
template <class ST>
class CoreSpread {
public:
  CoreSpread();

  // all-to-all diffuse; can change array sizes
  void diffuse_all(const std::array<Vector<ST>,3>&,
                   Vector<ST>&,
                   const ST,
                   const CoreType);

  void from_json(const nlohmann::json);
  void add_to_json(nlohmann::json&) const;

private:

};

// primary constructor
template <class ST>
CoreSpread<ST>::CoreSpread() {}


//
// Apply the random vortex method to the particles
//
template <class ST>
void CoreSpread<ST>::diffuse_all(const std::array<Vector<ST>,3>& pos,
                                 Vector<ST>& rad,
                                 const ST h_nu,
                                 const CoreType core_func) {

  // make sure all vector sizes are identical
  assert(pos[0].size()==pos[1].size() && "Input arrays are not uniform size");
  assert(pos[0].size()==rad.size() && "Input arrays are not uniform size");
  const size_t n = rad.size();

  std::cout << "  Running CoreSpread with n " << n << std::endl;

  // start timer
  auto start = std::chrono::system_clock::now();

  const ST core_second_mom = get_core_second_mom<ST>(core_func);

  // for each particle (can parallelize this part)
  // note that an OpenMP loop here will need to use int32_t as the counter variable type
  for (size_t i=0; i<n; ++i) {
    // increase the core radius
    rad[i] = std::sqrt( std::pow(rad[i], 2) + (2.0/core_second_mom)*std::pow(h_nu, 2));
  }

  // finish timer and report
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  printf("    corespread.diffuse_all:\t[%.4f] seconds\n", (float)elapsed_seconds.count());
}

//
// read/write parameters to json
//

// create and write a json object for all diffusion parameters
template <class ST>
void CoreSpread<ST>::from_json(const nlohmann::json simj) {

  /*
  if (simj.find("CoreSpread") != simj.end()) {
    nlohmann::json j = simj["CoreSpread"];

    if (j.find("ignoreBelow") != j.end()) {
      ignore_thresh = j["ignoreBelow"];
      std::cout << "  setting ignore_thresh= " << ignore_thresh << std::endl;
    }

    if (j.find("relativeThresholds") != j.end()) {
      thresholds_are_relative = j["relativeThresholds"];
      std::cout << "  setting thresholds_are_relative= " << thresholds_are_relative << std::endl;
    }
  }
  */
}

// create and write a json object for all diffusion parameters
template <class ST>
void CoreSpread<ST>::add_to_json(nlohmann::json& simj) const {

  /*
  // set corespread-specific parameters
  nlohmann::json j;
  j["ignoreBelow"] = ignore_thresh;
  j["relativeThresholds"] = thresholds_are_relative;
  simj["CoreSpread"] = j;
  */
}

