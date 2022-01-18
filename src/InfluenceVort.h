/*
 * InfluenceVort.h - Non-class influence calculations
 *
 * (c)2020 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include "Omega3D.h"
#include "VectorHelper.h"
#include "Kernels.h"
#include "Points.h"
#include "Surfaces.h"
#include "ResultsType.h"
#include "ExecEnv.h"

#ifdef USE_VC
#include <Vc/Vc>
#endif

#include <iostream>
#include <vector>
#include <memory>
#include <optional>
#include <chrono>
#include <cmath>
#include <cassert>


//
// Vc and x86 versions of Points/Particles affecting Points/Particles
//
template <class S, class A>
void points_affect_points_vorticity (const Points<S>& src, Points<S>& targ, const ExecEnv& env) {

  std::cout << "    in ptptvort with" << env.to_string() << std::endl;

  auto start = std::chrono::system_clock::now();
  float flops = (float)targ.get_n();

  // get references to use locally
  const std::array<Vector<S>,Dimensions>& sx = src.get_pos();
  const Vector<S>&                        sr = src.get_rad();
  const std::array<Vector<S>,Dimensions>& ss = src.get_str();

  const std::array<Vector<S>,Dimensions>& tx = targ.get_pos();
  std::array<Vector<S>,Dimensions>&       tw = targ.get_vort();

  // We need 2 different loops here, for the options:
  //   Vc or no Vc

  //
  // targets are field points, with no core radius ===============================================
  //
    std::cout << "    0v_0p compute influence of" << src.to_string() << " on" << targ.to_string() << std::endl;
    // targets are field points

#ifdef USE_VC
    if (env.get_instrs() == cpu_vc) {

      // define vector types for Vc
      typedef Vc::Vector<S> StoreVec;
      typedef Vc::SimdArray<A, Vc::Vector<S>::size()> AccumVec;

      // initialize float_v versions of the source vectors
      const Vc::Memory<StoreVec> sxv  = stdvec_to_vcvec<S>(sx[0], 0.0);
      const Vc::Memory<StoreVec> syv  = stdvec_to_vcvec<S>(sx[1], 0.0);
      const Vc::Memory<StoreVec> szv  = stdvec_to_vcvec<S>(sx[2], 0.0);
      const Vc::Memory<StoreVec> srv  = stdvec_to_vcvec<S>(sr,    1.0);
      const Vc::Memory<StoreVec> ssxv = stdvec_to_vcvec<S>(ss[0], 0.0);
      const Vc::Memory<StoreVec> ssyv = stdvec_to_vcvec<S>(ss[1], 0.0);
      const Vc::Memory<StoreVec> sszv = stdvec_to_vcvec<S>(ss[2], 0.0);

        #pragma omp parallel for
        for (int32_t i=0; i<(int32_t)targ.get_n(); ++i) {
          const StoreVec txv = tx[0][i];
          const StoreVec tyv = tx[1][i];
          const StoreVec tzv = tx[2][i];
          AccumVec accumw = 0.0;
          for (size_t j=0; j<sxv.vectorsCount(); ++j) {
            const AccumVec dx = txv - sxv.vector(j);
            const AccumVec dy = tyv - syv.vector(j);
            const AccumVec dz = tzv - szv.vector(j);
            const AccumVec oor2 = my_recip(srv.vector(j)*srv.vector(j));
            const AccumVec distsq = (dx*dx + dy*dy + dz*dz) * oor2;
            AccumVec toadd = StoreVec(0.0);
            toadd(distsq < StoreVec(16.0)) = StoreVec(2.0) * ssxv.vector(j) * my_exp(-distsq) * oor2;
            accumw += toadd;
          }
          tw[0][i] += accumw.sum();
        }
        //std::cout << "pt " << i << " has new vel " << tu[0][i] << " " << tu[1][i] << std::endl;
        flops *= 1.0 + 13.0 * (float)src.get_n();
    } else
#endif  // no Vc
    {
        #pragma omp parallel for
        for (int32_t i=0; i<(int32_t)targ.get_n(); ++i) {
          A accumw = 0.0;
          for (size_t j=0; j<src.get_n(); ++j) {
            const A dx = tx[0][i] - sx[0][j];
            const A dy = tx[1][i] - sx[1][j];
            const A dz = tx[2][i] - sx[2][j];
            const A distsq = (dx*dx + dy*dy + dz*dz) / std::pow(sr[j], 2);
            // use the Gaussian function
            accumw += S(2.0) * ss[0][j] * std::exp(-distsq) / std::pow(sr[j], 2);
          }
          tw[0][i] += accumw;
          //std::cout << "pt " << i << " at " << tx[0][i] << " " << tx[1][i] << " has vort " << tw[i] << std::endl;
        }
        flops *= 1.0 + 13.0 * (float)src.get_n();
    }

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  const float gflops = 1.e-9 * flops / (float)elapsed_seconds.count();
  printf("    points_affect_points_vorticity: [%.4f] seconds at %.3f GFlop/s\n", (float)elapsed_seconds.count(), gflops);
}


//
// Do the same calculation but save the matrix - use this to solve the equation
//
template <class S>
Vector<S> points_on_points_vort_coeff (Points<S> const& src, Points<S>& targ, const ExecEnv& env) {

  assert(false && "InfluenceVort::points_on_points_vort_coeff not implemented");

}
