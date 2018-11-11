#pragma once

#define _USE_MATH_DEFINES
#include <cmath>


// velocity influence functions

template <class S, class A>
static inline void kernel_0_0s (const S sx, const S sy, const S sz,
                                const S sr,
                                const S ssx, const S ssy, const S ssz,
                                const S tx, const S ty, const S tz,
                                const S tr,
                                A* __restrict__ tu) {
  // 30 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr + tr*tr;
  r2 = 1.0 / (r2*std::sqrt(r2));
  const S dxxw = dz*ssy - dy*ssz;
  const S dyxw = dx*ssz - dz*ssx;
  const S dzxw = dy*ssx - dx*ssy;
  tu[0] += r2 * dxxw;
  tu[1] += r2 * dyxw;
  tu[2] += r2 * dzxw;
}

template <class S, class A>
static inline void kernel_0_0v (const S* __restrict__ sx, const S __restrict__ sr, const S* __restrict__ ss,
                                const S* __restrict__ tx, const S __restrict__ tr, A* __restrict__ tu) {
  // 30 flops
  const S dx = tx[0] - sx[0];
  const S dy = tx[1] - sx[1];
  const S dz = tx[2] - sx[2];
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr + tr*tr;
  r2 = 1.0 / (r2*std::sqrt(r2));
  const S dxxw = dz*ss[1] - dy*ss[2];
  const S dyxw = dx*ss[2] - dz*ss[0];
  const S dzxw = dy*ss[0] - dx*ss[1];
  tu[0] += r2 * dxxw;
  tu[1] += r2 * dyxw;
  tu[2] += r2 * dzxw;
}

template <class S, class A>
static inline void kernel_0_0sg (const S sx, const S sy, const S sz,
                                 const S sr,
                                 const S ssx, const S ssy, const S ssz,
                                 const S tx, const S ty, const S tz,
                                 const S tr,
                                 A* __restrict__ tu) {
  // 30 flops
  const S dx = tx - sx;
  const S dy = ty - sy;
  const S dz = tz - sz;
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr + tr*tr;
  r2 = 1.0 / (r2*std::sqrt(r2));
  S dxxw = dz*ssy - dy*ssz;
  S dyxw = dx*ssz - dz*ssx;
  S dzxw = dy*ssx - dx*ssy;
  tu[0] += r2 * dxxw;
  tu[1] += r2 * dyxw;
  tu[2] += r2 * dzxw;

  // HACK - you need to figure out what this term is
  const S bbb = r2 / std::sqrt(r2);
  // continuing with grads - this section is 33 flops
  dxxw *= bbb;
  dyxw *= bbb;
  dzxw *= bbb;
  tu[3] += dx*dxxw;
  tu[4] += dx*dyxw + ssz*r2;
  tu[5] += dx*dzxw - ssy*r2;
  tu[6] += dy*dxxw - ssz*r2;
  tu[7] += dy*dyxw;
  tu[8] += dy*dzxw + ssx*r2;
  tu[9] += dz*dxxw + ssy*r2;
  tu[10] += dz*dyxw - ssx*r2;
  tu[11] += dz*dzxw;
}

template <class S, class A>
static inline void kernel_0_0vg (const S* __restrict__ sx, const S __restrict__ sr, const S* __restrict__ ss,
                                 const S* __restrict__ tx, const S __restrict__ tr, A* __restrict__ tu) {
  // 30 flops
  const S dx = tx[0] - sx[0];
  const S dy = tx[1] - sx[1];
  const S dz = tx[2] - sx[2];
  S r2 = dx*dx + dy*dy + dz*dz + sr*sr + tr*tr;
  r2 = 1.0 / (r2*std::sqrt(r2));
  S dxxw = dz*ss[1] - dy*ss[2];
  S dyxw = dx*ss[2] - dz*ss[0];
  S dzxw = dy*ss[0] - dx*ss[1];
  tu[0] += r2 * dxxw;
  tu[1] += r2 * dyxw;
  tu[2] += r2 * dzxw;

  // HACK - you need to figure out what this term is
  const S bbb = r2 / std::sqrt(r2);
  // continuing with grads - this section is 33 flops
  dxxw *= bbb;
  dyxw *= bbb;
  dzxw *= bbb;
  tu[3] += dx*dxxw;
  tu[4] += dx*dyxw + ss[2]*r2;
  tu[5] += dx*dzxw - ss[1]*r2;
  tu[6] += dy*dxxw - ss[2]*r2;
  tu[7] += dy*dyxw;
  tu[8] += dy*dzxw + ss[0]*r2;
  tu[9] += dz*dxxw + ss[1]*r2;
  tu[10] += dz*dyxw - ss[0]*r2;
  tu[11] += dz*dzxw;
}


//
// analytic influence of 2d linear constant-strength vortex panel on target point
//   ignoring the 1/2pi factor, which will be multiplied later
//   40 flops average
//
template <class S, class A>
static inline void kernel_1_0v (const S* __restrict__ sx0, const S* __restrict__ sx1, const S str,
                                const S* __restrict__ tx, A* __restrict__ tu) {

  // side lengths of the triangle s0, s1, t
  const S rij2  = std::pow(tx[0]-sx0[0],2) + std::pow(tx[1]-sx0[1],2);
  const S rij   = std::sqrt(rij2);
  const S rij12 = std::pow(tx[0]-sx1[0],2) + std::pow(tx[1]-sx1[1],2);
  const S rij1  = std::sqrt(rij12);
  //std::cout << "rij is " << rij << " and rijp1 is " << rij1 << std::endl;
  const S vstar = std::log(rij/rij1);
  S ustar = std::atan2(tx[0]-sx1[0], tx[1]-sx1[1]) - std::atan2(tx[0]-sx0[0], tx[1]-sx0[1]);
  //std::cout << "ustar started off as " << ustar << std::endl;
  if (ustar < -M_PI) ustar += 2.*M_PI;
  if (ustar > M_PI) ustar -= 2.*M_PI;
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

  const S px    = sx1[0]-sx0[0];
  const S py    = sx1[1]-sx0[1];
  //std::cout << "px is " << px << " and py is " << py << std::endl;

  // finally, rotate back into global coordinates
  const S velx  = ustar*px - vstar*py;
  const S vely  = ustar*py + vstar*px;
  //std::cout << "velx is " << velx << " and vely is " << vely << std::endl;
  const S mult  = str / std::sqrt(std::pow(px,2) + std::pow(py,2));
  //std::cout << "finalx is " << (mult*velx) << " and finaly is " << (mult*vely) << std::endl;

  // and multiply by vortex sheet strength
  tu[0] += mult*velx;
  tu[1] += mult*vely;
}

template <class S, class A>
static inline void kernel_1_0s (const S sx0, const S sy0, const S sz0,
                                const S sx1, const S sy1, const S sz1,
                                const S ssx, const S ssy, const S ssz,
                                const S tx, const S ty, const S tz,
                                A* __restrict__ tu) {

  // side lengths of the triangle s0, s1, t
  const S rij2  = std::pow(tx-sx0,2) + std::pow(ty-sy0,2);
  const S rij   = std::sqrt(rij2);
  const S rij12 = std::pow(tx-sx1,2) + std::pow(ty-sy1,2);
  const S rij1  = std::sqrt(rij12);
  //std::cout << "rij is " << rij << " and rijp1 is " << rij1 << std::endl;
  const S vstar = std::log(rij/rij1);
  S ustar = std::atan2(tx-sx1, ty-sy1) - std::atan2(tx-sx0, ty-sy0);
  //std::cout << "ustar started off as " << ustar << std::endl;
  if (ustar < -M_PI) ustar += 2.*M_PI;
  if (ustar > M_PI) ustar -= 2.*M_PI;
  //std::cout << "ustar is " << ustar << " and vstar is " << vstar << std::endl;

  const S px    = sx1-sx0;
  const S py    = sy1-sy0;
  //std::cout << "px is " << px << " and py is " << py << std::endl;

  // finally, rotate back into global coordinates
  const S velx  = ustar*px - vstar*py;
  const S vely  = ustar*py + vstar*px;
  //std::cout << "velx is " << velx << " and vely is " << vely << std::endl;
  const S mult  = ssx / std::sqrt(std::pow(px,2) + std::pow(py,2));
  //std::cout << "finalx is " << (mult*velx) << " and finaly is " << (mult*vely) << std::endl;

  // and multiply by vortex sheet strength
  tu[0] += mult*velx;
  tu[1] += mult*vely;
}

