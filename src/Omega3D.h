#pragma once

const size_t Dimensions = 3;

// templated solution type

// templated on the accumulator precision
//template <class A>
//class Solver {
//public:
//private:
//};

enum solution_t {
  direct_cpu   = 1,
  direct_vc    = 2,
  direct_glsl  = 3,
  treecode_cpu = 4
};

enum elem_t {
  active   = 1,
  reactive = 2,
  inert    = 3
};

enum move_t {
  lagrangian = 1,
  bodybound  = 2,
  fixed      = 3
};

