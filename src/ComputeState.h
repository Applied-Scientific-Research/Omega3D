/*
 * ComputeState.h - An atomic enum to allow interthread coordination
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

// solver order and method
enum compute_state_t {
  no_compute    = 0,
  begin_compute = 1,
  computing     = 2,
  compute_done  = 3
};

