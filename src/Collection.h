/*
 * Collection.h - definition of variant type for element collections
 *
 * (c)2018-9 Applied Scientific Research, Inc.
 *           Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Points.h"
#include "Surfaces.h"

#include <variant>

// alias for any type of collection of elements
// eventually will have Lines<float> and Volumes<float> here

using Collection = std::variant<Points<float>,
                                Surfaces<float>>;

