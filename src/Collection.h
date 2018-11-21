#pragma once

#include "Points.h"
#include "Panels.h"

#include <variant>

// alias for any type of collection of elements
using Collection = std::variant<Points<float>,
				Panels<float>>;

