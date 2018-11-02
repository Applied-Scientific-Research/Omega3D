#pragma once
#include <iostream>
#include <variant>
#include "Points.h"
#include "Panels.h"


// helper alias for any type of collection of elements
using Collection = std::variant<Points<float>,
				Panels<float>>;

// helper function
std::string to_string(Collection& _c) {
  return std::visit([=](auto& elem) { return elem.to_string(); }, _c);
}

