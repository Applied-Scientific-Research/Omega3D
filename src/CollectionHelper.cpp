//#include "CollectionHelper.h"
#include "Collection.h"

#include <iostream>
#include <variant>

// helper function
std::string to_string(Collection& _c) {
  return std::visit([=](auto& elem) { return elem.to_string(); }, _c);
}

