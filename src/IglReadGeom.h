/*
 * IglReadGeom.h - Call into igl to read a triangle mesh
 *
 * (c)2019-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "ElementPacket.h"

#include <string>

ElementPacket<float> read_geometry_file(const std::string);
