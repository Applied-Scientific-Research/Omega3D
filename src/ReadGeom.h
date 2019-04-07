/*
 * ReadGeom.h - Call into igl to read a triangle mesh
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega3D.h"

#include <string>

ElementPacket<float> read_geometry_file(const std::string);
