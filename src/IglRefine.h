/*
 * IglRefine.h - Call into igl to refine a triangle mesh
 *
 * (c)2019 Applied Scientific Research, Inc.
 *         Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "Omega3D.h"

void refine_geometry(ElementPacket<float>&);
