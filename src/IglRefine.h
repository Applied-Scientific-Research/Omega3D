/*
 * IglRefine.h - Call into igl to refine a triangle mesh
 *
 * (c)2019-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "ElementPacket.h"

void refine_geometry(ElementPacket<float>&);
