/*
 * IglDecimate.h - Call into igl to decimate a triangle mesh
 *
 * (c)2021 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "ElementPacket.h"

void decimate_geometry(ElementPacket<float>&, const size_t);
