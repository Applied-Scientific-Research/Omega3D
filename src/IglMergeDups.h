/*
 * IglMergeDups.h - Call into igl to merge duplicate verts
 *
 * (c)2021 Applied Scientific Research, Inc.
 *         Mark J Stock <markjstock@gmail.com>
 */

#pragma once

#include "ElementPacket.h"

void mergeduplicates_geometry(ElementPacket<float>&, const float);
