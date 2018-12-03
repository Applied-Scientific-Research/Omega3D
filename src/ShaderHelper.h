/*
 * ShaderHelper.h - Methods for generating opengl programs
 *
 * These methods are from https://solarianprogrammer.com/2013/05/13/opengl-101-drawing-primitives/
 * and https://github.com/sol-prog/OpenGL-101
 *
 * (c)2017-8 Applied Scientific Research, Inc.
 *           Written by Mark J Stock <markjstock@gmail.com>
 */

#pragma once

// for glad
#ifdef _WIN32
  #define APIENTRY __stdcall
#endif
#include "glad.h"

#include <cstdint>

// Create a program from two shaders
GLuint create_particle_program();
//GLuint create_particlept_program();
//GLuint create_panel_program();

// Useful defines for templated OpenGL calls

template <class T>
struct GetTypeHelper {
};

template <> struct GetTypeHelper<uint8_t> { enum { type = GL_UNSIGNED_BYTE }; };
template <> struct GetTypeHelper<uint16_t> { enum { type = GL_UNSIGNED_SHORT }; };
template <> struct GetTypeHelper<uint32_t> { enum { type = GL_UNSIGNED_INT }; };
template <> struct GetTypeHelper<float> { enum { type = GL_FLOAT }; };
template <> struct GetTypeHelper<double> { enum { type = GL_DOUBLE }; };

template <class T>
constexpr auto get_gl_type = GetTypeHelper<T>::type;

