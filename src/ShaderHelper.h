//
// ShaderHelper.h
//
// These methods are from https://solarianprogrammer.com/2013/05/13/opengl-101-drawing-primitives/
//  https://github.com/sol-prog/OpenGL-101
//

#pragma once

// for glad
#ifdef _WIN32
  #define APIENTRY __stdcall
#endif
#include "glad.h"

// Create a program from two shaders
GLuint create_particle_program();
//GLuint create_particlept_program();
//GLuint create_panel_program();

// Useful defines for templated OpenGL calls
template <class T> GLenum get_gl_index_type() = delete;
template <> inline GLenum get_gl_index_type<uint8_t>() { return GL_UNSIGNED_BYTE; }
template <> inline GLenum get_gl_index_type<uint16_t>() { return GL_UNSIGNED_SHORT; }
template <> inline GLenum get_gl_index_type<uint32_t>() { return GL_UNSIGNED_INT; }

template <class T> GLenum get_gl_float_type() = delete;
template <> inline GLenum get_gl_float_type<float>() { return GL_FLOAT; }
template <> inline GLenum get_gl_float_type<double>() { return GL_DOUBLE; }
