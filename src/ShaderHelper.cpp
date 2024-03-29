/*
 * ShaderHelper.cpp - Methods for generating opengl programs
 *
 * These methods are from https://solarianprogrammer.com/2013/05/13/opengl-101-drawing-primitives/
 * and https://github.com/sol-prog/OpenGL-101
 *
 * (c)2017-20 Applied Scientific Research, Inc.
 *            Mark J Stock <markjstock@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "ShaderHelper.h"
#include "CoreFunc.h"

#include <GLFW/glfw3.h>

#include <string>
#include <vector>
#include <iostream>
//#include <fstream>

// clang-format off

const std::string point_vert_shader_source =
#include "shaders/point.vert"
;
const std::string point_frag_shader_source =
#include "shaders/point.frag"
;

const std::string part_vert_shader_source =
#include "shaders/particle.vert"
;
const std::string part_frag_shader_source =
#include "shaders/particle.frag"
;

const std::string panel_vert_shader_source =
#include "shaders/flattriangle.vert"
;
const std::string panel_frag_shader_source =
#include "shaders/flattriangle.frag"
;

//
// Select core function based on cpp defines in CoreFunc.h
//
#ifdef USE_RM_KERNEL
const std::string ptpt_velgrad_shader_source =
#include "shaders/ptptvelgradrm.comp"
;
#endif
#ifdef USE_EXPONENTIAL_KERNEL
const std::string ptpt_velgrad_shader_source =
#include "shaders/ptptvelgradexp.comp"
;
#endif
#ifdef USE_WL_KERNEL
const std::string ptpt_velgrad_shader_source =
#include "shaders/ptptvelgradwl.comp"
;
#endif
#ifdef USE_V2_KERNEL
const std::string ptpt_velgrad_shader_source =
#include "shaders/ptptvelgradv2.comp"
;
#endif
#ifdef USE_V3_KERNEL
const std::string ptpt_velgrad_shader_source =
#include "shaders/ptptvelgradrm.comp"
;
#endif

const std::string init_velgrad_shader_source =
#include "shaders/zeropointvels.comp"
;


// Compile a shader
GLuint load_and_compile_shader(const std::string shader_src, GLenum shaderType) {
  // load the string into the right data type
  const char* as_cstr = shader_src.c_str();
  //std::cout << "Shader source:" << as_cstr;

  // Compile the shader
  GLuint shader = glCreateShader(shaderType);
  glShaderSource(shader, 1, &as_cstr, nullptr);
  glCompileShader(shader);

  // Check the result of the compilation
  GLint test;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &test);
  if (!test) {
    std::cerr << "  Shader compilation failed with this message:" << std::endl;
    std::vector<char> compilation_log(512);
    glGetShaderInfoLog(shader, compilation_log.size(), nullptr, &compilation_log[0]);
    std::cerr << &compilation_log[0] << std::endl;
    glfwTerminate();
    exit(-1);
  } else {
    std::cerr << "  Shader compilation successful" << std::endl;
  }
  return shader;
}

// Create a program from two shaders to render particles as blobs
GLuint create_draw_blob_program() {
  // Load and compile the vertex and fragment shaders
  GLuint vertexShader = load_and_compile_shader(part_vert_shader_source, GL_VERTEX_SHADER);
  GLuint fragmentShader = load_and_compile_shader(part_frag_shader_source, GL_FRAGMENT_SHADER);

  // Attach the above shader to a program
  GLuint shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);

  // Flag the shaders for deletion
  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);

  // Link and use the program
  glLinkProgram(shaderProgram);
  glUseProgram(shaderProgram);

  return shaderProgram;
}


// Create a program from one shader to render particles as points
GLuint create_draw_point_program() {
  // Load and compile the vertex and fragment shaders
  GLuint vertexShader = load_and_compile_shader(point_vert_shader_source, GL_VERTEX_SHADER);
  GLuint fragmentShader = load_and_compile_shader(point_frag_shader_source, GL_FRAGMENT_SHADER);

  // Attach the above shader to a program
  GLuint shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);

  // Flag the shaders for deletion
  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);

  // Link and use the program
  glLinkProgram(shaderProgram);
  glUseProgram(shaderProgram);

  return shaderProgram;
}


// Create a program from two shaders to render panels
GLuint create_draw_surface_tri_prog() {
  // Load and compile the vertex and fragment shaders
  GLuint vertexShader = load_and_compile_shader(panel_vert_shader_source, GL_VERTEX_SHADER);
  GLuint fragmentShader = load_and_compile_shader(panel_frag_shader_source, GL_FRAGMENT_SHADER);

  // Attach the above shader to a program
  GLuint shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);

  // Flag the shaders for deletion
  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);

  // Link and use the program
  glLinkProgram(shaderProgram);
  glUseProgram(shaderProgram);

  return shaderProgram;
}


// Create a program from two shaders to render panels
GLuint create_vertfrag_prog(const std::string vert_shader_src,
                            const std::string frag_shader_src) {

  // Load and compile the vertex and fragment shaders
  GLuint vertexShader = load_and_compile_shader(vert_shader_src, GL_VERTEX_SHADER);
  GLuint fragmentShader = load_and_compile_shader(frag_shader_src, GL_FRAGMENT_SHADER);

  // Attach the above shader to a program
  GLuint shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);

  // Flag the shaders for deletion
  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);

  // Link and use the program
  glLinkProgram(shaderProgram);
  int success;
  char infoLog[512];
  glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
  if (!success) {
    glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
    std::cout << "ERROR: Shader program linking failed\n" << infoLog << std::endl;
  }
  glUseProgram(shaderProgram);

  return shaderProgram;
}

// Create a compute program from one shader
GLuint create_ptptvelgrad_program() {
  // Load and compile the vertex and fragment shaders
  GLuint computeShader = load_and_compile_shader(ptpt_velgrad_shader_source, GL_COMPUTE_SHADER);

  // Attach the above shader to a program
  GLuint shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, computeShader);

  // Flag the shaders for deletion
  glDeleteShader(computeShader);

  // Link and use the program
  glLinkProgram(shaderProgram);
  glUseProgram(shaderProgram);

  return shaderProgram;
}

// And a simple program to zero the vels and vel grads
GLuint create_initvelgrad_program() {
  // Load and compile the vertex and fragment shaders
  GLuint computeShader = load_and_compile_shader(init_velgrad_shader_source, GL_COMPUTE_SHADER);

  // Attach the above shader to a program
  GLuint shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, computeShader);

  // Flag the shaders for deletion
  glDeleteShader(computeShader);

  // Link and use the program
  glLinkProgram(shaderProgram);
  glUseProgram(shaderProgram);

  return shaderProgram;
}
