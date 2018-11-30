//
// ShaderHelper.cpp
//
// These methods are from
// https://solarianprogrammer.com/2013/05/13/opengl-101-drawing-primitives/
// https://github.com/sol-prog/OpenGL-101
//

#include "ShaderHelper.h"

#include <GLFW/glfw3.h>

#include <string>
#include <vector>
#include <iostream>
//#include <fstream>


// clang-format off

const std::string part_vert_shader_source =
#include "shaders/particle.vert"
;
const std::string part_frag_shader_source =
#include "shaders/particle.frag"
;

//const std::string panel_vert_shader_source =
//#include "shaders/panel.vert"
//;
//const std::string panel_frag_shader_source =
//#include "shaders/panel.frag"
//;


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
    std::cerr << "Shader compilation failed with this message:" << std::endl;
    std::vector<char> compilation_log(512);
    glGetShaderInfoLog(shader, compilation_log.size(), nullptr, &compilation_log[0]);
    std::cerr << &compilation_log[0] << std::endl;
    glfwTerminate();
    exit(-1);
  } else {
    std::cerr << "Shader compilation successful" << std::endl;
  }
  return shader;
}

// Create a program from two shaders to render particles as blobs
GLuint create_particle_program() {
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

