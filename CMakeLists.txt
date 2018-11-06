CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
PROJECT( Omega3D )

set(CMAKE_VERBOSE_MAKEFILE on)

set(CMAKE_EXPORT_COMPILE_COMMANDS on)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

IF (CMAKE_COMPILER_IS_GNUCXX)
  SET (CMAKE_CXX_FLAGS "-Wall -Wformat -std=c++17")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb")
  SET (CMAKE_CXX_FLAGS_RELEASE "-O3 -mavx2")
ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET (CMAKE_CXX_FLAGS "-O2 -Wall -Wformat -std=c++17 -stdlib=libc++")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb")
ELSEIF (MSVC)
  SET (CMAKE_CXX_FLAGS "/Ox")
ENDIF ()

IF( WIN32 )
  SET( GL_LIBS_HOME "C:/lib" CACHE PATH "Location of support libraries" )

  SET( GLFW3_ROOT "${GL_LIBS_HOME}/glfw-3.1.2.bin.WIN32" CACHE STRING "GLFW3 root" )

  INCLUDE_DIRECTORIES( "${GLFW3_ROOT}/include" )
  LINK_DIRECTORIES   ( "${GLFW3_ROOT}/lib-vc2015" )
  SET( FRAMEWORK_LIBS glfw3.lib )

  SET( PLATFORM_LIBS ${FRAMEWORK_LIBS} opengl32 )

ELSEIF( APPLE )
  INCLUDE_DIRECTORIES ( /System/Library/Frameworks )
  FIND_LIBRARY(COCOA_LIBRARY Cocoa)
  FIND_LIBRARY(IOKIT_LIBRARY IOKit)
  FIND_LIBRARY(COREVIDEO_LIBRARY CoreVideo)
  FIND_LIBRARY(OpenGL_LIBRARY OpenGL )
  MARK_AS_ADVANCED (COCOA_LIBRARY
                    IOKIT_LIBRARY
                    OpenGL_LIBRARY)

  #SET( FRAMEWORK_LIBS glfw3 )
  FIND_LIBRARY( FRAMEWORK_LIBS NAMES glfw3 glfw libglfw PATHS /usr/local/lib)

  INCLUDE_DIRECTORIES ( "/usr/local/Cellar/eigen/3.3.4/include/eigen3" )

  SET(PLATFORM_LIBS ${FRAMEWORK_LIBS} ${COCOA_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY} ${OpenGL_LIBRARY})

  FIND_PACKAGE( OpenMP )
  if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()


ELSEIF( UNIX )
  # glfw finds opengl
  # this works for CentOS 7
  #FIND_PACKAGE( glfw3 3.0 REQUIRED )
  # this works for Fedora 21
  #FIND_PACKAGE( glfw 3.0 REQUIRED )
  # this seems to be more general
  FIND_LIBRARY( FRAMEWORK_LIBS NAMES glfw3 glfw libglfw PATHS /usr/local/lib)

  # but glad needs dl, which cmake doesn't think glfw needs
  SET( FRAMEWORK_LIBS glfw dl )
  SET( PLATFORM_LIBS ${FRAMEWORK_LIBS} ${GLFW_LIBRARIES})

  INCLUDE_DIRECTORIES ( "/usr/include/eigen3" )

  FIND_PACKAGE( OpenMP REQUIRED)
  if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()

ENDIF()

set(GLADSRC "src/glad/glad.c" )
INCLUDE_DIRECTORIES( "src/glad" )
file(GLOB SOURCES "src/*.cpp" )
file(GLOB IMGUISOURCES "src/imgui/*.cpp" )

ADD_EXECUTABLE( ${PROJECT_NAME} ${GLADSRC} ${IMGUISOURCES} ${SOURCES} )
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${PROJECT_NAME}.bin")
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${PLATFORM_LIBS} )

