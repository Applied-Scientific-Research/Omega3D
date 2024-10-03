#
# Omega3D
# (c)2017-24 Applied Scientific Research, Inc.
#
CMAKE_MINIMUM_REQUIRED( VERSION 3.14 )
PROJECT( Omega3D )

include(FetchContent)

# User-serviceable variables here
SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE "Release")
ENDIF ()
SET (CMAKE_BUILD_TYPE "Release" CACHE STRING "Select which configuration to build")
SET (OMEGA3D_BUILD_GUI TRUE CACHE BOOL "Build the GUI version")
SET (OMEGA3D_BUILD_BATCH FALSE CACHE BOOL "Build the batch (no GUI) version")
SET (CMAKE_INSTALL_PREFIX CACHE PATH "Installation location for binaries, sample inputs, and licenses")
SET (OMEGA3D_USE_OMP TRUE CACHE BOOL "Use OpenMP multithreading")
SET (OMEGA3D_USE_VC TRUE CACHE BOOL "Use Vc for vector arithmetic")
SET (OMEGA3D_USE_OGL_COMPUTE FALSE CACHE BOOL "Use OpenGL compute shaders")
SET (OMEGA3D_USE_PLUGIN_AVRM FALSE CACHE BOOL "Enable adaptive VRM plugin")
SET (OMEGA3D_USE_PLUGIN_SIMPLEX FALSE CACHE BOOL "Enable simplex solver plugin for VRM")
SET (OMEGA3D_USE_ONBODY FALSE CACHE BOOL "Enable onbody velocity solver")
#SET (OMEGA3D_USE_HO_HYBRID FALSE CACHE BOOL "Enable external hybrid solver, if any")
SET (OMEGA3D_DO_PROFILE FALSE CACHE BOOL "Enable profiling via gprof")
SET (CMAKE_VERBOSE_MAKEFILE on)
SET (CMAKE_EXPORT_COMPILE_COMMANDS on)

# Output more info if debuging
IF (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  ADD_DEFINITIONS(-DVERBOSE=true)
ELSE ()
  ADD_DEFINITIONS(-DVERBOSE=false)
ENDIF ()

IF (APPLE)
  SET (CMAKE_INSTALL_PREFIX /usr/local/share)
ENDIF ()

# Finds directory to pre-made sims
IF (${CMAKE_BUILD_TYPE} STREQUAL "Release")
  ADD_DEFINITIONS(-DEXAMPLES_DIR="${CMAKE_INSTALL_PREFIX}/3Dexamples")
ELSE ()
  ADD_DEFINITIONS(-DEXAMPLES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/3Dexamples")
ENDIF ()


# Define compiler-specific options

IF (CMAKE_COMPILER_IS_GNUCXX)
  # support Core2 at a minimum, haswell ideally
  #SET (CMAKE_CXX_FLAGS "-Wall -Wformat -Wno-int-in-bool-context -std=c++17 -march=core2 -mtune=haswell")
  SET (CMAKE_CXX_FLAGS "-Wall -Wformat -std=c++17 -Wno-int-in-bool-context")
  #SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb3 -fsanitize=thread")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb3 -fsanitize=address")
  IF( NOT OMEGA3D_USE_VC )
    SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")
  ENDIF()
  SET (CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")
  SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -march=native -g -ggdb3")
  IF (OMEGA3D_DO_PROFILE)
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
  ENDIF ()
ELSEIF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  SET (CMAKE_CXX_FLAGS "-Wall -Wformat -std=c++17 -stdlib=libc++")
  SET (CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb")
  IF( NOT OMEGA3D_USE_VC )
    SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -D_GLIBCXX_DEBUG")
  ENDIF()
  SET (CMAKE_CXX_FLAGS_RELEASE "-O3")
  SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -ggdb")
ELSEIF (MSVC)
  SET (CMAKE_CXX_FLAGS "/std:c++17 /EHsc /D_USE_MATH_DEFINES /DNOMINMAX")
  SET (CMAKE_CXX_FLAGS_DEBUG "/Zi")
  SET (CMAKE_CXX_FLAGS_RELEASE "/O2")
ENDIF ()


# Fetch 3rd party libraries (Vc, onbody, etc,)
SET( BASE_LIBS ${CMAKE_THREAD_LIBS_INIT} )
SET( EXTERNAL_LIBS "" )

# mandatory libraries
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
SET( BASE_LIBS ${BASE_LIBS} Eigen3::Eigen )

# optional libraries
FetchContent_Declare(
  Vc
  URL https://github.com/VcDevel/Vc/archive/refs/tags/1.4.5.zip
  DOWNLOAD_EXTRACT_TIMESTAMP true
)
IF( OMEGA3D_USE_VC )
  FetchContent_MakeAvailable(Vc)
  SET (BASE_LIBS ${BASE_LIBS} Vc)
  SET (CPREPROCDEFS ${CPREPROCDEFS} -DUSE_VC)
ENDIF ()

# OpenMP for multithreading
IF( OMEGA3D_USE_OMP )
  FIND_PACKAGE (OpenMP REQUIRED)
  IF (OPENMP_FOUND)
    SET( BASE_LIBS ${BASE_LIBS} OpenMP::OpenMP_CXX )
  ENDIF ()
ENDIF()

FetchContent_Declare(
  onbody
  GIT_REPOSITORY https://github.com/Applied-Scientific-Research/onbody.git
  GIT_TAG 37204f21210321a59a69bdf63c88062ed862c96c
)
IF( OMEGA3D_USE_ONBODY )
  FetchContent_MakeAvailable(onbody)
  SET (EXTERNAL_LIBS ${EXTERNAL_LIBS} bh3dvortgrads)
  SET (CPREPROCDEFS ${CPREPROCDEFS} -DEXTERNAL_VEL_SOLVE)
ENDIF ()


# Define platform-specific libraries and options

IF( WIN32 )
  # because clang is stupid
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  ENDIF ()

  SET( GL_LIBS_HOME "C:/Lib" CACHE PATH "Location of support libraries" )
  SET( GLFW3_ROOT "${GL_LIBS_HOME}/glfw-3.3.bin.WIN32" CACHE STRING "GLFW3 root" )
  INCLUDE_DIRECTORIES( "${GLFW3_ROOT}/include" )
  LINK_DIRECTORIES   ( "${GLFW3_ROOT}/lib-vc2017" )
  SET( FRAMEWORK_LIBS glfw3.lib )

  SET( GUI_LIBS ${FRAMEWORK_LIBS} opengl32 )
  SET_PROPERTY(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Omega3D)

ELSEIF( APPLE )
  # because clang is stupid
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  ENDIF ()

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

  # turn OFF compute shaders because Apple doesn't support them
  SET( OMEGA3D_USE_OGL_COMPUTE FALSE )

  SET( CMAKE_MACOSX_RPATH 1 )
  SET( GUI_LIBS ${FRAMEWORK_LIBS} ${COCOA_LIBRARY} ${IOKIT_LIBRARY} ${COREVIDEO_LIBRARY} ${OpenGL_LIBRARY} )

ELSEIF( UNIX )
  # because clang is stupid
  IF (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
  ENDIF ()

  # glfw finds opengl
  # this works for CentOS 7
  #FIND_PACKAGE( glfw3 3.0 REQUIRED )
  # this works for Fedora 21
  #FIND_PACKAGE( glfw 3.0 REQUIRED )
  # this seems to be more general
  FIND_LIBRARY( FRAMEWORK_LIBS NAMES glfw3 glfw libglfw PATHS /usr/local/lib )

  # tell imgui to use glad
  SET (CPREPROCDEFS ${CPREPROCDEFS} -DIMGUI_IMPL_OPENGL_LOADER_GLAD)

  # but glad needs dl, which cmake doesn't think glfw needs
  SET( FRAMEWORK_LIBS glfw dl )

  # OpenMP for multithreading
  IF( NOT OMEGA3D_USE_OMP )
    # suppress the warning
    SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    # need to find pthreads (for std::async?)
    FIND_PACKAGE( Threads REQUIRED )
  ENDIF()

  SET( GUI_LIBS ${FRAMEWORK_LIBS} ${GLFW_LIBRARIES} )
ENDIF()


# Enable plugins

# adaptive VRM
IF( OMEGA3D_USE_PLUGIN_AVRM )
  SET (CPREPROCDEFS ${CPREPROCDEFS} -DPLUGIN_AVRM)
  INCLUDE_DIRECTORIES( "plugin-o3d-avrm" )
  # add VRMadaptive.h to HEADERS?
ENDIF()

# simplex VRM solver
IF( OMEGA3D_USE_PLUGIN_SIMPLEX )
  SET (CPREPROCDEFS ${CPREPROCDEFS} -DPLUGIN_SIMPLEX)
  INCLUDE_DIRECTORIES( "plugin-o2d-simplex" )
  # add simplex.h to HEADERS?
ENDIF()

# external fast summation routine
#IF( OMEGA3D_USE_EXTERNAL_SUM )
#  SET (CPREPROCDEFS ${CPREPROCDEFS} -DEXTERNAL_VEL_SOLVE)
#  SET (FASTSUM_LIBS "libfastsumm.a" CACHE STRING "Static library exposing velocity summation function call" )
#  #FIND_LIBRARY( FASTSUM_LIBS NAMES Vc PATHS "${VC_ROOT}/lib" )
#  SET( EXTERNAL_LIBS ${FASTSUM_LIBS} )
#ELSE()
#  SET (FASTSUM_LIBS "")
#  SET( EXTERNAL_LIBS ${FASTSUM_LIBS} )
#ENDIF()

# Enable submodules
#IF( OMEGA3D_USE_HO_HYBRID )
  #ADD_SUBDIRECTORY( extern/HO-3D )
  #INCLUDE_DIRECTORIES( "extern/HO-3D/src" )
  #INCLUDE_DIRECTORIES ( "${EIGEN_ROOT}" )
  #INCLUDE_DIRECTORIES ( "${AMGCL_ROOT}" )
  #SET( EXTERNAL_LIBS ${EXTERNAL_LIBS} HO-3D )
  #SET (CPREPROCDEFS ${CPREPROCDEFS} -DHO-3D)
#ENDIF()

# Define files and targets

ADD_DEFINITIONS (${CPREPROCDEFS})

INCLUDE_DIRECTORIES( "lib" )
# For some of the libraries who need glad
INCLUDE_DIRECTORIES( "lib/glad" )
# For some of the internal libraries who need nanoflann
INCLUDE_DIRECTORIES( "lib/nanoflann" )
# For some of the internal libraries who need nnls 
INCLUDE_DIRECTORIES( "lib/eigen-nnls" )
INCLUDE_DIRECTORIES( "src" )

IF (CMAKE_COMPILER_IS_GNUCXX)
  SET_SOURCE_FILES_PROPERTIES("lib/imgui/imgui_draw.cpp" PROPERTIES COMPILE_FLAGS -Wno-maybe-uninitialized)
ENDIF()

FILE(GLOB IMGUISOURCES "lib/imgui/*.cpp" )
SET(SOURCES "src/Simulation.cpp"
            "src/FlowFeature.cpp"
            "src/CollectionHelper.cpp"
            "src/BoundaryFeature.cpp"
            "src/MeasureFeature.cpp"
            "src/Body.cpp"
            "src/RenderParams.cpp"
            "src/JsonHelper.cpp"
            "src/GeomHelper.cpp"
            "src/IglReadGeom.cpp"
            "src/IglRefine.cpp"
            "src/IglDecimate.cpp"
            "src/IglMergeDups.cpp"
            "src/StatusFile.cpp"
            "lib/tinyxml2/tinyxml2.cpp"
            "lib/tinyexpr/tinyexpr.c" )
SET(GUI_SOURCES "lib/glad/glad.c"
                "src/ShaderHelper.cpp"
                "src/FeatureDraw.cpp"
                "lib/miniz/miniz.c"
                "lib/miniz/FrameBufferToImage.cpp" )

# create a binary for the GUI version
IF( OMEGA3D_BUILD_GUI )
  ADD_EXECUTABLE( ${PROJECT_NAME} ${IMGUISOURCES} ${SOURCES} ${GUI_SOURCES} "src/main_gui.cpp" )
  SET_TARGET_PROPERTIES( ${PROJECT_NAME} PROPERTIES OUTPUT_NAME "${PROJECT_NAME}.bin" )
  TARGET_LINK_LIBRARIES( ${PROJECT_NAME} ${BASE_LIBS} ${GUI_LIBS} ${EXTERNAL_LIBS} )
  TARGET_COMPILE_DEFINITIONS( ${PROJECT_NAME} PRIVATE "-DUSE_GL" )
  TARGET_COMPILE_DEFINITIONS( ${PROJECT_NAME} PRIVATE "-DUSE_IMGUI" )
  IF( OMEGA3D_USE_OGL_COMPUTE AND NOT APPLE)
    TARGET_COMPILE_DEFINITIONS( ${PROJECT_NAME} PRIVATE "-DUSE_OGL_COMPUTE" )
  ENDIF()
  INSTALL( TARGETS ${PROJECT_NAME} DESTINATION bin )
ENDIF()

# create a binary for the command-line version
IF( OMEGA3D_BUILD_BATCH )
  ADD_EXECUTABLE( "${PROJECT_NAME}batch" ${SOURCES} "src/main_batch.cpp" )
  SET_TARGET_PROPERTIES( "${PROJECT_NAME}batch" PROPERTIES OUTPUT_NAME "${PROJECT_NAME}batch.bin" )
  TARGET_LINK_LIBRARIES( "${PROJECT_NAME}batch" ${BASE_LIBS} ${EXTERNAL_LIBS} )
  INSTALL( TARGETS "${PROJECT_NAME}batch" DESTINATION bin )
ENDIF()

INSTALL( DIRECTORY 3Dexamples/ DESTINATION 3Dexamples )
INSTALL( FILES LICENSE DESTINATION LICENSE )

