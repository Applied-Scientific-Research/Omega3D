# Omega3D
High-performance 3D computational fluid dynamics solver with easy GUI

## Overview

Omega3D aims to be an accurate combined Lagrangian-Eulerian fluid flow solver for unsteady flows with complex boundaries, with a greatly reduced reliance on meshing and parameter tuning for achieving accuracy and speed.

## Build and run
This code uses some C++17 features, so should compile on GCC 7, Clang 4, and MSVC 19.10 compilers.

#### Prerequisites
Both the GUI and batch versions require CMake to compile.  To build the GUI version, users will also need GLFW3. These can be installed in Red Hat/Fedora with

    sudo yum install cmake glfw3-devel

or on Ubuntu with

    sudo apt-get install cmake glfw3-dev

or on OSX with

    brew install cmake glfw3

#### Optional libraries
[Vc](https://github.com/VcDevel/Vc) is a vectorization library, and Omega3D uses it to greatly accelerate the velocity evaluations. This package can be built and installed external to Omega3D with

    git clone https://github.com/VcDevel/Vc.git
    cd Vc
    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/opt/Vc -DBUILD_TESTING=OFF ..
    make -j 4
    sudo make install
    cd ../..

#### Compilation
Upon installation of the prerequisites, the following commands should build Omega3D.

    git clone git@github.com:Applied-Scientific-Research/Omega3D.git
    cd Omega3D
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=ON -DUSE_VC=OFF ..
    make

If you were able to build and install Vc, then you should set `-DUSE_VC=ON` in the above `cmake` command.

Finally, you can run the GUI program with

    ./Omega3D.bin

and the batch program with

    ./Omega3Dbatch.bin

## To do
Tasks to consider or implement:

* Add second order convection - in a class?
* Use the actual core function to draw the blobs - but what is the real core function?
* Add arcball rotation to the viewport - make my brain hurt
* Use [libigl](https://github.com/libigl/libigl/) or [OpenMesh](http://openmesh.org/intro/) to load geometry files for boundaries
* Start fresh GUI main file, look for first run and splash a help window
* Add other repos as submodules, like [Vc](https://github.com/VcDevel/Vc) and [nlohmann/json](https://github.com/nlohmann/json) and [libigl](https://github.com/libigl/libigl/), or just by copying? `submodule add https://...xxx.git thirdparty/xxx`
* ~~Rework the particle shader to read from all the various arrays - goal is to draw an inviscid sim~~
* ~~Rework the updateGL code to move all the appropriate arrays to the GPU~~
* ~~Make the reset button actally reset~~
* ~~Get new features to just add particles to the current Points object~~
* ~~Have the GUI set up a vortex ring object~~
* ~~Get it to run a vortex ring without drawing anything~~
* ~~Copy the simulation time step logic from Omega2D (using std::async)~~
* ~~Move time step code into Simulation/Convection files, like Omega2D~~
* ~~Make separate batch and GUI main files and binaries~~
* ~~Support CMake with optional OpenMP and optional Vc~~
* ~~Wrap solver in a GUI - use existing code from [Omega2D](https://github.com/Applied-Scientific-Research/Omega2D). This means supporting only, say, vortex rings and stretch, but no diffusion.~~
* ~~Rework into struct of arrays to allow Vc to work properly (this requires more work for compute and draw shaders :one buffer per array), but should be nice and flexible - and easier than all the `4*idx+3` crap~~
* ~~Test speed and accuracy with Vc as the float type~~
* ~~move rand() to std::random (see Omega2D)~~

## Thanks

Many thanks to NBL for valuable discussions of architecture and C++ syntax and idioms.
