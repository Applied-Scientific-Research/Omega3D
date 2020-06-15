# Omega3D
3D flow solver with GUI using vortex particle and boundary element methods


## Overview
[Computational Fluid Dynamics (CFD)](https://en.wikipedia.org/wiki/Computational_fluid_dynamics) encompasses a wide variety of methods to aid in the numerical simulation of fluid flows on digital computers. Most methods rely on the subdivision of the fluid domain into small, stationary cells, such as tetrahedra, and solve the [Navier-Stokes equations](https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations) on each Eulerian (not moving) cell. In contrast, vortex methods rely on a Lagrangian (moving with the flow) description of the only the [vorticity](https://en.wikipedia.org/wiki/Vorticity)-containing region of the fluid domain and any solid boundaries present. This eliminates many of the difficulties present in traditional CFD. In addition, the form of the equations used also removes the pressure term from the Navier-Stokes equations, which is a large source of instability and extra effort in traditional CFD. This is why many new flow solvers for unsteady momentum-dominated flows (non-microscopic in scale) are implemented using vortex methods.

Omega3D aims to be an accurate combined Lagrangian-Eulerian fluid flow solver for unsteady flows with complex boundaries, with a greatly reduced reliance on meshing and parameter tuning for achieving accuracy and speed.


## Build the software
This code uses some C++17 features, so should compile on GCC 7, Clang 4, and MSVC 19.10 (Visual Studio 15 2017) compilers.

#### Prerequisites
Both the GUI and batch versions require CMake and Eigen to compile. To build the GUI version, users will also need GLFW3. These can be installed in Red Hat/Fedora with

    sudo dnf install cmake glfw3-devel eigen3-devel

or on Ubuntu with

    sudo apt-get install cmake glfw3-dev libeigen3-dev

or on OSX via [Homebrew](https://docs.brew.sh/Installation) with

    brew install cmake glfw eigen

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

The above commands should work verbatim on Linux and OSX. Don't ask about Windows - there's a calling convention issue preventing this from working.

#### Compile
Upon installation of the prerequisites, the following commands should build Omega3D.

    git clone git@github.com:Applied-Scientific-Research/Omega3D.git
    cd Omega3D
    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DUSE_OMP=ON -DUSE_VC=OFF ..
    make

If you were able to build and install Vc, then you should set `-DUSE_VC=ON` in the above `cmake` command.

On OSX, to get OpenMP parallelization of the solver, you may need to install GCC with brew (as above), and add a few more arguments to the `cmake` command:

    brew install gcc
    cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc-9 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-9 ..


## Run a simulation in the GUI
If you were able to build the software, you should be able to run

    ./Omega3D.bin

Upon running Omega2D, you will see a GUI floating over a black field. Using the *Select a simulation...* pull-down menu, you can quickly load and run a preset simulation. Let's load "flow over sphere".

At any time you can press *PAUSE* to pause the simulation or *Reset* to go back to the original conditions. At any time, you can left-click and drag on the flow field to move your point of view, or use the scroll wheel to zoom and unzoom. Space bar also pauses and unpauses the simulation. Note that some simulations quickly become large enough to take several seconds between updates. Don't worry: when you pause, the current simulation step will finish.

There are several collapsible headers which you can open to modify this simulation, those include *Simulation globals* such as viscosity and flow speed, *Flow structures* such as solid bodies, vortex blobs, and tracers, and *Rendering parameters*. Some changes you make in these fields will affect the simulation immediately, but most will require you to *Reset*.

### Run a batch job
If you already have an input file in JSON format, or you exported one from the GUI, you can run a batch (no GUI) simulation with

    ./Omega3Dbatch.bin input.json

Output will be written to the terminal and files to the working directory.


## To do
Tasks to consider or implement:

* Add GUI checkbox for saving the per-step status file
* Add measure feature for 3D grid of values
* Add option to only compute vels on fldpts when a vtu file is written - saves time
* Make the freestream a formula-entry system as well
* Standardize the core function selection: put it in the GUI, update the moments for aVRM, update flops, use new shader
* Get compute shader evaluation to work when called within vtu-writing
* Get BEM working for moving bodies
* Allow general rectangles - use the ips to scale panel sizes
* Use Eigen's Quaternion to represent rotations internally, but load them in as axis-angle
* Consider using [libigl](https://libigl.github.io/tutorial/#closest-points) to find point-mesh closest point queries
* Add some pics, maybe aGIF, to this readme
* Instead of manipulating the projection matrix, have the mouse change the view matrix (assume model matrix is unity), see [here](https://solarianprogrammer.com/2013/05/22/opengl-101-matrices-projection-view-model/) for a nice write-up on the three OpenGL matrices
* Add arcball rotation to the viewport - see [here](https://www.3dgep.com/understanding-the-view-matrix/) for some glm code
* Add other repos as submodules, or just by copying? `submodule add https://...xxx.git thirdparty/xxx`

## Thanks
This project is funded by the [National Institutes of Health (NIH)](https://www.nih.gov/) under grant number 1 R01 EB022180-01A1 ("A Fast High-Order CFD for Turbulent Flow Simulation in Cardio-Devices").

Thanks to [Omar Cornut](http://www.miracleworld.net/) for his [dear imgui](https://github.com/ocornut/imgui) library, file browser dialogs from [Imgui-IGS-Snippets](https://github.com/gileoo/Imgui-IGS-Snippets), png writing from [stb\_image\_write](https://github.com/nothings/stb/blob/master/stb_image_write.h), sol-prog's [OpenGL Tutorials](https://github.com/sol-prog/OpenGL-101), and Jim Susinno's [OpenGL-Boilerplate](https://github.com/jimbo00000/OpenGL-Boilerplate). Reading triangle mesh files was easy with [libigl](https://github.com/libigl/libigl/).

VRM code is functional thanks to jlblancoc for [Nanoflann](https://github.com/jlblancoc/nanoflann) (a header-only tree search library), and to all of the developers of [Eigen](http://eigen.tuxfamily.org/) (a C++ matrix/vector library). The BEM code also relies heavily on [Eigen](http://eigen.tuxfamily.org/). We also love [Vc](https://github.com/VcDevel/Vc), an excellent SIMD library by Matthias Kretz.

JSON reading and writing is thanks to [JSON for Modern C++](https://github.com/nlohmann/json) by [Niels Lohmann](http://nlohmann.me). XML output to VTK files is done using [tinyxml2](https://github.com/leethomason/tinyxml2) and [cppcodec](https://github.com/tplgy/cppcodec) for base64 encoding. And mathematical expression parsing came from [Lewis Van Winkle](https://codeplea.com/)'s [tinyexpr](https://github.com/codeplea/tinyexpr).

Many thanks to NBL for valuable discussions of architecture and C++ syntax and idioms.
