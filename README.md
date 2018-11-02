# omega3d
High-performance 3D vortex methods solver with easy GUI

## Build and run

    git clone git@github.com:Applied-Scientific-Research/Omega3D.git
    cd Omega3D/src
    make
    ./Omega3D.bin

## To do

* Rework into struct of arrays to allow Vc to work properly
* This requires more work for compute and draw shaders (one buffer per array)
* But should be nice and flexible - and easier than all the `4*idx+3` crap


