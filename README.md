# Intro.

SoftPOSIT is an algorithm propose by David etc. <sup>[1](#myfootnote1)</sup>
a model based camera algorithm for determining the pose (position and orientation) of a 3D object from a single 2D image when correspondences between object points and image points are not known. This algorithm integrates the Softassign technique for computing correspondences and the POSIT technique for computing object pose. The method finds the rotation and translation parameters of the camera with respect to an object.

# Existing Implementations

## matlab imp.

The author provides matlab implementation here: http://www.cfar.umd.edu/~daniel/Site_2/Code.html

## C (fortran) imp.

Someone else with a c/fortran implementation: https://github.com/hanshuo/softposit

## C++ imp.

This repo. provides my cpp implementation.


# How to Build

## Dependency
1. Eigen3

It mainly depends on the Eigen3 linear algebra library.


## Build Steps
The project is managed with CMake.

Suggested steps:

```
git clone https://github.com/wqhot/softposit-cpp
cd softposit
git checkout removearma
mkdir build
cd build
cmake ..
make
./softpositdemo
```

# TODO

1. could someone tell me how to do a serious (yet simple) logging ?  Boost Logging seems just too heavy for such a small lib.
2. bindings for Python? is this needed?

# Reference

<a name="myfootnote1">[1]</a>: David, Philip and Dementhon, Daniel and Duraiswami, Ramani and Samet, Hanan, "{SoftPOSIT: Simultaneous pose and correspondence determination}", Int. J. Comput. Vis., 59:259--284 (2004)
