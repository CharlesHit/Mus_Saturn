README.md

# Lighting, Shading and Hidder

## User Interface

THIS PROGRAM WILL COST A LOT OF MEMORY, CAREFULLY DECIDED WHETER TO OPEN IT!!!

Delete main.cpp - main( ) - glutIdleFunc(idle); to stop rotation.

Active main.cpp - main( ) - scalarization/rotation/translation to change camera

Go to data.h to change windows' size.


## Compiler

Using cmake -G "$Xcode" to compiler in xcode.

For vs2017, reading this:https://devblogs.microsoft.com/cppblog/cmake-support-in-visual-studio/#build-cmake

For CLion, since CMake is its default way to construct a project, things are the same.

## Program structure

data.h   -|
matrix.h -|
          |
      camera.h
          |
        draw.h
          |
      render.cpp

data.h is the basic config of this program. Including position of eyes, gaze point, and direction of camera.

matrix.h contains common matrix operations.

camera.h includes the construction of camera, and projection functions.

draw.h includes line, pixel and polygon filling algorithm, plus Torus, Sphere, Cone's parametric rendering.

render.cpp is the drawing process.

