cmake_minimum_required(VERSION 3.13)
project(Staurn)
add_executable(Staurn render.cpp)


set(CMAKE_CXX_STANDARD 14)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

target_link_libraries(Staurn ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )
