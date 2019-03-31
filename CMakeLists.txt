cmake_minimum_required(VERSION 3.13)
project(SATURN)
add_executable(SATURN render.cpp)

set(CMAKE_CXX_STANDARD 14)
set(SATURN_INSTALL_DIR "${SATURN_SOURCE_DIR}/bin")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(X11 REQUIRED)

include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} ${X11_INCLUDE_DIR})

target_link_libraries(SATURN ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${X11_LIBRARY})

#cc -o "$1" "$1".c -I/opt/X11/include -L/opt/X11/lib -lX11