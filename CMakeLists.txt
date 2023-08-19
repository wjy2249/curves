cmake_minimum_required(VERSION 3.20)
project(curves)
set(SOURCES  "test.cpp" "NURBS3D/NURBS3D.hpp")
add_executable(curves ${SOURCES})