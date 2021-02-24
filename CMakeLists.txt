cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_CXX_STANDARD 14)
project(gauss2d)
# Set source directory
set(SOURCE_DIR "src")
set(HEADER_DIR "include")
include_directories(${HEADER_DIR})
file (GLOB SOURCE_FILES "src/*.cc")
file (GLOB HEADER_FILES "include/*.h")
# Generate Python module
find_package(pybind11 REQUIRED)
pybind11_add_module(_gauss2d ${HEADER_FILES} ${SOURCE_FILES} "gauss2d/pybind11.cc")
