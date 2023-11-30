cmake_minimum_required(VERSION 3.11)

# Project name
project(Epidemic-Simulation-GUI)

set(CMAKE_CXX_STANDARD 20)

# Add an executable compiled from app.cpp
add_executable(Epidemic-Simulation-GUI app.cpp)

# The package PkgConfig detects needed GTK+ headers/library files
find_package(PkgConfig REQUIRED)
pkg_check_modules(GTKMM REQUIRED gtkmm-3.0)

# Setup CMake to use GTK+ to tell the compiler where to look for headers
# And the linker where to look for libraries
include_directories(${GTKMM_INCLUDE_DIRS})
link_directories(${GTKMM_LIBRARY_DIRS})

# Flags for the compiler
add_definitions(${GTKMM_CFLAGS_OTHER})

# Link the target to the GTK+ libraries
target_link_libraries(Epidemic-Simulation-GUI ${GTKMM_LIBRARIES})