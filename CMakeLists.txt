cmake_minimum_required(VERSION 3.15)
project(verlet)

set(CMAKE_CXX_STANDARD 14)

add_executable(verlet
        main.cpp Particle.cpp Particle.h Utils.cpp Utils.h)
