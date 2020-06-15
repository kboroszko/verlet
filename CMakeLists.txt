cmake_minimum_required(VERSION 3.15)
project(verlet)



#set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD 11)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(verlet
        main_seq.cpp Particle.cpp Particle.h Utils.cpp Utils.h)

add_executable(verlet_par
        main_par.cpp Particle.cpp Particle.h Utils.cpp Utils.h)


include_directories(SYSTEM ${MPI_INCLUDE_PATH})

target_link_libraries(verlet_par ${MPI_CXX_LIBRARIES})
set_target_properties(verlet_par PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} -g -O3 ")