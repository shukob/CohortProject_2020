cmake_minimum_required(VERSION 3.16)

set(CMAKE_Fortran_COMPILER /usr/local/Cellar/gcc/10.1.0/bin/gfortran)
project(Task2Code CXX C Fortran)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_definitions(-DDUSHQUICK -DBLAS -DLINUX)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -framework Accelerate")

#find_library(gfortran
#        NAMES gfortran
#        HINTS "/usr/local/Cellar/gcc/10.1.0/lib/gcc/10")
#
#if (NOT gfortran)
#    message(STATUS "gfortran missing")
#endif()

add_executable(Task2Code
#        FC.cpp
        FC.h
        FC2.cpp
        matvec.cpp
        matvec.h
        prod.f90)

link_directories("/usr/local/Cellar/gcc/10.1.0/lib/gcc/10")
target_link_libraries(Task2Code gfortran stdc++)
