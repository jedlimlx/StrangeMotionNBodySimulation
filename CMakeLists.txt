cmake_minimum_required(VERSION 3.21)
project(SDESolver)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3 -pthread -lboost_math_tr1")

add_executable(SDESolver main.cpp TreeNode.cpp TreeNode.h sde_solver_constants.h)
