cmake_minimum_required(VERSION 3.9.2)
project(SDESolver)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")

add_executable(${PROJECT_NAME} main.cpp sde_solver_constants.h TreeNode.cpp TreeNode.h)
