cmake_minimum_required(VERSION 3.21)
project(SDESolver)

set(CMAKE_CXX_STANDARD 14)

add_executable(SDESolver main.cpp sde_solver_constants.h TreeNode.cpp TreeNode.h)
