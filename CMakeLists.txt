cmake_minimum_required(VERSION 3.9.2)
project(SDESolver)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")

find_package(Boost 1.78 COMPONENTS math_tr1 REQUIRED)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    add_executable(${PROJECT_NAME} main.cpp sde_solver_constants.h TreeNode.cpp TreeNode.h)
    MESSAGE("${Boost_INCLUDE}")
    target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
endif()
