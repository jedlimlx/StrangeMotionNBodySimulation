CC=g++
CFLAGS=-I. -std=c++14 -O3 -pthread

default: SDESolver

SDESolver: SDESolver
	$(CC) -o SDESolver main.cpp sde_solver_constants.h TreeNode.h TreeNode.cpp $(CFLAGS)


