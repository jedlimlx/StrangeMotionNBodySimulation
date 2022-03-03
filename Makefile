CC = gcc
opts = -O3 -Wall -g -lpthread -lboost_math_tr1 -I. -std=c++17

default: SDESolver

SDESolver: main.cpp TreeNode.cpp
	$(CC) -o sdesolver main.cpp TreeNode.cpp $(opts)

clean:
