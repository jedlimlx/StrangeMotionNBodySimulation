//
// Created by qiuzi on 3/1/2022.
//

#ifndef SDESOLVER_TREENODE_H
#define SDESOLVER_TREENODE_H
#include "sde_solver_constants.h"
#include "vector"
using namespace std;

struct particle{
    long double x;
    long double y;
    long double vx;
    long double vy;
    long double ax;
    long double ay;
};

class TreeNode {
public:
    vector<TreeNode> children{}; // direct children of this cell
    TreeNode* parent{}; // direct parent of this cell
    int n_particles{}; // number of particles in this cell (recursive)
    vector<particle*> particles; // all particles in this cell (recursive)
    long double x_start{}; // bottom left
    long double y_start{}; // bottom left
    long double size{}; // size of cell
    long double x_com;
    long double y_com;

    TreeNode(long double x_start, long double y_start, long double size);

    void insert(particle* p); //recursive

    void calculate_force(particle* source, long double* force);
};


#endif //SDESOLVER_TREENODE_H
