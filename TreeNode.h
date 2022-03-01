//
// Created by qiuzi on 3/1/2022.
//

#ifndef SDESOLVER_TREENODE_H
#define SDESOLVER_TREENODE_H
#include "sde_solver_constants.h"

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
    TreeNode* children; // direct children of this cell
    TreeNode* parent; // direct parent of this cell
    int n_particles; // number of particles in this cell (recursive)
    particle* particles; // all particles in this cell (recursive)
    long double x_start;
    long double y_start;

    TreeNode(){

    }
};


#endif //SDESOLVER_TREENODE_H
