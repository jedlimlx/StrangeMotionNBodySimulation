//
// Created by qiuzi on 3/1/2022.
//


#include <vector>
#ifndef SDESOLVER_TREENODE_H
#define SDESOLVER_TREENODE_H
#include "sde_solver_constants.h"
#include <memory>


struct particle{
    long double x;
    long double y;
    long double vx;
    long double vy;
    long double ax;
    long double ay;
    particle(long double, long double);
};

class TreeNode {
public:
    std::vector<TreeNode> children; // direct children of this cell TODO use array
    int n_particles; // number of particles in this cell (recursive)
    std::vector<particle*> particles; // all particles in this cell (recursive)
    long double x_start; // bottom left
    long double y_start; // bottom left
    long double size; // size of cell
    long double x_com;
    long double y_com;
    long double** interp;

    TreeNode(long double x_start, long double y_start, long double size, long double** interp);

    void insert_particle(particle* p); //recursive

    void calculate_force(particle* source, long double* force);

    void clear();
};

long double eval_interp(long double **interp, long double point);

#endif //SDESOLVER_TREENODE_H
