//
// Created by qiuzi on 3/1/2022.
//
#include "TreeNode.h"
#include <iostream>
long double eval_interp(long double **interp, long double point);

//using namespace std;
TreeNode::TreeNode(long double x_start, long double y_start, long double size, long double** interp) {
    children = new std::vector<TreeNode>;
    children->reserve(4);
    this->x_start = x_start;
    this->y_start = y_start;
    this->size = size;
    this->interp = interp;
    x_com = 0;
    y_com = 0;
    n_particles = 0;
}

void TreeNode::insert(particle* p) {
    particles.push_back(p);
    x_com = (x_com * n_particles + p->x) / (n_particles + 1);
    y_com = (y_com * n_particles + p->y) / (n_particles + 1);
    n_particles++;
    //resolve collisions
    long double vx_ave = 0;
    long double vy_ave = 0;
    if (size < 2 * SDESOLVER_RADIUS) {
        for (int i = 0; i < n_particles; i++) {
            vx_ave += particles[i]->vx;
            vy_ave += particles[i]->vy;
        }
        for (int i = 0; i < n_particles; i++) {
            particles[i]->vx = vx_ave / n_particles;
            particles[i]->vy = vy_ave / n_particles;
        }
    }
    // terminal node
    if (n_particles == 1) {
        return;
        // generate children
    } else if (n_particles == 2) {
        children->emplace_back(x_start, y_start, size / 2, interp);
        children->emplace_back(x_start + size / 2, y_start, size / 2, interp);
        children->emplace_back(x_start, y_start + size / 2, size / 2, interp);
        children->emplace_back(x_start + size / 2, y_start + size / 2, size / 2, interp);
        for (int i = 0; i < 2; i++) {
            if (particles[i]->x < x_start + size / 2) {
                if (particles[i]->y < y_start + size / 2) {
                    children.insert(particles[i]);
                } else {
                    children[2].insert(particles[i]);
                }
            } else {
                if (particles[i]->y < y_start + size / 2) {
                    children[1].insert(particles[i]);
                } else {
                    children[3].insert(particles[i]);
                }
            }
        }
    } else {
        if (p->x < x_start + size / 2) {
            if (p->y < y_start + size / 2) {
                children[0].insert(p);
            } else {
                children[2].insert(p);
            }
        }else{
            if (p->y < y_start + size / 2) {
                children[1].insert(p);
            } else {
                children[3].insert(p);
            }
        }
    }
}

void TreeNode::calculate_force(particle *source, long double* force) {
    long double r = sqrt(pow(x_com - source->x, 2) + pow(y_com - source->y, 2));
    if(size / r > SDESOLVER_MAX_ANGLE){
        for(int i = 0; i < 4; i++){
            children[i].calculate_force(source, force);
        }
    }
    long double f = eval_interp(interp, r) / r;
    force[0] += f * (x_com - source->x);
    force[1] += f * (y_com - source->y);
}

void TreeNode::clear() {
    children.clear();
    particles.clear();
    n_particles = 0;
}

TreeNode::~TreeNode() {
    delete children;
}

#define MAX_VALUE (-0.000123193)
long double eval_interp(long double **interp, long double point){
    if(point < 0.001) return MAX_VALUE;
    long double *c0 = interp[0], *c1 = interp[1], *bd = interp[2];
    unsigned i=0;
    while(c0[i]>0){
        if(point<bd[i+1])break;
        ++i;
    } // TODO: use some binary search
    if(c0[i]==0)--i;
    return c0[i]+c1[i]*point;
}

particle::particle(long double x, long double y) {
    this->x = x;
    this->y = y;
    vx = 0;
    vy = 0;
    ax = 0;
    ay = 0;
}
