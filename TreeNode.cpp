//
// Created by qiuzi on 3/1/2022.
//

#include <cstdlib>
#include <cmath>
#include <boost/math/tr1.hpp>
#include "TreeNode.h"
//using namespace std;
TreeNode::TreeNode(long double x_start, long double y_start, long double size) {
    children.reserve(4);
    this->x_start = x_start;
    this->y_start = y_start;
    this->size = size;
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
    if(size < 2 * RADIUS){
        for(int i = 0; i < n_particles; i++){
            vx_ave += particles[i]->vx;
            vy_ave += particles[i]->vy;
        }
        for(int i = 0; i < n_particles; i++){
            particles[i]->vx = vx_ave / n_particles;
            particles[i]->vy = vy_ave / n_particles;
        }
    }
    // terminal node
    if(n_particles == 1) {
        return;
        // generate children
    }else if(n_particles == 2){
        children.emplace_back(x_start, y_start, size / 2);
        children.emplace_back(x_start + size / 2, y_start, size / 2);
        children.emplace_back(x_start, y_start + size / 2, size / 2);
        children.emplace_back(x_start + size / 2, y_start + size / 2, size / 2);
        for(int i = 0; i < 4; i++){
            children[i].parent = this;
        }
        for (int i = 0 ; i < 2; i++){
            if(particles[i]->x < x_start + size / 2){
                if(particles[i]->y < y_start + size / 2){
                    children[0].insert(particles[i]);
                }else{
                    children[2].insert(particles[i]);
                }
            }else{
                if(particles[i]->y < y_start + size / 2){
                    children[1].insert(particles[i]);
                }else{
                    children[3].insert(particles[i]);
                }
            }
        }
    }if(p->x < x_start + size / 2){
        if(p->y < y_start + size / 2){
            children[0].insert(p);
        }else{
            children[2].insert(p);
        }
    }else{
        if(p->y < y_start + size / 2){
            children[1].insert(p);
        }else{
            children[3].insert(p);
        }
    }
}

void TreeNode::calculate_force(particle *source, long double *force) {
    long double r = sqrt(pow(x_com - source->x, 2) + pow(y_com - source->y, 2));
    if(size / r > MAX_ANGLE){
        for(int i = 0; i < 4; i++){
            children[i].calculate_force(source, force);
        }
    }
    long double f = COEFF * n_particles * QK * QK * boost::math::tr1::cyl_bessel_k(0.0, r * CAPL_LENGTH) / r;
    force[0] += f * (x_com - source->x);
    force[1] += f * (y_com - source->y);
}

