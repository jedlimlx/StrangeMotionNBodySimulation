#include "TreeNode.h"
#include <iostream>
#include <math.h>
long double eval_interp(long double **interp, long double point);

//using namespace std;
TreeNode::TreeNode(long double x_start, long double y_start, long double size, long double** interp, int layer, int printLayers) {
    children.reserve(4);

    this -> x_start = x_start;
    this -> y_start = y_start;
    this -> size = size;
    this -> interp = interp;
    this -> layer = layer;

    x_com = 0;
    y_com = 0;
    n_particles = 0;
    stop = 1;

    this -> printLayers = printLayers;
}

void TreeNode::insert_particle(particle* p) {
    particles.push_back(p);
    if (printLayers)
        std::cout << layer << std::endl;
    if (isnan((x_com * n_particles + p -> x) / (n_particles + 1))){
        std::cout << "x_com is nan" << x_com << ", " << n_particles << ", " << p->x << std::endl;
        std::exit(1);
    }

    x_com = (x_com * n_particles + p -> x) / (n_particles + 1);
    y_com = (y_com * n_particles + p -> y) / (n_particles + 1);
    n_particles++;

    // terminal node
    if (n_particles == 1 || size < (2 * SDESOLVER_RADIUS)) {
        return;
        // generate children
    } else if (n_particles == 2) {
        stop = 0;
        children.emplace_back(x_start, y_start, size / 2, interp, layer + 1, printLayers);
        children.emplace_back(x_start + size / 2, y_start, size / 2, interp, layer + 1, printLayers);
        children.emplace_back(x_start, y_start + size / 2, size / 2, interp, layer + 1, printLayers);
        children.emplace_back(x_start + size / 2, y_start + size / 2, size / 2, interp, layer + 1, printLayers);
        for (int i = 0; i < 2; i++) {
            if (particles[i] -> x < x_start + size / 2) {
                if (particles[i] -> y < y_start + size / 2) {
                    children[0].insert_particle(particles[i]);
                } else {
                    children[2].insert_particle(particles[i]);
                }
            } else {
                if (particles[i] -> y < y_start + size / 2) {
                    children[1].insert_particle(particles[i]);
                } else {
                    children[3].insert_particle(particles[i]);
                }
            }
        }
    } else {
//        std::cout << "multiple" << std::endl;
        if (p -> x < x_start + size / 2) {
            if (p -> y < y_start + size / 2) {
                children[0].insert_particle(p);
            } else {
                children[2].insert_particle(p);
            }
        } else {
            if (p -> y < y_start + size / 2) {
                children[1].insert_particle(p);
            } else {
                children[3].insert_particle(p);
            }
        }
    }
}

void TreeNode::calculate_force(particle* source, long double* force) {
    if (n_particles == 0) return;

    long double r = 1e-12 + sqrt(pow(x_com - source -> x, 2) + pow(y_com - source -> y, 2));
    if (size / r > SDESOLVER_MAX_ANGLE && !stop) {
        for (int i = 0; i < 4; i++) {
            children[i].calculate_force(source, force);
        }
    } else {
        long double f;
        if (r < 2 * SDESOLVER_RADIUS) {
            f = 0;
        } else {
            f = eval_interp(interp, r) / r;
            if (isnan(f)) {
                std::cout << x_com << std::endl;
            }
        }

        // std::cout << f * r << " " << r << std::endl;
        // std::cout << "x: " << x_com << " " << source -> x << std::endl;
        // std::cout << "y: " << y_com << " " << source -> y << std::endl;
        force[0] += n_particles * f * (x_com - source->x);
        force[1] += n_particles * f * (y_com - source->y);
    }
}

void TreeNode::clear() {
    children.clear();
    particles.clear();
    n_particles = 0;
}

void TreeNode::resolve_collisions() {
    // resolve collisions
    if (size > 3.5 * SDESOLVER_RADIUS && !stop) {
        for (int i = 0; i < 4; i++) {
            children[i].resolve_collisions();
        }
    } else {
        long double vx_ave = 0;
        long double vy_ave = 0;
        if (size < (3 * SDESOLVER_RADIUS)) {
            for (int i = 0; i < n_particles; i++) {
                vx_ave += particles[i] -> vx;
                vy_ave += particles[i] -> vy;
            }

            for (int i = 0; i < n_particles; i++) {
                particles[i] -> vx = vx_ave / n_particles;
                particles[i] -> vy = vy_ave / n_particles;
            }
        }
    }
}

void TreeNode::get_particles() {
    if (n_particles == 0) return;
    if (!stop) {
        if (children.empty()) return;
        for (int i = 0; i < 4; i++) {
            children[i].get_particles();
        }
    } else {
        std::cout << x_com << "," << y_com << std::endl;
        //std::cout << particles[0]->x << " " << particles[0]->y << std::endl;
    }
}

#define MAX_VALUE (-0.000123193)
long double eval_interp(long double **interp, long double point) {
    if (point < 0.001) return MAX_VALUE;

    unsigned i = 0;
    long double *c0 = interp[1], *c1 = interp[0], *bd = interp[2];
    while (c0[i] > 0) {
        if (point < bd[i + 1]) break;
        ++i;
    } // TODO: use some binary search

    if (c0[i] == 0) --i;
    return c0[i] + c1[i] * point;
}

particle::particle(long double x, long double y) {
    this -> x = x;
    this -> y = y;

    vx = 0;
    vy = 0;
    ax = 0;
    ay = 0;
}
