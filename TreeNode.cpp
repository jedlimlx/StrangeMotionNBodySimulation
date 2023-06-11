#include "TreeNode.h"
#include <iostream>
#include <math.h>
#include <set>
#include <queue>
#include <unordered_set>

long double eval_interp(long double **interp, long double point);
void resolve_collisions_array(std::vector<particle*> particles);

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
    if (!stop) {
        for (int i = 0; i < 4; i++) {
            children[i].clear();
        }
    }

    children.clear();
    particles.clear();
    n_particles = 0;
}

void TreeNode::resolve_collisions() {
    if (n_particles <= 1) return;

    // resolve collisions
    if (size > 64 * SDESOLVER_RADIUS && !stop) {
        for (int i = 0; i < 4; i++) {
            children[i].resolve_collisions();
        }
    } else {
        if (size < 64 * SDESOLVER_RADIUS) {
            // std::cout << "resolving collisions..." << std::endl;

            /*
            for (int i = 0; i < n_particles; i++) {
                vx_ave += particles[i] -> vx;
                vy_ave += particles[i] -> vy;
            }

            for (int i = 0; i < n_particles; i++) {
                particles[i] -> vx = vx_ave / n_particles;
                particles[i] -> vy = vy_ave / n_particles;
            }
            */

            resolve_collisions_array(particles);
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

#define MAX_VALUE (0)
long double eval_interp(long double **interp, long double point) {
    if (point < 0.0001) return MAX_VALUE;

    unsigned left = 0;
    unsigned right = 1999; // todo insert correct number
    unsigned middle = (left + right) / 2;
    long double *c0 = interp[1], *c1 = interp[0], *bd = interp[2];

    while (c0[middle] > 0 && left < right) {
        // if (point < bd[i + 1]) break;
        // ++i;
        middle = (left + right) / 2;
        if (point > bd[middle]) {
            left = middle + 1;
        } else if (point < bd[middle]) {
            right = middle - 1;
        }

        // std::cout << middle << " " << bd[middle] << " " << point << std::endl;
    }

    // if (c0[middle] == 0) --middle;
    return c0[middle] + c1[middle] * point;
}

void resolve_collisions_array(std::vector<particle*> particles) {
    int i, j, k;
    long double max_v = 0, v;
    for (i = 0; i < particles.size(); i++) {
        v = pow(particles[i]->vx, 2) + pow(particles[i]->vy, 2);
        if (max_v < v) max_v = v;
    }

    std::vector<std::vector<int>> adj_list;

    int num_sub_timesteps = SDESOLVER_COLLISION_TOLERANCE * sqrt(max_v) * SDESOLVER_DT / SDESOLVER_RADIUS;
    if (num_sub_timesteps > 1) {
        for (i = 0; i < particles.size(); i++) {
            particles[i]->x -= particles[i]->vx * SDESOLVER_DT;
            particles[i]->y -= particles[i]->vy * SDESOLVER_DT;

            std::vector<int> adj;
            adj_list.push_back(adj);
        }

        bool collided = false;
        for (i = 0; i < num_sub_timesteps; i++) {
            // building adjacency list
            for (j = 0; j < particles.size(); j++) {
                for (k = 0; k < particles.size(); k++) {
                    if (j >= k) continue;
                    if (pow(particles[j]->x - particles[k]->x, 2) + pow(particles[j]->y - particles[k]->y,
                                                                        2) <= pow(SDESOLVER_RADIUS + 1e-6, 2)) {
                        adj_list[j].push_back(k);
                        adj_list[k].push_back(j);

                        collided = true;
                    }

                    /*
                    if (pow(particles[j]->x - particles[k]->x, 2) + pow(particles[j]->y - particles[k]->y,
                                                                        2) <= pow(SDESOLVER_RADIUS + 1e-7, 2)) {
                        adj_list[j].push_back(k);
                        adj_list[k].push_back(j);

                        long double vx_diff = particles[j]->vx - particles[k]->vx;
                        long double vy_diff = particles[j]->vy - particles[k]->vy;
                        long double x_diff = particles[j]->x - particles[k]->x;
                        long double y_diff = particles[j]->y - particles[k]->y;

                        if (vx_diff == 0 || vy_diff == 0) continue;

                        collided = true;

                        long double time1 = (-x_diff*vx_diff-y_diff*vy_diff
                                            -pow(pow(SDESOLVER_RADIUS,2)*(pow(vx_diff,2)+pow(vy_diff,2))
                                            -pow(vy_diff*x_diff-vx_diff*y_diff,2),0.5))/(pow(vx_diff,2)+pow(vy_diff,2));
                        long double time2 = (-x_diff*vx_diff-y_diff*vy_diff
                                            -pow(pow(SDESOLVER_RADIUS,2)*(pow(vx_diff,2)+pow(vy_diff,2))
                                            +pow(vy_diff*x_diff-vx_diff*y_diff,2),0.5))/(pow(vx_diff,2)+pow(vy_diff,2));

                        long double time;
                        if (time1 > time2 && time1 < 0) time = time1;
                        else if (time2 > time1 && time2 < 0) time = time2;
                        else time = 0;

                        if (time > SDESOLVER_DT) continue;

                        particles[j]->x += particles[j]->vx * time;
                        particles[j]->y += particles[j]->vy * time;

                        particles[k]->x += particles[k]->vx * time;
                        particles[k]->y += particles[k]->vy * time;
                    }
                     */
                }
            }

            if (collided) {
                // store velocity sums
                long double vx_sum = 0;
                long double vy_sum = 0;

                // perform BFS
                std::queue<int> queue;
                std::set<int> visited;
                std::vector<particle*> lst;

                int updated[particles.size()] = { 0 };
                for (j = 0; j < particles.size(); j++) {
                    if (updated[j] == 1) continue;

                    if (adj_list[j].size() == 0) continue;

                    vx_sum = 0;
                    vy_sum = 0;

                    queue.push(j);
                    while (!queue.empty()) {
                        int curr = queue.front();
                        visited.insert(curr);
                        queue.pop();

                        for (k = 0; k < adj_list[curr].size(); k++) {
                            if (visited.count(adj_list[curr][k]) == 0) {
                                queue.push(adj_list[curr][k]);
                            }
                        }

                        updated[curr] = 1;
                        lst.push_back(particles[curr]);

                        vx_sum += particles[curr]->vx;
                        vy_sum += particles[curr]->vy;
                    }

                    for (j = 0; j < lst.size(); j++) {
                        lst[j]->vx = vx_sum / lst.size();
                        lst[j]->vy = vy_sum / lst.size();
                    }

                    lst.clear();
                    visited.clear();
                }

                for (j = 0; j < adj_list.size(); j++) {
                    adj_list[j].clear();
                }
            }

            // perform timestep
            for (j = 0; j < particles.size(); j++) {
                particles[j]->x += particles[j]->vx * SDESOLVER_DT / num_sub_timesteps;
                particles[j]->y += particles[j]->vy * SDESOLVER_DT / num_sub_timesteps;
            }
        }

        // collision detection one last time
        collided = false;

        // building adjacency list
        for (j = 0; j < particles.size(); j++) {
            for (k = 0; k < particles.size(); k++) {
                if (j >= k) continue;
                if (pow(particles[j]->x - particles[k]->x, 2) + pow(particles[j]->y - particles[k]->y,
                                                                    2) <= pow(SDESOLVER_RADIUS, 2)) {
                    adj_list[j].push_back(k);
                    adj_list[k].push_back(j);

                    collided = true;
                }

                /*
                if (pow(particles[j]->x - particles[k]->x, 2) + pow(particles[j]->y - particles[k]->y,
                                                                    2) <= pow(SDESOLVER_RADIUS + 1e-7, 2)) {
                    adj_list[j].push_back(k);
                    adj_list[k].push_back(j);

                    long double vx_diff = particles[j]->vx - particles[k]->vx;
                    long double vy_diff = particles[j]->vy - particles[k]->vy;
                    long double x_diff = particles[j]->x - particles[k]->x;
                    long double y_diff = particles[j]->y - particles[k]->y;

                    if (vx_diff == 0 || vy_diff == 0) continue;

                    collided = true;

                    long double time1 = (-x_diff*vx_diff-y_diff*vy_diff
                                        -pow(pow(SDESOLVER_RADIUS,2)*(pow(vx_diff,2)+pow(vy_diff,2))
                                        -pow(vy_diff*x_diff-vx_diff*y_diff,2),0.5))/(pow(vx_diff,2)+pow(vy_diff,2));
                    long double time2 = (-x_diff*vx_diff-y_diff*vy_diff
                                        -pow(pow(SDESOLVER_RADIUS,2)*(pow(vx_diff,2)+pow(vy_diff,2))
                                        +pow(vy_diff*x_diff-vx_diff*y_diff,2),0.5))/(pow(vx_diff,2)+pow(vy_diff,2));

                    long double time;
                    if (time1 > time2 && time1 < 0) time = time1;
                    else if (time2 > time1 && time2 < 0) time = time2;
                    else time = 0;

                    if (time > SDESOLVER_DT) continue;

                    particles[j]->x += particles[j]->vx * time;
                    particles[j]->y += particles[j]->vy * time;

                    particles[k]->x += particles[k]->vx * time;
                    particles[k]->y += particles[k]->vy * time;
                }
                 */
            }
        }

        if (collided) {
            // store velocity sums
            long double vx_sum = 0;
            long double vy_sum = 0;

            // perform BFS
            std::queue<int> queue;
            std::set<int> visited;
            std::vector<particle*> lst;

            int updated[particles.size()] = { 0 };
            for (j = 0; j < particles.size(); j++) {
                if (updated[j] == 1) continue;

                if (adj_list[j].size() == 0) continue;

                vx_sum = 0;
                vy_sum = 0;

                queue.push(j);
                while (!queue.empty()) {
                    int curr = queue.front();
                    visited.insert(curr);
                    queue.pop();

                    for (k = 0; k < adj_list[curr].size(); k++) {
                        if (visited.count(adj_list[curr][k]) == 0) {
                            queue.push(adj_list[curr][k]);
                        }
                    }

                    updated[curr] = 1;
                    lst.push_back(particles[curr]);

                    vx_sum += particles[curr]->vx;
                    vy_sum += particles[curr]->vy;
                }

                for (j = 0; j < lst.size(); j++) {
                    lst[j]->vx = vx_sum / lst.size();
                    lst[j]->vy = vy_sum / lst.size();
                }


                lst.clear();
                visited.clear();
            }

            for (j = 0; j < adj_list.size(); j++) {
                adj_list[j].clear();
            }
        }

        // perform timestep
        for (j = 0; j < particles.size(); j++) {
            particles[j]->x += particles[j]->vx * SDESOLVER_DT / num_sub_timesteps;
            particles[j]->y += particles[j]->vy * SDESOLVER_DT / num_sub_timesteps;
        }
    }
}

particle::particle(long double x, long double y) {
    this -> x = x;
    this -> y = y;

    vx = 0;
    vy = 0;
    ax = 0;
    ay = 0;
}
