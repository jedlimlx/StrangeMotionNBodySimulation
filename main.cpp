#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <thread>
#include <chrono>
#include "sde_solver_constants.h"
#include "boost/math/tr1.hpp"


using namespace std;

//import coefficients for the force polynomial
int import_coeffs(long double coeffs[]){
    ifstream infile;
    infile.open(FORCE_COEFFS_FILENAME);
    //error checking
    if(! infile.is_open()){
        printf("Failed to open file");
        return 1;
    }

    string word;
    //load the data into the arrays
    for(int i = 0; i < TERMS; i++) {
        getline(infile, word);
        coeffs[i] = stod(word);
    }
    infile.close();
    return 0;
}

//polynomial for the force
long double force(long double r, const long double coeffs[]){
    long double sum = 0;
    for (int i = 0; i < TERMS; i++){
        sum += coeffs[i] * pow(r, i);
    }
    return sum;
}

//reads in the initial positions
int import_initial_data(double positions[]){
    ifstream infile;
    infile.open(INITIAL_DATA_FILENAME);
    //error checking
    if(! infile.is_open()){
        cout << "Failed to open file" << endl;
        return 1;
    }
    //read in the data
    string word;
    for(int i = 0; i < INITIAL_DATA_LENGTH; i ++){
        getline(infile, word);
        positions[i] = stod(word);
    }
    infile.close();
    return 0;
}



//initialize positions of particles
int generate_initial_positions(const double initial_data[], long double* positions[], bool**mesh){
    default_random_engine generator(PARTICLES);
    uniform_real_distribution<long double> angle_distribution(0, 2 * M_PI);
    uniform_int_distribution<int> index_distribution(0,  INITIAL_DATA_LENGTH - 1);

    long double r, theta, x, y;
    int index;
    for(int i = 0; i < PARTICLES; i++){
        index = index_distribution(generator);
        r = initial_data[index];
        theta = angle_distribution(generator);
        x = r * cos(theta);
        y = r * sin(theta);
        while(mesh[1500 + ((int) (x / 75E-6))][1500 + ((int) (y / 75E-6))]) {
            theta = angle_distribution(generator);
            x = r * cos(theta);
            y = r * sin(theta);
        }
        positions[i][0] = x;
        positions[i][1] = y;
    }
    return 0;
}

//rng time
long double get_random(){
    static default_random_engine generator(time(NULL));
    normal_distribution<long double> distribution(0, sqrt(DT));
    return RANDOM_COEFFICIENT * distribution(generator);
}

//split up the loop for sde solving
void loop_for_particles(int start, int end, long double* positions[], long double* velocities[], bool* mesh[], long double coeffs[], int terms,
                        long double cd, long double random_coeff, long double mass, long double dt){
    long double x_original, y_original, x, y, v_x_i, v_y_i, r, v_x_f, v_y_f, force_val, rand_x, rand_y;
//    if(start == 0) {
//        cout << get_random(random_coeff, dt) << endl;
//    }
    for (int i = start; i < end; ++i) {
        x_original = positions[i][0];
        y_original = positions[i][1];
        v_x_i = velocities[i][0];
        v_y_i = velocities[i][1];
        mesh[1500 + (int) (x_original / 75E-6)][1500 + (int) (y_original / 75E-6)] = false;
        //collision detection
        r = sqrt(pow(x_original, 2) + pow(y_original, 2));

        force_val = force(r, coeffs);
        rand_x = get_random();
        rand_y = get_random();
        v_x_f = v_x_i + (force_val * (x_original / r) * dt -
                         cd * v_x_i * dt +
                         rand_x) / mass;
        v_y_f = v_y_i + (force_val * (y_original / r) * dt -
                         cd * v_y_i * dt +
                         rand_y) / mass;
//        if (i == 0) {
//            cout << y << ", " << (force(r, coeffs, terms) * (y / r) * dt -
//                                  cd * v_y_i * dt +
//                                  0.001 * get_random(random_coeff, dt)) / mass << ", " << v_y_i << endl;
//        }
        velocities[i][0] = v_x_f;
        velocities[i][1] = v_y_f;
        x = x_original + v_x_f * dt;
        y = y_original + v_y_f * dt;
        if (x > 0.09 || x < -0.09 || y > 0.09 || y < -0.09 ||
            mesh[1500 + (int) (x / 75E-6)][1500 + (int) (y / 75E-6)]) {
            velocities[i][0] = 0;
            velocities[i][1] = 0;
            mesh[1500 + (int) (x_original / 75E-6)][1500 + (int) (y_original / 75E-6)] = true;
        } else {
            positions[i][0] = x;
            positions[i][1] = y;
            mesh[1500 + (int) (x / 75E-6)][1500 + (int) (y / 75E-6)] = true;
        }
    }
}

//euler's method
int solve_sde(long double* positions[], long double* velocities[], int particles, long double coeffs[] , bool* mesh[],
long double cd, long double random_coeff, long double mass, int n_threads){
    long double dt = (T_END - T_START)/N;

    auto start = chrono::high_resolution_clock::now();

    ofstream outfile;
    outfile.open("position.csv");


    float length = ((float) particles) / n_threads; //number of particles for each thread to process
    for (int t = 0; t < N; ++t) {
        thread threads[n_threads];
        for(int i = 0; i < n_threads; i++){
            threads[i] = thread(loop_for_particles, (i * length), ((i + 1) * length), positions, velocities, mesh, coeffs, TERMS, cd, random_coeff, mass, dt);
        }
        for(auto& thread: threads){
            thread.join();
        }
        outfile << sqrt(pow(positions[0][0], 2) + pow(positions[0][1], 2))<< endl;

        if(t % 1000 == 0){
            cout << chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - start).count() << endl;
        }
    }
    outfile.close();
    return 0;
}


int main() {
    cout << MASS << endl;
    return 0;
    auto *coeffs = (long double*) malloc(TERMS * sizeof(long double));
    //import polynomial coefficients
    if(import_coeffs(coeffs)){
        return 1;
    }
    auto *initial_data = (double*) malloc(INITIAL_DATA_LENGTH * sizeof(double));
    //import initial data
    if(import_initial_data(initial_data)){
        return 1;
    }

    //initialize 2d array for positions
    auto** positions = (long double**) malloc(PARTICLES * sizeof(long double*));
    for(int i = 0; i < PARTICLES; i++){
        positions[i] = (long double*) malloc(2 * sizeof(long double));
    }
    //initialize 2d array for collision detection
    auto** mesh = (bool**) malloc(MESH_FINENESS * sizeof(bool*));
    for(int i = 0; i < MESH_FINENESS; i++){
        mesh[i] = (bool*) malloc(MESH_FINENESS * sizeof(bool));
    }
    //initialize to false
    for (int i = 0; i < MESH_FINENESS; i++) {
        for (int j = 0; j < MESH_FINENESS; j++) {
            mesh[i][j] = false;
        }
    }
    //initialize positions
    generate_initial_positions(initial_data, positions, mesh);
    free(initial_data); // free up memory since its not needed any more

    //allocate memory for initial velocities array
    auto** velocities = (long double**) malloc(PARTICLES * sizeof(long double*));
    for(int i = 0; i < PARTICLES; i++){
        velocities[i] = (long double*) malloc(2 * sizeof(long double));
    }
    //set initial velocities to 0
    for (int i = 0; i < PARTICLES; i++) {
        velocities[i][0] = 0;
        velocities[i][1] = 0;
    }

    solve_sde(positions, velocities, PARTICLES, coeffs, mesh, CD, RANDOM_COEFFICIENT, MASS, N_THREADS);

    ofstream outfile;
    outfile.open("final_positions.csv");
    for (int i = 0; i < PARTICLES; ++i) {
        outfile << positions[i][0] << "," << positions[i][1] << endl;
    }

    outfile.close();

    free(velocities);
    free(coeffs);
    free(positions);
    free(mesh);
    return 0;
}

//int read_csv(const string& filename, int length, double r[], double f[]){//read in the csv for force interpolation
//    ifstream infile;
//    infile.open(filename);
//    //error checking
//    if(! infile.is_open()){
//        printf("Failed to open file");
//        return 1;
//    }
//
//    string line, word;
//    for(int i = 0; i < length; i++){
//        //read one line
//        getline(infile, line);
//        stringstream s(line);
//        //load the data into the arrays
//        getline(s, word, ',');
//        r[i] = stof(word);
//        getline(s, word, ',');
//        f[i] = stof(word);
//    }
//    return 0;
//}