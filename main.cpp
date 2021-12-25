#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <thread>
#include <chrono>

using namespace std;

//import coefficients for the force polynomial
int import_coeffs(const string& filename, int terms, long double coeffs[]){
    ifstream infile;
    infile.open(filename);
    //error checking
    if(! infile.is_open()){
        printf("Failed to open file");
        return 1;
    }

    string word;
    //load the data into the arrays
    for(int i = 0; i < terms; i++) {
        getline(infile, word);
        coeffs[i] = stod(word);
    }
    infile.close();
    return 0;
}

//polynomial for the force
long double force(long double r, const long double coeffs[], int terms){
    long double sum = 0;
    for (int i = 0; i < terms; i++){
        sum += coeffs[i] * pow(r, i);
    }
    return sum;
}

//reads in the initial positions
int import_initial_data(const string& filename, int length, double positions[]){
    ifstream infile;
    infile.open(filename);
    //error checking
    if(! infile.is_open()){
        cout << "Failed to open file" << endl;
        return 1;
    }
    //read in the data
    string word;
    for(int i = 0; i < length; i ++){
        getline(infile, word);
        positions[i] = stod(word);
    }
    infile.close();
    return 0;
}

//initialize positions of particles
int generate_initial_positions(const double initial_data[], long double* positions[], int particles, bool**mesh, int initial_data_length){
    default_random_engine generator(particles);
    uniform_real_distribution<long double> angle_distribution(0, 2 * M_PI);
    uniform_int_distribution<int> index_distribution(0,  initial_data_length - 1);

    long double r, theta, x, y;
    int index;
    for(int i = 0; i < particles; i++){
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
long double get_random(long double random_coeff, long double dt){
    static default_random_engine generator(time(NULL));
    normal_distribution<long double> distribution(0, sqrt(dt));
    return random_coeff * distribution(generator);
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

        force_val = force(r, coeffs, terms);
        rand_x = get_random(random_coeff, dt);
        rand_y = get_random(random_coeff, dt);
        v_x_f = v_x_i + (force_val * (x_original / r) * dt -
                         cd * v_x_i * dt +
                         0.1 *rand_x) / mass;
        v_y_f = v_y_i + (force_val * (y_original / r) * dt -
                         cd * v_y_i * dt +
                         0.1 * rand_y) / mass;
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
int solve_sde(long double* positions[], long double* velocities[], int N, int particles,long double t_start, long double t_end, long double coeffs[] , int terms, bool* mesh[],
              long double cd, long double random_coeff, long double mass, int n_threads){

    long double dt = (t_end - t_start)/N;

    auto start = chrono::high_resolution_clock::now();

    ofstream outfile;
    outfile.open("position.csv");

    float length = ((float) particles) / n_threads; //number of particles for each thread to process
    for (int t = 0; t < N; ++t) {
        thread threads[n_threads];
        for(int i = 0; i < n_threads; i++){
            threads[i] = thread(loop_for_particles, (i * length), ((i + 1) * length), positions, velocities, mesh, coeffs, terms, cd, random_coeff, mass, dt);
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
    const string FORCE_COEFFS_FILENAME = "forcecoeffs.csv"; //coefficients for force polynomial
    const int TERMS = 26; //number of terms in force polynomial
    const int INITIAL_DATA_LENGTH = 750357; //number of initial r values
    const string INITIAL_DATA_FILENAME = "initial_data.csv"; //initial r values
    const int PARTICLES = 20000; //number of particles to simulate
    const int MESH_FINENESS = 3000; //dimensions of mesh (MESH_FINENESS * MESH_FINENESS)
    const int N = 100000; //number of timesteps
    const int N_THREADS = 6;

    const long double VISCOSITY = 0.0010518; //dynamic viscosity of water
    const long double RADIUS = 75e-6; //radius of particle
    const long double DENSITY = 8960; //density of particles
    const long double MASS = (4.0/3) * DENSITY * M_PI * pow(RADIUS, 3);
    const long double CD = 6 * M_PI * VISCOSITY * RADIUS; //stokes drag
    const long double TEMPERATURE = 28 + 273.15; //temperature
    const long double KB = 1.38064852e-23; //boltzmann's constant
    const long double RANDOM_COEFFICIENT = sqrt(2 * KB * TEMPERATURE / CD); //coefficient in front of the dW term
    const int T_START = 0;
    const int T_END = 1000;


    auto *coeffs = (long double*) malloc(TERMS * sizeof(long double));
    //import polynomial coefficients
    if(import_coeffs(FORCE_COEFFS_FILENAME, TERMS, coeffs)){
        return 1;
    }
    auto *initial_data = (double*) malloc(INITIAL_DATA_LENGTH * sizeof(double));
    //import initial data
    if(import_initial_data(INITIAL_DATA_FILENAME, INITIAL_DATA_LENGTH, initial_data)){
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
    generate_initial_positions(initial_data, positions, PARTICLES, mesh, INITIAL_DATA_LENGTH);
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
    cout << force(0.0653622022084068693415, coeffs, TERMS) << endl;
    solve_sde(positions, velocities, N, PARTICLES, T_START, T_END, coeffs, TERMS, mesh, CD, RANDOM_COEFFICIENT, MASS, N_THREADS);

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