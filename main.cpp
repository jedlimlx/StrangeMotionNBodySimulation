#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <thread>
#include <chrono>
#include "sde_solver_constants.h"
#include "TreeNode.h"


using namespace std;

//import coefficients for the force polynomial
int import_coeffs(long double coeffs[]){
    ifstream infile;
    infile.open(SDESOLVER_FORCE_COEFFS_FILENAME);
    //error checking
    if(! infile.is_open()){
        printf("Failed to open file");
        return 1;
    }

    string word;
    //load the data into the arrays
    for(int i = 0; i < SDESOLVER_TERMS; i++) {
        getline(infile, word);
        coeffs[i] = stod(word);
    }
    infile.close();
    return 0;
}

//polynomial for the force
long double force(long double r, const long double coeffs[]){
    long double sum = 0;
    for (int i = 0; i < SDESOLVER_TERMS; i++){
        sum += coeffs[i] * pow(r, i);
    }
    return sum;
}

//reads in the initial positions
int import_initial_data(double positions[]){
    ifstream infile;
    infile.open(SDESOLVER_INITIAL_DATA_FILENAME);
    //error checking
    if(! infile.is_open()){
        cout << "Failed to open file" << endl;
        return 1;
    }
    //read in the data
    string word;
    for(int i = 0; i < SDESOLVER_INITIAL_DATA_LENGTH; i ++){
        getline(infile, word);
        positions[i] = stod(word);
    }
    infile.close();
    return 0;
}



//initialize positions of particles
int generate_initial_positions(const double initial_data[], long double* positions[], bool**mesh){
    default_random_engine generator(SDESOLVER_PARTICLES);
    uniform_real_distribution<long double> angle_distribution(0, 2 * M_PI);
    uniform_int_distribution<int> index_distribution(0,  SDESOLVER_INITIAL_DATA_LENGTH - 1);

    long double r, theta, x, y;
    int index;
    for(int i = 0; i < SDESOLVER_PARTICLES; i++){
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
    normal_distribution<long double> distribution(0, sqrt(SDESOLVER_DT));
    return SDESOLVER_RANDOM_COEFFICIENT * distribution(generator);
}

//split up the loop for sde solving
void loop_for_particles(int start, int end, vector<particle*> particles, long double* coeffs, TreeNode base){
    long double x_original, y_original, x, y, v_x_i, v_y_i, r, v_x_f, v_y_f, force_val, rand_x, rand_y, *gravity;
//    if(start == 0) {
//        cout << get_random(random_coeff, dt) << endl;
//    }
    gravity = (long double*) malloc(2 * sizeof(long double));
    for (int i = start; i < end; ++i) {
        gravity[0] = 0;
        gravity[1] = 0;
        x_original = particles[i]->x;
        y_original = particles[i]->y;
        v_x_i = particles[i]->vx;
        v_y_i = particles[i]->vx;

        r = sqrt(pow(x_original, 2) + pow(y_original, 2));

        force_val = force(r, coeffs);
        rand_x = get_random();
        rand_y = get_random();
        base.calculate_force(particles[i], gravity);
        v_x_f = v_x_i + (force_val * (x_original / r) * SDESOLVER_DT -
                SDESOLVER_CD * v_x_i * SDESOLVER_DT +
                         rand_x + gravity[0] * SDESOLVER_DT) / SDESOLVER_MASS;
        v_y_f = v_y_i + (force_val * (y_original / r) * SDESOLVER_DT -
                SDESOLVER_CD * v_y_i * SDESOLVER_DT +
                         rand_y + gravity[0] * SDESOLVER_DT) / SDESOLVER_MASS;
//        if (i == 0) {
//            cout << y << ", " << (force(r, coeffs, terms) * (y / r) * dt -
//                                  cd * v_y_i * dt +
//                                  0.001 * get_random(random_coeff, dt)) / mass << ", " << v_y_i << endl;
//        }
        particles[i]->vx = v_x_f;
        particles[i]->vy = v_y_f;
        x = x_original + v_x_f * SDESOLVER_DT;
        y = y_original + v_y_f * SDESOLVER_DT;
        if (x > 0.09 || x < -0.09 || y > 0.09 || y < -0.09){
            particles[i]->x = 0;
            particles[i]->y = 0;
        } else {
            particles[i]->x = x;
            particles[i]->y = y;
        }
    }
}

//euler's method
vector<particle*> solve_sde(long double* positions[], long double coeffs[], long double** interp){
    auto start = chrono::high_resolution_clock::now();

    ofstream outfile;
    outfile.open("position.csv");

    TreeNode base(-0.05, -0.05, 0.1, interp);
    struct particle *p;
    vector<particle*> particles;
    cout << SDESOLVER_PARTICLES << endl;
    for(int i = 0; i < SDESOLVER_PARTICLES; i++) {
        p = new struct particle(positions[i][0], positions[i][1]);
        particles.push_back(p);
    }
    float length = ((float) SDESOLVER_PARTICLES) / SDESOLVER_N_THREADS; //number of particles for each thread to process
    for (int t = 0; t < SDESOLVER_N; ++t) {
        base.clear();
        for(int i = 0; i < SDESOLVER_PARTICLES; i++){
            cout << i << endl;
            base.insert(particles[i]);
        }
        thread threads[SDESOLVER_N_THREADS];
        for(int i = 0; i < SDESOLVER_N_THREADS; i++){
            threads[i] = thread(loop_for_particles, (i * length),  ((i + 1) * length), particles, coeffs, base);
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
    return particles;
}

struct arrsize{
    void *a;
    unsigned n;
};

struct arrsize* readforce(){ // double**
    FILE* f = fopen(SDESOLVER_BESSELINTERP, "r");
    unsigned nlines = 0,i;
    char c;
    long double **forcedata;
    forcedata = (long double**) malloc(2*sizeof(long double *));
    for(c=getc(f); c!=EOF; c=getc(f)){
        if(c=='\n')++nlines;
    }
    if(c!='\n')++nlines;
    rewind(f);
    forcedata[0] = (long double*) malloc((nlines+1)*sizeof(long double));
    forcedata[1] = (long double*) malloc((nlines+1)*sizeof(long double));
    //printf("%d\n",nlines);
    for(i=0;i<nlines;++i){
        fscanf(f,"%Lf,%Lf",&forcedata[0][i],&forcedata[1][i]);
        //printf("%d: %.9e,%.9e\n",i,forcedata[0][i],forcedata[1][i]);
    }
    struct arrsize *r = (arrsize*) malloc(sizeof(struct arrsize));
    r->n = nlines;
    r->a = forcedata;
    return r;
}

long double** gen_lin_interp(long double **data,unsigned n){
    long double *x=data[0], *y=data[1],
            **coeff;
    unsigned int i;
    coeff = (long double**) malloc(3*sizeof(long double *));
    coeff[0] = (long double*) malloc(n*sizeof(long double *)); // b
    coeff[1] = (long double*) malloc(n*sizeof(long double *)); // a (bx+a linear interp)
    coeff[2] = (long double*) malloc(n*sizeof(long double *)); // the boundary
    for(i=0;i<n-1;++i)coeff[1][i]=(y[i+1]-y[i])/(x[i+1]-x[i]);
    for(i=0;i<n-1;++i)coeff[0][i]=y[i+1]-coeff[1][i]*x[i+1];
    for(i=0;i<n;++i)coeff[2][i]=x[i];
    coeff[0][n-1]=0;
    coeff[1][n-1]=0;
    //for(i=0;i<n;++i)printf("%d: %.9e,%.9e,%.9e\n",i,coeff[0][i],coeff[1][i],coeff[2][i]);
    return coeff;
}

int main() {
    auto *coeffs = (long double*) malloc(SDESOLVER_TERMS * sizeof(long double));
    //import polynomial coefficients
    if(import_coeffs(coeffs)){
        return 1;
    }
    auto *initial_data = (double*) malloc(SDESOLVER_INITIAL_DATA_LENGTH * sizeof(double));
    //import initial data
    if(import_initial_data(initial_data)){
        return 1;
    }

    //initialize 2d array for positions
    auto** positions = (long double**) malloc(SDESOLVER_PARTICLES * sizeof(long double*));
    for(int i = 0; i < SDESOLVER_PARTICLES; i++){
        positions[i] = (long double*) malloc(2 * sizeof(long double));
    }
    //initialize 2d array for collision detection
    auto** mesh = (bool**) malloc(SDESOLVER_MESH_FINENESS * sizeof(bool*));
    for(int i = 0; i < SDESOLVER_MESH_FINENESS; i++){
        mesh[i] = (bool*) malloc(SDESOLVER_MESH_FINENESS * sizeof(bool));
    }
    //initialize to false
    for (int i = 0; i < SDESOLVER_MESH_FINENESS; i++) {
        for (int j = 0; j < SDESOLVER_MESH_FINENESS; j++) {
            mesh[i][j] = false;
        }
    }
    //initialize positions
    generate_initial_positions(initial_data, positions, mesh);
    free(initial_data); // free up memory since its not needed any more

    //allocate memory for initial velocities array
    auto** velocities = (long double**) malloc(SDESOLVER_PARTICLES * sizeof(long double*));
    for(int i = 0; i < SDESOLVER_PARTICLES; i++){
        velocities[i] = (long double*) malloc(2 * sizeof(long double));
    }
    //set initial velocities to 0
    for (int i = 0; i < SDESOLVER_PARTICLES; i++) {
        velocities[i][0] = 0;
        velocities[i][1] = 0;
    }

    struct arrsize* forcedata = readforce();
    long double** bessel_interp = gen_lin_interp((long double**) forcedata->a, forcedata->n);

    vector<particle*> particles = solve_sde(positions, coeffs, bessel_interp);

    ofstream outfile;
    outfile.open("final_positions.csv");
    for (int i = 0; i < SDESOLVER_PARTICLES; ++i) {
        outfile << particles[i]->x << "," << particles[i]->y << endl;
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