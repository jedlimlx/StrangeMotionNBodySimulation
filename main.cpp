#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <thread>
#include <chrono>
#include "sde_solver_constants.h"
#include "TreeNode.h"

using namespace std;

// polynomial for the force
long double force(long double r, const long double coeffs[]) {
    long double sum = 0;
    for (int i = 0; i < SDESOLVER_TERMS; i++){
        sum += coeffs[i] * pow(r, i);
    }
    return sum;
}

// reads in the initial positions
int import_initial_data(double positions[]) {
    ifstream infile;
    infile.open(SDESOLVER_INITIAL_DATA_FILENAME);

    // error checking
    if (!infile.is_open()) {
        cout << "Failed to open file" << endl;
        return 1;
    }

    // read in the data
    string word;
    for(int i = 0; i < SDESOLVER_INITIAL_DATA_LENGTH; i ++){
        getline(infile, word);
        positions[i] = stod(word);
    }

    infile.close();
    return 0;
}

// initialize positions of particles
int generate_initial_positions(const double initial_data[], long double* positions[], bool**mesh) {
    default_random_engine generator(time(NULL));
    uniform_real_distribution<long double> angle_distribution(0, 2 * M_PI);
    uniform_int_distribution<int> index_distribution(0,  SDESOLVER_INITIAL_DATA_LENGTH - 1);

    long double r, theta, x, y;
    int index;
    for (int i = 0; i < SDESOLVER_PARTICLES; i++) {
        index = index_distribution(generator);
        r = initial_data[index];

        theta = angle_distribution(generator);
        x = r * cos(theta);
        y = r * sin(theta);

        while (mesh[6000 + ((int) (x / SDESOLVER_RADIUS))][6000 + ((int) (y / SDESOLVER_RADIUS))]) {
            theta = angle_distribution(generator);
            x = r * cos(theta);
            y = r * sin(theta);
        }

        if (isnan(x) || isnan(y)) {
            cout << "test" << endl;
        }

        positions[i][0] = x;
        positions[i][1] = y;
    }

    return 0;
}

// rng time
long double get_random(){
    static default_random_engine generator(time(NULL));
    normal_distribution<long double> distribution(0, sqrt(SDESOLVER_DT));
    return SDESOLVER_RANDOM_COEFFICIENT * distribution(generator);
}

long double eval_interp1(long double **interp, long double point) {
    long double *c0 = interp[1], *c1 = interp[0], *bd = interp[2];

    unsigned start = 0, end = 920, i;

    while (start + 1 < end) {
        i = (start + end) / 2;
        // cout << start << " " << end << " " << i << " " << endl;

        if (point > bd[i]) start = i;
        else if (point < bd[i]) end = i;
    }

     /*
    unsigned i = 0;
    while (c0[i] < INFINITY) {
        if (point < bd[i+1]) break;
        ++i;
    } // TODO: use some binary search
     */

    if (isnan(c0[i])) --i;
    return c0[i] + c1[i] * point;
}

// split up the loop for sde solving
void loop_for_particles(int start, int end, struct particle** particles, long double** force_interp,
        TreeNode base, struct particle*** k, long double** memo_gravity, int stage, bool recalculate_gravity) {
    long double x_original, y_original, x, y, v_x_i, v_y_i, r, ax, ay, force_val, *gravity, total_a, total_v, cd;
    gravity = (long double*) malloc(2 * sizeof(long double));
    for (int i = start; i < end; ++i) {
        gravity[0] = 0;
        gravity[1] = 0;

        x_original = particles[i]->x;
        y_original = particles[i]->y;
        v_x_i = particles[i]->vx;
        v_y_i = particles[i]->vy;

        r = sqrt(pow(x_original, 2) + pow(y_original, 2));
        force_val = eval_interp1(force_interp, r);

        if (isnan(force_val)) {
            cout << "jed " << r << endl;
            exit(1);
        }

        if (recalculate_gravity) {
            base.calculate_force(particles[i], gravity);
            memo_gravity[i][0] = gravity[0];
            memo_gravity[i][1] = gravity[1];
        } else {
            gravity[0] = memo_gravity[i][0];
            gravity[1] = memo_gravity[i][1];
        }

        // cout << force_val * (x_original / r) << " " << force_val * (y_original / r) << " "
        // << gravity[0] << " " << gravity[1] << " " << endl;

        ax = force_val * (x_original / r) + gravity[0];
        ay = force_val * (y_original / r) + gravity[1];

        // Updating things in k
        k[stage][i]->vx = v_x_i;
        k[stage][i]->vy = v_y_i;
        k[stage][i]->ax = ax / SDESOLVER_MASS;
        k[stage][i]->ay = ay / SDESOLVER_MASS;

        // Updating for the next time-step
        cd = SDESOLVER_CD / SDESOLVER_MASS;
        if (true || stage < SDESOLVER_RK_STAGES - 1) {
            // Update x-direction first
            total_a = 0;
            for (int j = 0; j < stage + 1; j++)
                total_a += COEFFICIENTS[stage][j] * k[j][i]->ax;

            total_v = 0;
            for (int j = 0; j < stage + 1; j++)
                total_v += COEFFICIENTS[stage][j] * k[j][i]->vx;

            x = exp(-cd * SDESOLVER_DT) / pow(cd, 2) * (total_a -
                    cd * total_v + exp(cd * SDESOLVER_DT) * (total_a *
                    (cd * SDESOLVER_DT - 1) + cd * (total_v + cd * x_original)));
            particles[i]->vx = exp(-cd * SDESOLVER_DT) / cd * (total_a *
                    (exp(cd * SDESOLVER_DT) - 1) + cd * total_v);

            // Update y-direction now
            total_a = 0;
            for (int j = 0; j < stage + 1; j++)
                total_a += COEFFICIENTS[stage][j] * k[j][i]->ay;

            total_v = 0;
            for (int j = 0; j < stage + 1; j++)
                total_v += COEFFICIENTS[stage][j] * k[j][i]->vy;

            y = exp(-cd * SDESOLVER_DT) / pow(cd, 2) * (total_a -
                    cd * total_v + exp(cd * SDESOLVER_DT) * (total_a *
                    (cd * SDESOLVER_DT - 1) + cd * (total_v + cd * y_original)));
            particles[i]->vy = exp(-cd * SDESOLVER_DT) / cd * (total_a *
                    (exp(cd * SDESOLVER_DT) - 1) + cd * total_v);

            /*
            particles[i]->vy = v_y_i;
            for (int j = 0; j < stage + 1; j++)
                particles[i]->vy += COEFFICIENTS[stage][j] * k[j][i]->ay * SDESOLVER_DT;

            x = x_original;
            for (int j = 0; j < stage + 1; j++)
                x += COEFFICIENTS[stage][j] * k[j][i]->vx * SDESOLVER_DT;

            y = y_original;
            for (int j = 0; j < stage + 1; j++)
                y += COEFFICIENTS[stage][j] * k[j][i]->vy * SDESOLVER_DT;
            */
        }

        // Check if the particle hit the side of the bowl
        if (sqrt(pow(x, 2) + pow(y, 2)) > SDESOLVER_CONTAINER_RADIUS) {
            particles[i]->vx = 0;
            particles[i]->vy = 0;

            particles[i]->x = (SDESOLVER_CONTAINER_RADIUS - 0.001) / sqrt(pow(x, 2) + pow(y, 2)) * x;
            particles[i]->y = (SDESOLVER_CONTAINER_RADIUS - 0.001) / sqrt(pow(x, 2) + pow(y, 2)) * y;
        } else {
            particles[i]->x = x;
            particles[i]->y = y;
        }
    }

    free(gravity);
}

// prince dormand method :DDD
struct particle** solve_sde(long double* positions[], long double** force_interp, long double** interp) {
    auto start = chrono::high_resolution_clock::now();

    TreeNode base(-0.1, -0.1, 0.2, interp, 0, 0);

    // Stores their original positions
    auto **particles = (particle**) malloc(SDESOLVER_PARTICLES * sizeof(struct particle*));

    // Stores their "newer" positions (temporary)
    auto **particles2 = (particle**) malloc(SDESOLVER_PARTICLES * sizeof(struct particle*));

    // runge-kutta constants
    auto **k = (particle***) malloc(SDESOLVER_RK_STAGES * SDESOLVER_PARTICLES * sizeof(particle**));
    for (int i = 0; i < SDESOLVER_RK_STAGES; i++) k[i] = (particle**) malloc(SDESOLVER_PARTICLES * sizeof(struct particle*));

    // memorise the strength of gravity
    long double** memo_gravity;
    memo_gravity = (long double**) malloc(SDESOLVER_PARTICLES * sizeof(long double *));
    for (int i = 0; i < SDESOLVER_PARTICLES; i++) {
        memo_gravity[i] = (long double*) malloc(2*sizeof(long double));
        memo_gravity[i][0] = 0;
        memo_gravity[i][1] = 0;
    }

    for (int i = 0; i < SDESOLVER_PARTICLES; i++) {
        particles[i] = new struct particle(positions[i][0], positions[i][1]);

        for (int j = 0; j < SDESOLVER_RK_STAGES; j++)
            k[j][i] = new struct particle(0, 0);

        if (isnan(particles[i] -> x)) {
            cout << i << endl;
            exit(1);
        }
    }

    float length = ((float) SDESOLVER_PARTICLES) / SDESOLVER_N_THREADS; // number of particles for each thread to process
    for (int t = 0; t < SDESOLVER_N; ++t) {
        base.clear();
        for (int i = 0; i < SDESOLVER_PARTICLES; i++) {
            base.insert_particle(particles[i]);

            // Making deep copy of particles
            particles2[i] = new struct particle(particles[i]->x, particles[i]->y);
            particles2[i]->vx = particles[i]->vx;
            particles2[i]->vy = particles[i]->vy;
        }

        // Runge-kutta stages
        for (int j = 0; j < SDESOLVER_RK_STAGES; j++) {
            thread threads[SDESOLVER_N_THREADS];
            for (int i = 0; i < SDESOLVER_N_THREADS; i++) {
                threads[i] = thread(loop_for_particles, (i * length), ((i + 1) * length), particles2, force_interp,
                                    base, k, memo_gravity, j, j == 0 && t % 5 == 0);
            }

            for (auto &thread: threads) {
                thread.join();
            }
        }

        // Update everything :D
        for (int i = 0; i < SDESOLVER_PARTICLES; i++) {
            /*
            for (int j = 0; j < SDESOLVER_RK_STAGES; j++)
                particles[i]->vx += COEFFICIENTS[SDESOLVER_RK_STAGES - 1][j] * k[j][i]->ax * SDESOLVER_DT;

            for (int j = 0; j < SDESOLVER_RK_STAGES; j++)
                particles[i]->vy += COEFFICIENTS[SDESOLVER_RK_STAGES - 1][j] * k[j][i]->ay * SDESOLVER_DT;

            for (int j = 0; j < SDESOLVER_RK_STAGES; j++)
                particles[i]->x += COEFFICIENTS[SDESOLVER_RK_STAGES - 1][j] * k[j][i]->vx * SDESOLVER_DT;

            for (int j = 0; j < SDESOLVER_RK_STAGES; j++)
                particles[i]->y += COEFFICIENTS[SDESOLVER_RK_STAGES - 1][j] * k[j][i]->vy * SDESOLVER_DT;
            */

            particles[i]->x = particles2[i]->x;
            particles[i]->y = particles2[i]->y;
            particles[i]->vx = particles2[i]->vx;
            particles[i]->vy = particles2[i]->vy;

            /*
            if (sqrt(pow(particles[i]->x, 2) + pow(particles[i]->y, 2)) > SDESOLVER_CONTAINER_RADIUS) {
                particles[i]->vx = 0;
                particles[i]->vy = 0;

                particles[i]->x = (SDESOLVER_CONTAINER_RADIUS - 0.001) / sqrt(pow(particles[i]->x, 2) +
                        pow(particles[i]->y, 2)) * particles[i]->x;
                particles[i]->y = (SDESOLVER_CONTAINER_RADIUS - 0.001) / sqrt(pow(particles[i]->x, 2) +
                        pow(particles[i]->y, 2)) * particles[i]->y;
            }
            */
        }

        base.resolve_collisions();

        if (t % SDESOLVER_EXPORT_FREQUENCY == 0) {
            cout << chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - start).count() << endl;

            ofstream outfile;
            outfile.open("output/final_positions" + to_string(t) + ".csv");

            // base.get_particles();

            for (int i = 0; i < SDESOLVER_PARTICLES; ++i) {
                outfile << particles[i]->x << "," << particles[i]->y << endl;
            }

            outfile.close();
        }
    }

    return particles;
}

struct arrsize {
    void *a;
    unsigned n;
};

struct arrsize* readforce(char* filename) { // double**
    FILE* f = fopen(filename, "r");
    unsigned nlines = 0, i;
    char c;
    long double **forcedata;
    forcedata = (long double**) malloc(2*sizeof(long double *));
    for (c = getc(f); c != EOF; c = getc(f)){
        if (c == '\n') ++nlines;
    }

    if (c != '\n') ++nlines;
    rewind(f);

    forcedata[0] = (long double*) malloc((nlines+1)*sizeof(long double));
    forcedata[1] = (long double*) malloc((nlines+1)*sizeof(long double));

    for(i = 0; i < nlines; ++i) {
        fscanf(f,"%Lf,%Lf",&forcedata[0][i],&forcedata[1][i]);
    }

    auto *r = (arrsize*) malloc(sizeof(struct arrsize));
    r->n = nlines;
    r->a = forcedata;
    return r;
}

long double** gen_lin_interp(long double **data, unsigned n) {
    long double *x=data[0], *y=data[1], **coeff;
    unsigned int i;
    coeff = (long double**) malloc(3*sizeof(long double *));
    coeff[0] = (long double*) malloc(n*sizeof(long double)); // b
    coeff[1] = (long double*) malloc(n*sizeof(long double)); // a (bx+a linear interp)
    coeff[2] = (long double*) malloc(n*sizeof(long double)); // the boundary

    for (i = 0; i < n - 1; ++i)
        coeff[1][i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);

    for (i = 0; i < n - 1; ++i)
        coeff[0][i] = y[i + 1] - coeff[1][i] * x[i + 1];

    for (i = 0; i < n; ++i)
        coeff[2][i] = x[i];

    coeff[0][n - 1] = 0;
    coeff[1][n - 1] = 0;
    //for(i=0;i<n;++i)printf("%d: %.9e,%.9e,%.9e\n",i,coeff[0][i],coeff[1][i],coeff[2][i]);
    return coeff;
}

int main() {
    cout << SDESOLVER_MASS << endl;

    auto *initial_data = (double*) malloc(SDESOLVER_INITIAL_DATA_LENGTH * sizeof(double));

    // import initial data
    if (import_initial_data(initial_data)) {
        return 1;
    }

    // initialize 2d array for positions
    auto** positions = (long double**) malloc(SDESOLVER_PARTICLES * sizeof(long double*));
    for(int i = 0; i < SDESOLVER_PARTICLES; i++) {
        positions[i] = (long double*) malloc(2 * sizeof(long double));
    }

    // initialize 2d array for collision detection
    auto** mesh = (bool**) malloc(SDESOLVER_MESH_FINENESS * sizeof(bool*));
    for (int i = 0; i < SDESOLVER_MESH_FINENESS; i++) {
        mesh[i] = (bool*) malloc(SDESOLVER_MESH_FINENESS * sizeof(bool));
    }

    // initialize to false
    for (int i = 0; i < SDESOLVER_MESH_FINENESS; i++) {
        for (int j = 0; j < SDESOLVER_MESH_FINENESS; j++) {
            mesh[i][j] = false;
        }
    }

    // initialize positions
    generate_initial_positions(initial_data, positions, mesh);
    free(initial_data); // free up memory since its not needed any more

    struct arrsize* bessel_data = readforce(SDESOLVER_BESSELINTERP);
    long double** bessel_interp = gen_lin_interp((long double**) bessel_data->a, bessel_data->n);

    struct arrsize* force_data = readforce(SDESOLVER_FORCE_DATA_FILENAME);
    long double** force_interp = gen_lin_interp((long double**) force_data->a, force_data->n);

    struct particle** particles = solve_sde(positions, force_interp, bessel_interp);

    ofstream outfile;
    outfile.open("final_positions.csv");
    for (int i = 0; i < SDESOLVER_PARTICLES; ++i) {
        outfile << particles[i]->x << "," << particles[i]->y << endl;
    }

    outfile.close();

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
//import coefficients for the force polynomial