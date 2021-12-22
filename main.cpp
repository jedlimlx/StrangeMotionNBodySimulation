#include <iostream>
#include <fstream>
#include <cmath>
#include <random>

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

int main() {
    const string FORCE_COEFFS_FILENAME = "forcecoeffs.csv";
    const int TERMS = 26;
    const int INITIAL_DATA_LENGTH = 750357;
    const string INITIAL_DATA_FILENAME = "initial_data.csv";
    const int PARTICLES = 20000;
    const int MESH_FINENESS = 3000;

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

    ofstream outfile;
    outfile.open("cppinitial.csv");
    for(int i = 0; i < PARTICLES; i++){
        outfile << sqrt(pow(positions[i][0],2) + pow(positions[i][1] , 2)) << ",";
    }
    outfile.close();

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