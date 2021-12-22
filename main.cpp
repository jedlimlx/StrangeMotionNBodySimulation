#include <iostream>
#include <fstream>
#include <vector>

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
    return 0;
}

int main() {
    const string FORCE_COEFFS_FILENAME = "forcecoeffs.csv";
    const int TERMS = 26;
    const int INITIAL_DATA_LENGTH = 750357;
    const string INITIAL_DATA_FILENAME = "initial_data.csv";
    const int PARTICLES = 20000;

    long double coeffs[TERMS];
    //import polynomial coefficients
    if(import_coeffs(FORCE_COEFFS_FILENAME, TERMS, coeffs)){
        return 1;
    }

    double initial_data[INITIAL_DATA_LENGTH];
    //import initial data
    if(import_initial_data(INITIAL_DATA_FILENAME, INITIAL_DATA_LENGTH, initial_data)){
        return 1;
    }

    //randomly generate initial positions
    long double positions[PARTICLES][2];
    long double initial_r[PARTICLES];

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