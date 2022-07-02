#ifndef SDESOLVER_SDE_SOLVER_CONSTANTS_H
#define SDESOLVER_SDE_SOLVER_CONSTANTS_H

#define SDESOLVER_FORCE_DATA_FILENAME "forcedata.csv" // coefficients for force polynomial
#define SDESOLVER_TERMS  26 // number of terms in force polynomial
#define SDESOLVER_INITIAL_DATA_LENGTH  1000000 // number of initial r values
#define SDESOLVER_INITIAL_DATA_FILENAME  "initial_data.csv" // initial r values
#define SDESOLVER_PARTICLES  10000  // number of particles to simulate
#define SDESOLVER_MESH_FINENESS  12000 // dimensions of mesh (MESH_FINENESS * MESH_FINENESS)
#define SDESOLVER_N 15000000  // number of timesteps
#define SDESOLVER_COLLISION_TOLERANCE 0.1  // smaller is more accurate collision detection
#define SDESOLVER_N_THREADS  7

#define SDESOLVER_VISCOSITY  (0.04) // dynamic viscosity of water
#define SDESOLVER_RADIUS  (5e-5) // radius of particle
#define SDESOLVER_DENSITY  (2266) // density of particles
#define SDESOLVER_MASS  ((4.0/3) * SDESOLVER_DENSITY * M_PI * pow(SDESOLVER_RADIUS, 3))
#define SDESOLVER_CD  (6 * M_PI * SDESOLVER_VISCOSITY * SDESOLVER_RADIUS) // stokes drag
#define SDESOLVER_TEMPERATURE  (25 + 273.15) // temperature
#define SDESOLVER_KB  (1.38064852e-23) // boltzmann's #defineant
#define SDESOLVER_RANDOM_COEFFICIENT  (sqrt(2 * SDESOLVER_KB * SDESOLVER_TEMPERATURE / SDESOLVER_CD)) //coefficient in front of the dW term
#define SDESOLVER_CONTAINER_RADIUS 0.0925 // radius of the container
#define SDESOLVER_T_START  0.0
#define SDESOLVER_T_END  1000.0
#define SDESOLVER_CAPL_LENGTH 368.3751861591355
#define SDESOLVER_MAX_ANGLE 0.5
#define SDESOLVER_COEFF 166.649
#define SDESOLVER_DT ((SDESOLVER_T_END - SDESOLVER_T_START) / SDESOLVER_N)
#define SDESOLVER_BESSELINTERP "besselinterp.csv"

/*
#define SDESOLVER_RK_STAGES 7
const double COEFFICIENTS[][SDESOLVER_RK_STAGES] = {
        {0.2},
        {3.0 / 40.0, 9.0 / 40.0},
        {44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0},
        {19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0},
        {9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0},
        {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0},
        {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0}
};
*/

/*
#define SDESOLVER_RK_STAGES 6
const double COEFFICIENTS[][SDESOLVER_RK_STAGES] = {
        {0.2},
        {3.0 / 40.0, 9.0 / 40.0},
        {3.0 / 10.0, -9.0 / 10.0, 6.0 / 5.0},
        {-11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0},
        {1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0, 253.0 / 4096.0},
        {37.0 / 378.0, 0.0, 250.0 / 612.0, 125.0 / 594.0, 0.0, 512.0 / 1771.0},
};
 */


#define SDESOLVER_RK_STAGES 6
const double COEFFICIENTS[][SDESOLVER_RK_STAGES] = {
        {0.25},
        {3.0 / 32.0, 9.0 / 32.0},
        {1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0},
        {439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0},
        {-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0},
        {16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0},
};

/*
#define SDESOLVER_RK_STAGES 3
const double COEFFICIENTS[][SDESOLVER_RK_STAGES] = {
        {1.0},
        {0.25, 0.25},
        {1.0/6.0, 1.0/6.0, 2.0/3.0}
};
*/

/*
#define SDESOLVER_RK_STAGES 4
const double COEFFICIENTS[][SDESOLVER_RK_STAGES] = {
        {0.5},
        {0, 0.5},
        {0, 0, 1.0},
        {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0}
};
*/

/*
#define SDESOLVER_RK_STAGES 4
const double COEFFICIENTS[][SDESOLVER_RK_STAGES] = {
        {1.0 / 3.0},
        {-1.0 / 3.0, 1.0},
        {1.0, -1.0, 1.0},
        {1.0 / 8.0, 3.0 / 8.0, 3.0 / 8.0, 1.0 / 8.0}
};
*/

#endif //SDESOLVER_SDE_SOLVER_CONSTANTS_H
