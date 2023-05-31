#ifndef SDESOLVER_SDE_SOLVER_CONSTANTS_H
#define SDESOLVER_SDE_SOLVER_CONSTANTS_H

#define SDESOLVER_FORCE_DATA_FILENAME "forcegraphite6mm.csv" // coefficients for force polynomial
#define SDESOLVER_TERMS  26 // number of terms in force polynomial
#define SDESOLVER_INITIAL_DATA_LENGTH  1000000 // number of initial r values
#define SDESOLVER_INITIAL_DATA_FILENAME  "initial_data.csv" // initial r values
#define SDESOLVER_N  50000  // number of timesteps
#define SDESOLVER_PARTICLES  10000  // number of particles to simulate
#define SDESOLVER_MESH_FINENESS  12000 // dimensions of mesh (MESH_FINENESS * MESH_FINENESS)
#define SDESOLVER_COLLISION_TOLERANCE 0.05  // smaller is more accurate collision detection
#define SDESOLVER_N_THREADS  7  // the number of threads to use

#define SDESOLVER_VISCOSITY  0.005 // dynamic viscosity of water
#define SDESOLVER_RADIUS  4e-05 // radius of particle
#define SDESOLVER_DENSITY  2260 // density of particles
#define SDESOLVER_MASS  6.058666152203036e-10
#define SDESOLVER_CD  3.7699111843077525e-06 // stokes drag

#define SDESOLVER_TEMPERATURE  298.15 // temperature
#define SDESOLVER_KB  1.38064852e-23 // boltzmann's #defineant
#define SDESOLVER_RANDOM_COEFFICIENT  4.673135901846724e-08  //coefficient in front of the dW term

#define SDESOLVER_CONTAINER_RADIUS 0.0925 // radius of the container

#define SDESOLVER_T_START  0.0
#define SDESOLVER_T_END  1000.0
#define SDESOLVER_MAX_ANGLE 0.5

#define SDESOLVER_DT 0.05

#define SDESOLVER_BESSELINTERP "cheeriosgraphite.csv"

#define SDESOLVER_EXPORT_FREQUENCY 1

#define SDESOLVER_RK_STAGES 1
const double COEFFICIENTS[][SDESOLVER_RK_STAGES] = {
        {1}
};

#endif //SDESOLVER_SDE_SOLVER_CONSTANTS_H
