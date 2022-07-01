#ifndef SDESOLVER_SDE_SOLVER_CONSTANTS_H
#define SDESOLVER_SDE_SOLVER_CONSTANTS_H

#define SDESOLVER_FORCE_DATA_FILENAME "forcedata.csv" // coefficients for force polynomial
#define SDESOLVER_TERMS  26 // number of terms in force polynomial
#define SDESOLVER_INITIAL_DATA_LENGTH  1000000 // number of initial r values
#define SDESOLVER_INITIAL_DATA_FILENAME  "initial_data.csv" // initial r values
#define SDESOLVER_PARTICLES  10000  // number of particles to simulate
#define SDESOLVER_MESH_FINENESS  12000 // dimensions of mesh (MESH_FINENESS * MESH_FINENESS)
#define SDESOLVER_N 800000  // number of timesteps
#define SDESOLVER_COLLISION_TOLERANCE 0.1  // smaller is more accurate collision detection
#define SDESOLVER_N_THREADS  7

#define SDESOLVER_VISCOSITY  (0.001) // dynamic viscosity of water
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

#endif //SDESOLVER_SDE_SOLVER_CONSTANTS_H
