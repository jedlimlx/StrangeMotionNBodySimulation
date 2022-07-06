#ifndef SDESOLVER_SDE_SOLVER_CONSTANTS_H
#define SDESOLVER_SDE_SOLVER_CONSTANTS_H

#define SDESOLVER_FORCE_DATA_FILENAME "forcedata.csv" // coefficients for force polynomial
#define SDESOLVER_TERMS  26 // number of terms in force polynomial
#define SDESOLVER_INITIAL_DATA_LENGTH  1000000 // number of initial r values
#define SDESOLVER_INITIAL_DATA_FILENAME  "initial_data.csv" // initial r values
#define SDESOLVER_N  20000  // number of timesteps
#define SDESOLVER_PARTICLES  100000  // number of particles to simulate
#define SDESOLVER_MESH_FINENESS  12000 // dimensions of mesh (MESH_FINENESS * MESH_FINENESS)
#define SDESOLVER_COLLISION_TOLERANCE 0.05  // smaller is more accurate collision detection
#define SDESOLVER_N_THREADS  7  // the number of threads to use

#define SDESOLVER_VISCOSITY  0.04 // dynamic viscosity of water
#define SDESOLVER_RADIUS  5e-05 // radius of particle
#define SDESOLVER_DENSITY  8940 // density of particles
#define SDESOLVER_MASS  4.680973053848793e-09
#define SDESOLVER_CD  3.769911184307752e-05 // stokes drag

#define SDESOLVER_TEMPERATURE  298.15 // temperature
#define SDESOLVER_KB  1.38064852e-23 // boltzmann's #defineant
#define SDESOLVER_RANDOM_COEFFICIENT  1.4777753265340708e-08  //coefficient in front of the dW term

#define SDESOLVER_CONTAINER_RADIUS 0.0925 // radius of the container

#define SDESOLVER_T_START  0.0
#define SDESOLVER_T_END  1000.0
#define SDESOLVER_MAX_ANGLE 0.5

#define SDESOLVER_DT 0.05

#define SDESOLVER_BESSELINTERP "besselinterp.csv"

#define SDESOLVER_EXPORT_FREQUENCY 20

#define SDESOLVER_RK_STAGES 1
const double COEFFICIENTS[][SDESOLVER_RK_STAGES] = {
        {1}
};

#endif //SDESOLVER_SDE_SOLVER_CONSTANTS_H
