//
// Created by qiuzi on 3/1/2022.
//

#ifndef SDESOLVER_SDE_SOLVER_CONSTANTS_H
#define SDESOLVER_SDE_SOLVER_CONSTANTS_H

#define SDESOLVER_FORCE_COEFFS_FILENAME "forcecoeffs.csv" //coefficients for force polynomial
#define SDESOLVER_TERMS  26 //number of terms in force polynomial
#define  SDESOLVER_INITIAL_DATA_LENGTH  399460 //number of initial r values
#define  SDESOLVER_INITIAL_DATA_FILENAME  "initial_positions.csv" //initial r values
#define  SDESOLVER_PARTICLES  20000 //number of particles to simulate
#define  SDESOLVER_MESH_FINENESS  3000 //dimensions of mesh (MESH_FINENESS * MESH_FINENESS)
#define  SDESOLVER_N  100000 //number of timesteps
#define  SDESOLVER_N_THREADS  6

#define   SDESOLVER_VISCOSITY  0.0010518 //dynamic viscosity of water
#define   SDESOLVER_RADIUS  48e-6 //radius of particle
#define   SDESOLVER_DENSITY  2260 //density of particles
#define   SDESOLVER_MASS  (4.0/3) * SDESOLVER_DENSITY * M_PI * pow(SDESOLVER_RADIUS, 3)
#define   SDESOLVER_CD  6 * M_PI * SDESOLVER_VISCOSITY * SDESOLVER_RADIUS //stokes drag
#define   SDESOLVER_TEMPERATURE  28 + 273.15 //temperature
#define   SDESOLVER_KB  1.38064852e-23 //boltzmann's #defineant
#define   SDESOLVER_RANDOM_COEFFICIENT  sqrt(2 * SDESOLVER_KB * SDESOLVER_TEMPERATURE / SDESOLVER_CD) //coefficient in front of the dW term
#define   SDESOLVER_T_START  0.0
#define   SDESOLVER_T_END  1000.0
#define SDESOLVER_CAPL_LENGTH 368.3751861591355
#define SDESOLVER_MAX_ANGLE 0.1
#define SDESOLVER_COEFF 166.649
#define SDESOLVER_DT (SDESOLVER_T_END - SDESOLVER_T_START)/SDESOLVER_N
#define SDESOLVER_QK 0.000553862
#define SDESOLVER_BESSELINTERP "besselinterp.csv"

#endif //SDESOLVER_SDE_SOLVER_CONSTANTS_H
