import os
import math
import threading

INITIAL_DATA_LENGTH = 1000000  # number of initial r values
INITIAL_DATA_FILENAME = "initial_data.csv"  # initial r values
PARTICLES = 100000  # number of particles to simulate
MESH_FINENESS = 12000  # dimensions of mesh (MESH_FINENESS * MESH_FINENESS)
DT = 0.05  # size of timesteps
COLLISION_TOLERANCE = 0.05  # smaller is more accurate collision detection
N_THREADS = 7

VISCOSITY = 0.005  # possibly fitted dynamic viscosity of water
RADIUS = 5e-5  # radius of particle
DENSITY = 8940  # density of particles
MASS = ((4.0/3) * DENSITY * math.pi * RADIUS ** 3)  # mass of the particle
CD = (6 * math.pi * VISCOSITY * RADIUS)  # stokes drag

TEMPERATURE = (25 + 273.15)  # temperature
KB = 1.38064852e-23  # boltzmann's #defineant
RANDOM_COEFFICIENT = (2 * KB * TEMPERATURE / CD) ** 0.5  # coefficient in front of the dW term

CONTAINER_RADIUS = 0.0925  # radius of the container

T_START = 0.0
T_END = 1000.0
MAX_ANGLE = 0.5

BESSELINTERP = "besselinterp.csv"

EXPORT_FREQUENCY = 20


def write_constants():
    with open("sde_solver_constants.h", "w+") as file:
        file.write(f"""#ifndef SDESOLVER_SDE_SOLVER_CONSTANTS_H
#define SDESOLVER_SDE_SOLVER_CONSTANTS_H

#define SDESOLVER_FORCE_DATA_FILENAME "forcedata.csv" // coefficients for force polynomial
#define SDESOLVER_TERMS  26 // number of terms in force polynomial
#define SDESOLVER_INITIAL_DATA_LENGTH  {INITIAL_DATA_LENGTH} // number of initial r values
#define SDESOLVER_INITIAL_DATA_FILENAME  "{INITIAL_DATA_FILENAME}" // initial r values
#define SDESOLVER_N  {int((T_END - T_START) / DT)}  // number of timesteps
#define SDESOLVER_PARTICLES  {PARTICLES}  // number of particles to simulate
#define SDESOLVER_MESH_FINENESS  {MESH_FINENESS} // dimensions of mesh (MESH_FINENESS * MESH_FINENESS)
#define SDESOLVER_COLLISION_TOLERANCE {COLLISION_TOLERANCE}  // smaller is more accurate collision detection
#define SDESOLVER_N_THREADS  {N_THREADS}  // the number of threads to use

#define SDESOLVER_VISCOSITY  {VISCOSITY} // dynamic viscosity of water
#define SDESOLVER_RADIUS  {RADIUS} // radius of particle
#define SDESOLVER_DENSITY  {DENSITY} // density of particles
#define SDESOLVER_MASS  {MASS}
#define SDESOLVER_CD  {CD} // stokes drag

#define SDESOLVER_TEMPERATURE  {TEMPERATURE} // temperature
#define SDESOLVER_KB  {KB} // boltzmann's #defineant
#define SDESOLVER_RANDOM_COEFFICIENT  {RANDOM_COEFFICIENT}  //coefficient in front of the dW term

#define SDESOLVER_CONTAINER_RADIUS {CONTAINER_RADIUS} // radius of the container

#define SDESOLVER_T_START  {T_START}
#define SDESOLVER_T_END  {T_END}
#define SDESOLVER_MAX_ANGLE {MAX_ANGLE}

#define SDESOLVER_DT {DT}

#define SDESOLVER_BESSELINTERP "{BESSELINTERP}"

#define SDESOLVER_EXPORT_FREQUENCY {EXPORT_FREQUENCY}

#define SDESOLVER_RK_STAGES 1
const double COEFFICIENTS[][SDESOLVER_RK_STAGES] = {{
        {{1}}
}};

#endif //SDESOLVER_SDE_SOLVER_CONSTANTS_H
""")


def run_programs(folders):
    for folder in folders:
        print(f"Running solver in {folder}...")
        os.system(f"cd output/{folder} && ./SDESolver")


recompile = input("Would you like to recompile the files? [y/n] ")

if recompile == "y":
    os.system("rm -r output")
    os.system("mkdir output")
    for n in range(5000, 100000+1, 5000):
        print(f"Compiling with {n} particles...")

        PARTICLES = n
        write_constants()

        os.system("rm ./SDESolver")
        os.system("cmake .")
        os.system("make")

        # print("\nRunning program...")
        # os.system("./SDESolver")

        print("\nShifting exectuables + data files...")
        os.system(f"mkdir output/{n}_particles")
        for file in ["SDESolver", "forcedata.csv", "initial_data.csv", "besselinterp.csv"]:
            os.system(fr"cp {file} output/{n}_particles")

        print("Done!\n")

# Now, run them (run 4 in parallel to use more CPU)
threading.Thread(target=lambda: run_programs([f"{n}_particles" for n in range(5000, 100000+1, 20000)])).start()
threading.Thread(target=lambda: run_programs([f"{n}_particles" for n in range(10000, 100000+1, 20000)])).start()
threading.Thread(target=lambda: run_programs([f"{n}_particles" for n in range(15000, 100000+1, 20000)])).start()
threading.Thread(target=lambda: run_programs([f"{n}_particles" for n in range(20000, 100000+1, 20000)])).start()
