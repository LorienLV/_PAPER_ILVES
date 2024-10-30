# ILVES

In molecular dynamics, the time step can be incremented through the imposition of constraints on internal degrees of freedom, which raises the total simulated time for equal computational effort. This enables researchers to simulate a wider collection of phenomena of interest, which frequently take place at large time scales. Customarily, the bond length between pair of atoms is constrained, although other constrains such as angles and dihedral angles can be imposed to further increase the time step.

The most popular constraint solvers are SHAKE and LINCS. SHAKE converges locally and linearly and operates serially, limiting parallel performance. LINCS, along with its parallel variant P-LINCS, converges linearly at best and, in certain cases, fails to converge. Additionally, LINCS is suited only for constraining bond lengths. As a result, achieving machine-precision solutions with these methods is time-consuming, and constraining angles or dihedral angles efficiently remains a challenge.

ILVES is a new constraint solver that solves the same system of differential-algebraic equations as SHAKE using Newton's method. ILVES converges quadratically and can solve nonlinear constraint equations in parallel, providing faster and more accurate solutions than current methods. Notably, ILVES is the first constraint solver to enable applying angle constraints in parallel.

This repository contains two versions of the ILVES algorithm, integrated into the GROMACS (v2021.0) molecular dynamics package:
- **ILVES**: The base version, solving the same equations as SHAKE.
- **ILVES-FAST**: A Quasi-Newton variation that achieves up to twice the performance of ILVES by solving a modified set of equations.

Additionally, this repository includes the molecular systems evaluated in the paper "Accurate and Efficient Constrained Molecular Dynamics through Direct Solvers: The ILVES Algorithm."

## ILVES in GROMACS

The `gromacs-2021` folder contains the 2021.0 version of GROMACS with the addition of the two versions of ILVES. The original P-LINCS solver has been modified to allow users to set a custom tolerance for the algorithm, identically to SHAKE's `shake-tol` parameter in the `.mdp` file. To reach the chosen tolerance, P-LINCS' repeats its correction phase, controlled by the `lincs-iter` parameter in the original implementation, until the desired error is met.

### Installation

To install GROMACS with ILVES, follow the standard [GROMACS installation guide](https://manual.gromacs.org/documentation/2021/install-guide/index.html). ILVES functionality will be available out of the box.

**Note**: It’s recommended to compile with `GCC v10` as later versions may cause compilation errors in parts of the package.

Example of an MPI/OpenMP installation in double-precision mode:
```
mkdir build
cd build

cmake ../gromacs-2021 \
-DGMX_BUILD_OWN_FFTW=on \
-DGMX_MPI=on \
-DGMX_DOUBLE=on \
-DCMAKE_BUILD_TYPE=Release \
-DREGRESSIONTEST_DOWNLOAD=off \
-DBUILD_TESTING=off \
-DGMX_GPU=off \
-DGMX_OPENMP=on \
-DGMX_X11=off \
-DGMXAPI=off

make
make install
```

### Usage

To use ILVES as the constraint solver in your simulation, specify the following parameters in the `.mdp` file:

```
constraint-algorithm=ilves  ; Use ILVES
; OR
constraint-algorithm=ilvesf ; Use ILVES-FAST
ilves-tol=0.0001            ; Relative tolerance for ILVES, analogous to shake-tol in SHAKE.
ilves-mpi-neigh=4"          ; With domain decomposition, the cell size is limited by the distance spanned by ilves-mpi-neigh+1 constraints, similar to lincs-order.
```

Our modified version of P-LINCS introduces changes to the algorithm's parameters:
```
lincs-tol=0.0001  ; Relative tolerance for LINCS, analogous to shake-tol and ilves-tol. This is a new parameter.
lincs-iter=100    ; Maximum allowed iterations in the LINCS correction phase (fixed number of iterations in the original P-LINCS).
```

## Reproduce our Results

The `simulations` folder includes the molecular systems evaluated in the ILVES manuscript. Additionally, it provides scripts to generate necessary files for each system's production run, along with a script to perform production runs with different tolerances, parallel configurations, and constraint solvers.

- `simulations/setup.sh`: Sets up the environment for running our simulations.
- `simulations/ffs`: Contains the force fields used in our simulations.
- `simulations/systems`: Contains a subfolder for each molecular system studied. Each system folder includes a setup script, directories with required input files (e.g., `.mdp` and `.pdb` files) for simulation preparation, and automation scripts for generating production-ready files. Additionally, the `files-for-prod` folder provides pre-generated input files for immediate use in production runs.
- `simulations/prod_run.sh`: Executes production runs with various tolerances, parallel configurations, and constraint solvers using the generated production files.

### Generating Production Files for a System

Each system folder includes a sequence of scripts, numbered from 1 to N, to automatically generate the required production files:

```
cd simulations
source setup.sh

cd systems/DESIRED_SYSTEM
source setup.sh

bash 1_ONE_SCRIPT
...
bash N_extract_prod_files.sh #  The files-for-prod folder will contain all required files for a production run.
```

You may skip the setup and directly use the pre-generated files in `files-for-prod`.

### Running a Production Simulation

To execute a production run, you can use `simulations/prod_run.sh`, which automates `.tpr` file generation for various tolerances, parallel setups, and solvers. Check out the script to verify which simulations it will execute. Run the script as follows:

```
source ./simulations/setup.sh
source ./simulations/systems/DESIRED_SYSTEM/setup.sh

bash prod_run.sh
```

## Cite Us

TODO