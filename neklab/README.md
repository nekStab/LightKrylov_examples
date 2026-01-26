# JOSS LightKrylov - neklab example

This folder contains the source code for the numerical examples presented in the reference below, showcasing the integration of the [LightKrylov](https://github.com/nekStab/LightKrylov) abstract linear algebra toolbox with the massively parallel open-source CFD code [Nek5000](https://github.com/Nek5000/Nek5000).

# Example list
* `Newton-GMRES`: Newton-Krylov iteration to find the steady fixed-point of the nonlinear Navier-Stokes equations for the supercritical cylinder flow at `Re = 100`.
* `eigs`: Krylov-Schur eigenvalue iteration to find the leading eigenpair of the exponential propagator of the linearized Navier-Stokes operator around the steady fixed-point of the supercritical cylinder flow at `Re = 100`.

## Usage

All examples contain necessary files for code compilation and running the example. Some additional files can be directly recreated using standard `Nek5000` tools. Moreover, to reduce a number of binary files in the repository, we do not include multiple copies of the mesh file `1cyl.re2`. This mesh file must be linked or copied into the run directory prior to execution.

Files provided with each example.
* setup source file `###.usr`.
* runtime parameters file `###.par`.
* required `SIZE` file containing definitions of static arrays dimensions.
* Solution fields for different Reynolds numbers `BF_1cyl*.fld`. These fields are used either as initial conditions for the Newton-GMRES fixed-point iteration or as baseflows for the eigenvalue computation.

To compile the code:

* Clone  `Nek5000` and apply the `LightKrylov`-specific changes. This is most easily achieved by running the `Nek5000_setup.sh` script in the `neklab` root directory
  ```bash
  bash Nek5000_setup.sh
  ```
  This script will clone a pinned version of the `Nek5000` repository that is known to work for this example. Furthermore, it can install some required dependencies and ensure the necessary `Nek5000` tools are compiled to proceed. It can also add some required environment variables to your `.bashrc`. If this is not desired, they can be exported directly
  ```bash
  export NEKLAB_SOURCE_ROOT=$(pwd)
  export NEK_SOURCE_ROOT="$NEKLAB_SOURCE_ROOT/Nek5000"
  export PATH=$NEK_SOURCE_ROOT/bin:$PATH"
  ```

* Clone `LightKrylov` install it on your system. This is most easily achieved by running the `LightKrylov_setup.sh` script in the `neklab` root directory
  ```bash
  bash LightKrylov_setup.sh
  ```
  We recommend running the test suite to check that everything works correctly.

* Navigate to the folder of the desired example and build it using the `app/makeneklab` script from the `neklab` root directory. This script requires the above mentioned environment variables, so ensure that they are correctly defined in your shell.

To run the case:

* Navigate to the folder of the example.
* Copy or link the mesh file `1cyl.re2` from the examples directory.
* Generate the processor map file `1cyl.ma2` using the tool `genmap` distributed together with `Nek5000`. Running this tool on the command-line will ask for the file name, where you should enter just `1cyl`.
* Run the code using the parallel executable provided in `Nek5000/bin`.

To run on 12 cores in the foreground with `nek5000` logfile output to stdout:
```
nekmpi 1cyl 12
```
To run on 12 cores in the background with `nek5000` logfile output to `logfile`:
```
nekbmpi 1cyl 12
```
The `nek5000` logfile contains runtime information for each timestep. The information about the progress of the algorithms from `LightKrylov` are output to `lightkrylov.log`.

**Note**: The executable that is compiled with the default settings can only be run on 12 cores (or more) since the static array dimensions are chosen appropriately. To run on fewer cores, change the line
```
parameter (lpmin=12)              ! min number of MPI ranks
```
in the `SIZE` file and recompile.

# References

Kern et al. (2025). LightKrylov: Lightweight implementation of Krylov subspace techniques in modern Fortran. Journal of Open Source Software, 2
¿VOL? (¿ISSUE?), ¿PAGE? https://doi.org/10.xxxxxx/draft.
