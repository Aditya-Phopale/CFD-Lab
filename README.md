<div align="center">
  <img width="466" height="492" src="FluidchenLogo.png">
</div>

Fluidchen is a CFD Solver developed for the CFD Lab taught at TUM Informatics, Chair of Scientific Computing in Computer Science.

After forking, use this `README.md` however you want: feel free to remove anything you don't need,
or add any additional details we should know to run the code.

## Working with fluidchen

You will extend this code step-by-step starting from a pure framework to a parallel CFD solver. Please follow these [instructions for work with git and submitting the assignments](docs/first-steps.md).

## Software Requirements

This code is known to work on all currently supported Ubuntu LTS versions (22.04, 20.04, 18.04).
In particular, you will need:

- A recent version of the GCC compiler. Other compilers should also work, but you may need to tweak the CMakeLists.txt file (contributions welcome). GCC 7.4, 9.3, and 11.2 are known to work. See `CMakeLists.txt` and `src/Case.cpp` for some compiler-specific code.
- CMake, to configure the build.
- The VTK library, to generate result files. libvtk7 and libvtk9 are known to work.
- OpenMPI (not for the skeleton, but when you implement parallelization).

Get the dependencies on Ubuntu:

```shell
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install build-essential cmake libvtk7-qt-dev openmpi-bin libopenmpi-dev
```
## Code Theory

For our 2D test case, we assume the fluid is viscous and follows the Navier-Stokes equations. Let,
`U` be the velocity in X-direction,
`V` be the velocity in Y-direction,
`p` be the pressure 

## Application of Boundary Conditions
For our test case, we assume the top wall to be an infinitely long lid moving at a constant velocity, with the rest of the 3 walls being fixed with 'no-slip' condition. The boundary conditions particular to each wall are defined in `Boundary.cpp` under the  `apply` method based on the position of fluid relative to that particular wall in `Case.cpp`.

## Timestep Calculation
In `Fields.cpp`, for adaptive timestep control, `calculate_dt` function was defined which calculates the timestep for next iteration, `dt`, which satisfies the **Courant-Friedrichs-Levi's** (CFL) conditions for `u`, `v`and `nu`.

## Discretization
The convection, diffusion as well as the pressure laplacian terms are discretized according to the finite difference formulation. These are implemented in the `Discretization.cpp`. The convection terms for `u` and `v` is calculated in unique functions, whereas the diffusion function is common to both.  

## Calculation of fluxes and velocity 
The fluxes `F`, `G` are calculated in the `Fields.cpp` using the Discretised form of convection and diffusion terms. The velocities are updated using the `calculate_velocities` function. Also the right side of the Pressure Poisson Equation is being calculated in `Fields.cpp` using the `calculate_rs` function.

## Calculation of pressure

For calculation of Pressure, **Succesive-Over-Relaxation** iterative solver is implemented in `PressureSolver.cpp`, with omega, `omg` as 1.7.

## Plotting Residuals
The functionality of pressure residuals plotting was added to enable the user to monitor the health of the simulation on the fly. To plot the residuals alongside the running simulation, 

1. Install gnuplot and copy the `Residuals.txt` from the `root` folder to the `build` folder. 
```shell
sudo apt-get update -y
sudo apt-get install -y gnuplot-qt
cp ../Residuals.txt .
```

2. After starting the simulation, open terminal and run the following from build.
```shell
gnuplot Residuals.txt
```

3. This would start plotting of the Residuals along side the simulation.  

## Building the code

```shell
git clone https://gitlab.lrz.de/oguzziya/GroupX_CFDLab.git
cd GroupX_CFDLab
mkdir build && cd build
cmake ..
make
```

After `make` completes successfully, you will get an executable `fluidchen` in your `build` directory. Run the code from this directory.
Note: Earlier versions of this documentation pointed to the option of `make install`. You may want to avoid this and directly work inside the repository while you develop the code.

### Build options

By default, **fluidchen** is installed in `DEBUG` mode. To obtain full performance, you can execute cmake as

```shell
cmake -DCMAKE_BUILD_TYPE=RELEASE ..
```

or

```shell
cmake -DCMAKE_CXX_FLAGS="-O3" ..
```

You can see and modify all CMake options with, e.g., `ccmake .` inside `build/` (Ubuntu package `cmake-curses-gui`).

A good idea would be that you setup your computers as runners for [GitLab CI](https://docs.gitlab.com/ee/ci/)
(see the file `.gitlab-ci.yml` here) to check the code building automatically every time you push.

## Running

In order to run **Fluidchen**, the case file should be given as input parameter. Some default case files are located in the `example_cases` directory. Navigate to the `build/` directory and run:

```shell
./fluidchen ../example_cases/LidDrivenCavity/LidDrivenCavity.dat
```

This will run the case file and create the output folder `../example_cases/LidDrivenCavity/LidDrivenCavity_Output`, which holds the `.vtk` files of the solution. 

If the input file does not contain a geometry file (added later in the course), fluidchen will run the lid-driven cavity case with the given parameters.

## Output

In the terminal window, we also output, `Timestep`, `Time`,`Residual`, and `Pressure Poisson Interpretation`. Until convergence, we also output error message. 

## 
### No rule to make target '/usr/lib/x86_64-linux-gnu/libdl.so'

We are investigating an [issue](https://gitlab.lrz.de/tum-i05/public/fluidchen-skeleton/-/issues/3) that appears on specific systems and combinations of dependencies.

## Results and Interpretation

After successfully running the simulation and generating the output files as a **VTK** output files we visualize the result in **Paraview** as follows,


Velocity,<br>
![Velocity Field](/docs/images/Vel.png)

Pressure,<br>
![Pressure](/docs/images/Pre.png)

Tubes with Glyph Cones, representing the direction of fluid flow<br>
![Tubes with Glyph](/docs/images/tubes.png)

### Residuals

For proper simulation paramters, expectedly, we observe convergence, but we also note that the peak value, i.e., initial residual, goes on reducing, indicating a good simulation
![Residuals](/docs/images/Residualconv.png)

For improper simulation paramters, expectedly, we observe divergence, the value of initial residual goes on increasing, until it can no longer be plotted. 
![Residuals](/docs/images/residualdiv.png)

## Problems

### 5. For the SOR solver, for different values of omega, `omg` the behaviour is,

#### `omg` = 0
The `residual` goes on increasing, i.e., the solution does not converge, despite the value of `itermax`.

#### 0 < `omg` < 1
The `residual` goes on decreasing, i.e., the solution does converge, however requires a lot of timesteps. For values of `omg`, closer to 0, convergence is slow and the residual fluctuates.

#### 1 < `omg` < 2
The `residual` goes on decreasing , i.e., the solution does converge, even for normal values of `itermax`.

### 6. Algorithm behaviour for varying dt 

Based on our tests, we conclude that the maximum dt required for stable simulation is **0.009**, while other parameters remained default.

### 7. For varying values of grid size at dt = 0.05

**For `imax` = `jmax` = 16**, we observed convergence for given value of `dt` and `nu = 0.01` the value of `u` at prescribed position is `0.18383 m/s`.

**For `imax` = `jmax` >= 32**, we observed that the solution does not converge for given `dt` and `nu=0.01`, however, upon reducing the value of `dt`, we observe that the solution converges.

### 8. For reducing values of kinematic viscosity `nu` to ( 0.01, 0.002, 0.0005, 0.0001).

We observe that, additional vortices are seen at **TOP LEFT** and **BOTTOM RIGHT**. 

**For `nu` = 0.01**
![nu1](/docs/images/nu1.png)

**For `nu` = 0.002**
![nu2](/docs/images/nu2.png)

**For `nu` = 0.0005**
![nu3](/docs/images/nu3.png)

**For `nu` = 0.0001**
![nu4](/docs/images/nu4.png)







