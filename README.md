# Computational Fluid Dynamics – course projects in MATLAB

Code from two assignments in the course 41317 Computational Fluid Dynamics at the Technical University of Denmark.

Original assignments: Bogusz Glaza and Giovanni Girardi, spring 2025.

Summary and code cleanup: Bogusz Glaza, November 2026.

A summary of both projects, with figures and the calculations explained, is in [CFD_Project_Summary.pdf](CFD_Project_Summary.pdf).

The solvers are written from scratch, starting from the conservation laws, and every part was checked against exact solutions or a benchmark.

## report1 – heat transport by a prescribed flow

Steady convection–diffusion of temperature in a square cavity with a given velocity field, solved with the finite volume method. The convective term uses central (CDS) or upwind (UDS) differencing, and boundary conditions are imposed with ghost points.

- `FiniteVolumeConvectionDiffusion1D.m` – 1D model problem, compared with the exact solution
- `FiniteVolumeConvectionDiffusion2D.m` – 2D solver for a uniform flow and a corner flow: convergence against the exact solution, global heat conservation and the heat flux through the western wall
- `FiniteVolumeConvectionDiffusion2D_UDS_pe.m` – heat flux through the western wall as a function of the Péclet number

## report2 – lid-driven cavity

Unsteady, incompressible Navier–Stokes equations on a staggered grid, with explicit Euler time stepping and a pressure Poisson equation. A discrete continuity correction keeps mass conserved to machine precision, and the Poisson matrix is LU-factorised once before the time loop.

- `NavierStokes2dMovingLidSquareCavityFlow.m` – main solver
- `NS2dHfunctions.m` – convective and diffusive terms
- `NS2dLaplaceMatrix.m` – Laplace matrix for the pressure equation
- `StaggeredMesh2dSquare.m` – coordinates of the staggered grid (provided with the course)
- `NS2dValidateHfunction.m`, `NS2dValidatePoissonSolver.m` – tests against exact solutions
- `NS2dRe1000ConvergenceAnalysis.m` – comparison with the benchmark of Bruneau & Saad (2006) at Re = 1000 on five grids

## Running

The scripts run in MATLAB without toolboxes. Open the folder, set the parameters at the top of a script and run it. Figures are saved in a `plots` folder. The scripts are built on the course templates, so some comments (`<< ... >>`) are the original template instructions. The convergence study at Re = 1000 takes several minutes on the finest grid.

## References

- H. B. Bingham, P. S. Larsen and V. A. Barker, *Computational Fluid Dynamics*, lecture notes for course 41317, DTU.
- C.-H. Bruneau and M. Saad, The 2D lid-driven cavity problem revisited, *Computers & Fluids* 35 (2006) 326–348.
