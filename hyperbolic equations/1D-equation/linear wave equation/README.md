# 1D linear wave equation solver

`linear_wave_problem.m` define a problem class

`update_solution.m` update solution at each time step

#### Initial conditions

`init_cond.m` initial conditions from book << Numerical methods for conservation law >>

`init_cond_riemann.m` initial conditions for Riemann problem

`init_cond_normal.m` initial conditions: normal distribution

#### Numerical flux

`laxfriedrich_flux.m` define Lax-Friedrich flux

`laxwendroff_flux.m` define Lax-Wendroff flux << Finite Volume Methods for Hyperbolic Problems >> Randall LeVeque

#### Error and Accuracy Order

`error_analysis.m` calculate error

`print_error.m` print and plot error and accuracy order result

