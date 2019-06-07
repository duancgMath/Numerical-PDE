# 1D heat conduction equation solver

`heat_conduct_1d.m` define a problem class: Intial-Boundary-Condition Problem, Dirichlet conditions

`solver.m` equation solver

`update_solution.m` update solution at each time step: forward, backward, Crank-Nickson

`main_parallel.m`  parallel computing; accuracy order result in **space**

`main_serial.m`  serial computing; accuracy order result in **space**

#### Error and Accuracy Order

`error_analysis.m` calculate error

`print_error.m` print and plot error and accuracy order result
