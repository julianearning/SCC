# Description  
  
Some exercises and results for the module Scientific Computing at the 
University of Applied Sciences Niederrhein. 

- Linear equation solvers: gauss_elimination, qr_decompositon, cholesky_decomposition, jacobi, gauss_seidel, conjugate_gradient  
- Differential equation solvers with the differential equation hardcoded (hunter-prey system): explicit_euler, heun, average, 4 stage runge_kutta  
- Several dimensional differential equations: Jacobi and Gauss-Seidel  
- Partial differential equations for warmth distribution across a stick (warmth.cpp) and oscillation. You can visualize the results of the oscillation using animate.r in R  


# Installation  
  
Uses Eigen3, so to compile use:

g++ -std=c++17 -I <path_to_eigen> <Input Cpp File> -o <Output name>

for example:

g++ -std=c++17 -I /usr/include/eigen3 linear_equations.cpp -o linear_equations

