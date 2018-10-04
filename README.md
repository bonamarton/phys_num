# phys_num
Numerical Methods for Physics

_______________________________________________

## Informations:

### From the project:

Various numerical methods, written in C language.
The project contains the solutions of an university course.
Tasks description can find [here](http://www.vo.elte.hu/~dobos/teaching/fiznum2018/default.aspx) (Hungarian)

### This project contains:

The project contains four main topic:

1. Matrix multipication (folder: matrix_multi):

The program can multiplicat row, column vectors, matrices and their combinations

2. Gauss-Jordan elimination (folder: gauss_jordan_elim):

Their are two version, first is the standard algorithmic solution,
second version contains singular-value decomposition (SVD) method.
For SVD, LAPACK needed. You can get from [here](http://icl.cs.utk.edu/lapack-for-windows/lapack/).
In simple case (in win ops sys), decompress lapacke.rar (in "gauss_jordan_elim") to "gauss_jordan_elim" folder, and it will work.

3. Polinomial fitting (folder: multi_dim_fitt):

The program make a general N variable polynomial fitting.
One small data for test added.
Plots and result in "predicted_values" folder

4. Numerical methods for ordinary differential equations:

4.1. Euler method (folder: euler_RK4_moon):
Calculates the solution of the Earth-Moon system with Euler method.

4.2. 4th order runge-kutta method (folder: euler_RK4_moon and RK4_attractors):

4.2.1. Calculates the solution of the Earth-Moon system with RK4 method.
(plots in euler_RK4_moon/plot folder)

4.2.2. Strange Attractors
Calculates the solution of special ODE systems with adaptive RK4 method.
(plots of the attractors find in RK4_attractors/plot folder,
dataset in /data for those who want make plots)

### Teaser:


![](../master/)


![](../master/)


![](../master/)
