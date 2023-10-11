# QS-MIPS-Numerics
> Numerically solves the governing equations in https://www.biorxiv.org/content/10.1101/2023.04.01.535124v2.

This code implements a Galerkin finite element method with quadratic Lagrange elements and BDF2 time stepping. The implementation uses the open-source library oomph-lib (see https://oomph-lib.github.io/oomph-lib/doc/html/index.html). 

![](header.png)

## Installation

After installing oomph-lib, place the folders 1D, 2D, 1Dsteady in the user_drivers folder and rerun the autogen.sh script in the top directory of the oomph-lib installation. The code in each folder can then be compiled with a call to make. 

## Usage examples

The folders 1D, 2D, and 1D steady each contain independent source files for different versions of the problem. The scripts run.sh in the 1D and 2D folders compute solutions of the time-dependent problem in 1 and 2 spatial dimensions and write snapshots to file. The run.sh script in the 1Dsteady folder computes a solution branch of the steady-state problem via pseudo arc-length continuation  as the parameter D'(u_*) is varied. It outputs solutions at points on the branch and records points on the bifurcation diagram (see preprint above). Once the code is compiled, the run.sh script can be executed to solve the equations. This script specifies all problem parameters, numerical parameters, and output/input file paths. The script can be modified to solve the problems for different parameter sets and folder structure by adjusting the parameter values within the script. Please refer to the bioRxiv preprint above for parameter descriptions. Differences in notation and numerical parameters are noted within the scripts.

The run.sh script writes all output to the result directory 'RESULT0/'. The on-screen output from each processor is written to 'output.XXX' where XXX is the processor number. The numerical solution for n(u,x,t) is written to 'soln_densityYYY_on_procXXX.dat' where YYY indexes the time for the time-dependent problems. The corresponding value of the time variable t is written into the YYY-th row of 't_out.dat'. For steady-state problem, YYY simply indexes the outputs along the solution branch. The numerical solution for c(x,t) and ubar(x,t) are written to soln_helmYYY_on_procXXX.dat. The outputs from each processor for the YYY-th time level can be combined by concatenating the appropriate files, e.g.

cat RESULT0/soln_densityYYY_on_proc* > soln_densityYYY.dat

The format of the soln_helm* and soln_density* files is as follows. From first column to final column, soln_density* contains (u,x,y,n), while soln_helm* contains (u,x,y,c,ubar). The column containing the y-coordinate is omitted in the 1D problems. The steady-state solver generates an additional output file 'bifurcation.dat' which contains the values (D'(u_*), rho(0), and rho(Lx)) along the computed solution branch. These values can be used to construct a bifurcation diagram as in figure 2 of the preprint. Please note: the maximum value of the arc-length step can only be modified from within the source code steady_solve.cc. It is specified by the double variable ds_max. It is recommended to adjust this variable if more or less resolution is desired in computing the bifurcation diagram. Decreasing ds_max can also resolve problems where the continuation 'jumps' between branches. Note also that the solver is not suitable for computing the trivial solution branch and will terminate with an error if this branch is found at any time during continuation. 

Each version of the code also periodically writes restart files. These files contain all of the information required to restart the program from the point where the restart file was generated (see oomph-lib documentation). An example restart script 'restart.sh' is provided for the 1D time-dependent code. For the steady-state problem, a restart can be achieved by uncommenting the marked line in the run.sh script. The frequency at which these files are generated can be adjusted in the unsteady_density.cc and steady_solve.cc source code via the unsigned variable steps_btw_restarts. Restart files are not provided here, but are generated during program execution. 


## CHANGELOG


## Meta

Distributed under the GNU LESSER GENERAL PUBLIC LICENSE. See ``LICENSE`` for more information.
