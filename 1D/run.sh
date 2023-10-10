#!/bin/bash
#Example usage script for steady-state code. This one assumes that we are running the code via submission to the SLURM scheduler
#SBATCH --ntasks=75
#SBATCH --array=0
#SBATCH --exclude=hpc[1-13]


#load relevant c++ module (modify as required)
module load gcc/gcc11

i=$SLURM_ARRAY_TASK_ID

# setup result directory
RESULT_DIR="RESULT${i}/"
mkdir $RESULT_DIR

#input file name
INPUT_FILE="${pwd}input${i}.dat"

#generate a new input file  
# output problem parameters to input file 
#**************** Problem Parameters *************************
printf "D2 = 0.1\n" > $INPUT_FILE	#D2=Dc in paper
printf "beta = 3.0\n" >> $INPUT_FILE
printf "alpha = 1\n" >> $INPUT_FILE
printf "D1_star = 0.1\n" >> $INPUT_FILE	#D1_star = D(u*) in paper
printf "D1_p_star = 2.0\n" >>  $INPUT_FILE	#D1_p_star = D'(u*) in paper
printf "D_inf_0 = 0.01\n" >> $INPUT_FILE	#D_inf_0 = D_\infty if D(u) decreasing, and D_inf_0 = D_0 if D(u) increasing
printf "Lu = 1.5\n" >>  $INPUT_FILE	#domain size in u-coordinate
printf "Lx = 8.0\n" >>  $INPUT_FILE	#domain size in x-coordinate
	
#**************** kinetic parameters ***********************
printf "a = 0.1\n" >>  $INPUT_FILE
printf "L = 20.0\n" >>  $INPUT_FILE
printf "K = 20.0\n" >>  $INPUT_FILE
printf "lambda = 1.0\n" >>  $INPUT_FILE
printf "rho_star=5.0\n" >>  $INPUT_FILE
printf "epsilon = 1.0e-4\n" >>  $INPUT_FILE
	
#****************** numerical parameters **********************
printf "epsilon_t = 9e-4\n" >>  $INPUT_FILE	#approximate time stepping error tolerance 
printf "t_max = 1e6\n" >>  $INPUT_FILE	#compute up to t =tmax
printf "dt_out = 0.1\n" >>  $INPUT_FILE	#max time interval for recording solution
printf "Nu_density = 600\n" >>  $INPUT_FILE	#number of elements in u-coordinate for n-mesh
printf "Nu_helm = 1\n" >>  $INPUT_FILE	#number of elements in u-coordinate for c-mesh
printf "Nx = 80\n" >>  $INPUT_FILE	#number of elements in x-coordinate for all meshes


mpirun ./unsteady_density $INPUT_FILE $RESULT_DIR 

