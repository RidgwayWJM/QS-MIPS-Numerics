#!/bin/bash
#SBATCH --ntasks=75
#SBATCH --array=0

module load gcc/gcc11

# Executes the unsteady_density code multiple times for a range of parameters

i=$SLURM_ARRAY_TASK_ID

# setup result directory
RESULT_DIR="RESULT${i}/"

#determine the index of the last restart file
NUMFILES=$(ls ${RESULT_DIR}/restar*_on_proc0.dat | wc -l)
RESTART_IND=$(($NUMFILES-1))

echo "Restart index = ${RESTART_IND}" 

RESTART_FILE="${RESULT_DIR}restart${RESTART_IND}"
#RESTART_FILE="HOPF_RESTART/restart147"
PARTITION_FILE="${RESULT_DIR}partition.dat"
#PARTITION_FILE="RESULT0/partition.dat"

#input file name
INPUT_FILE="${pwd}input${i}.dat"

mpirun ./unsteady_density $INPUT_FILE $RESULT_DIR --restart-file $RESTART_FILE --partition-file $PARTITION_FILE
