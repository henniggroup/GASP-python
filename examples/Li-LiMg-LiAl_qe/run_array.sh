#!/bin/bash
#SBATCH -J GASPeos # Job name
#SBATCH -n 100 # Number of total cores
#SBATCH -N 1 # Number of nodes
#SBATCH --time=00-02:00:00
#SBATCH -A cts180021p
#SBATCH -p RM
##SBATCH --nodelist=c024,c025
##SBATCH --mem=450GB # Memory pool for all cores in MB
#SBATCH -e outputs_slurm/%A_%a.err #change the name of the err file 
#SBATCH -o outputs_slurm/%A_%a.out # File to which STDOUT will be written %j is the job #
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=kianpu20@cmu.edu # Email to which notifications will be sent

echo "Job started on `hostname` at `date`" 

source ~/.bashrc
source ~/exe/qe.sh
# source ~/exe/conda_init.sh
 
# conda activate gasp-qe

pwifile=$(ls *.pwi)
pwiarr=($pwifile)

mpirun -np 100 $qe7_pw < ${pwiarr[$SLURM_ARRAY_TASK_ID]} > ${pwiarr[$SLURM_ARRAY_TASK_ID]}.pwo

echo " "
echo "Job Ended at `date`"