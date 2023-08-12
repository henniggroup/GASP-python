#!/bin/bash
#SBATCH -J GASPtest # Job name
#SBATCH -n 100 # Number of total cores
#SBATCH -N 1 # Number of nodes
#SBATCH --time=00-00:40:00
#SBATCH -A cts180021p
#SBATCH -p RM-shared
##SBATCH --nodelist=c024,c025
##SBATCH --mem=450GB # Memory pool for all cores in MB
#SBATCH -e outputs_slurm/err #change the name of the err file 
#SBATCH -o outputs_slurm/%x_%A.out # File to which STDOUT will be written %j is the job #
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=kianpu20@cmu.edu # Email to which notifications will be sent

echo "Job started on `hostname` at `date`" 

source ~/.bashrc
source ~/exe/qe.sh
source ~/exe/conda_init.sh

conda activate gasp-qe

# pwifile=$(ls tmp/*.pwi)
# pwiarr=($pwifile)

# echo ${pwiarr[$SLURM_ARRAY_TASK_ID]}
# python run.py ga_input.yaml
#mpirun -np 64 $qe7_pw < $pwiarr > ${pwiarr}.pwo

# echo " "
# echo "Job Ended at `date`"