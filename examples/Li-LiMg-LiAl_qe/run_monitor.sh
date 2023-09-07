#!/bin/bash
#SBATCH -J mainTEST # Job name
#SBATCH -n 4 # Number of total cores
#SBATCH -N 1 # Number of nodes
#SBATCH --time=00-05:00:00
#SBATCH -A cts180021p
#SBATCH -p RM-shared
##SBATCH --nodelist=c024,c025
#SBATCH --mem-per-cpu=2000MB # Memory pool for all cores in MB
#SBATCH -e outputs_slurm/slurm-%j.err #change the name of the err file 
#SBATCH -o ga_output.out # File to which STDOUT will be written %j is the job #
##SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=kianpu20@cmu.edu # Email to which notifications will be sent

echo "Job started on `hostname` at `date`" 

source ~/.bashrc
source ~/exe/qe.sh
source ~/exe/conda_init.sh

conda activate gasp-qe

python -u run_qe.py ga_input.yaml
