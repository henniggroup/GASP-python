#!/bin/bash
#SBATCH -J GASPrelax # Job name
#SBATCH -n 100 # Number of total cores
#SBATCH -N 1 # Number of nodes
#SBATCH --time=00-10:00:00
#SBATCH -A venkvis_gpu
#SBATCH -p idle
##SBATCH --nodelist=c024,c025
##SBATCH --mem=450GB # Memory pool for all cores in MB
#SBATCH -e outputs_slurm/err #change the name of the err file 
#SBATCH -o outputs_slurm/%x_%A.out # File to which STDOUT will be written %j is the job #
##SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
##SBATCH --mail-user=kianpu20@cmu.edu # Email to which notifications will be sent

echo "Job started on `hostname` at `date`" 

source ~/.bashrc
# source ~/exe/qe.sh
# source ~/exe/conda_init.sh
 
# conda activate gasp-qe

spack load quantum-espresso@7.0

pwifile=$(ls *_relax.pwi)
pwiarr=($pwifile)

srun --mpi=pmix pw.x < ${pwiarr[0]} > ${pwiarr[0]}.pwo

echo " "
echo "Job Ended at `date`"
