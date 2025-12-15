#!/bin/bash
#SBATCH -J aim1_baypass
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem 8G
#SBATCH -t 20:00:00
#SBATCH -p standard
#SBATCH -A berglandlab
#SBATCH --array=1-50
#SBATCH -o /scratch/cqh6wn/Class/baypass_project_updated/results/baypass_%a.out
#SBATCH -e /scratch/cqh6wn/Class/baypass_project_updated/results/baypass_%a.err


a=1
# Control file
OPTS=$(sed -n "${SLURM_ARRAY_TASK_ID}"p /scratch/cqh6wn/Class/baypass_project_updated/baypass_control_file.txt)
for i in $OPTS; do declare "opt$a=$i"; ((a++)); done
echo  $opt1
echo  $opt2
echo  $opt3 #  logging

module load gcc/11.4.0


baypass="/scratch/cqh6wn/baypass_public/sources/g_baypass"
cd /scratch/cqh6wn/Class/baypass_project_updated/results

#Running Baypass
$baypass -gfile /scratch/cqh6wn/Class/baypass_project_updated/inputs/subpool_"${opt1}".genobaypass \
-poolsizefile   /scratch/cqh6wn/Class/baypass_project_updated/inputs/subpool_"${opt1}".poolsize \
-outprefix "${opt3}" \
-nthreads 10  \
-contrastfile "${opt2}" \
-seed 12345

echo "Job Complete"
