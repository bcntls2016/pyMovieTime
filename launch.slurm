#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --ntasks-per-node=20
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --time 0:10:0

#SBATCH -J "pyMovie"

module purge
module load intel/16.1.3 intelmpi/5.1.3.210

export I_MPI_WAIT_MODE=1 
export OMP_NUM_THREADS=1
export WORKDIR=$(pwd)
cd ${WORKDIR}

date
echo "Execution starting"

mpirun chdb \
	--in-dir ../he-wfs \
	--in-files ../he-wfs/filter.txt \
	--out-dir Movie.${SLURM_JOBID} \
	--in-type dat \
	--command-line './movietime.sh 5 %name% %in-dir% %out-dir%' \
	--out-files %name% \
	--on-error=errors.${SLURM_JOBID} \
	--report=report.${SLURM_JOBID}

date
echo "Execution finished"

infoincidentjob
jobinfo ${SLURM_JOB_ID}
