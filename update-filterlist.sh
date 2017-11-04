#!/bin/bash

## CALL: ./update-filterlist.sh <path/to/wave-functions/directory> <start digit> <end digit>
## EXAMPLE: ./update-filterlist.sh ../he-wfs 001 123

DENPATH=${1}
START=${2}
FINISH=${3}
SLURM='launch.slurm'

rm -v ${DENPATH}/filter.txt
for ID in  $(seq -w ${START} 1 ${FINISH})
do
	echo "density.${ID}.dat" >> ${DENPATH}/filter.txt
done
echo "New filter list created"
N=$(cat ${DENPATH}/filter.txt | wc -l)
echo "${N} entries in filter list"
if ((N<180))
then
	if ((N%20>0))
	then
		sed -i "/^#SBATCH --nodes=/c\#SBATCH --nodes=$((N/20+1))" ${SLURM}
	else
		sed -i "/^#SBATCH --nodes=/c\#SBATCH --nodes=$((N/20))" ${SLURM}
	fi
	sed -i "/^#SBATCH --ntasks=/c\#SBATCH --ntasks=${N}" ${SLURM}
else
	sed -i "/^#SBATCH --nodes=/c\#SBATCH --nodes=9" ${SLURM}
	sed -i "/^#SBATCH --ntasks=/c\#SBATCH --ntasks=180" ${SLURM}
fi
echo "'${SLURM}' updated."
