#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=250gb,walltime=24:00:00
#PBS -M mmuscare@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/Peralta_IL_Wetlands
module load gcc/4.9.2
module load mothur/1.36.1
mothur IL_mothur_A.batch
qsub IL_BacBatch_B.sh
