#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=16,vmem=400gb,walltime=72:00:00
#PBS -M mmuscare@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/Peralta_IL_Wetlands
module load gcc
module load mothur/1.32.1
mothur IL_mothur_B.batch
