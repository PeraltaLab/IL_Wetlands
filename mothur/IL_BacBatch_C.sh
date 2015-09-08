#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=400gb,walltime=12:00:00
#PBS -M mmuscare@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/Peralta_IL_Wetlands
module load gcc
PATH=$PATH:/N/dc2/projects/Lennon_Sequences/mothur/mothur-1.35.1/source
mothur IL_mothur_B.batch
