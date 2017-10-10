#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=500gb,walltime=120:00:00
#PBS -M wisnoski@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/ul-seedbank
module load gcc/6.3.0
module load boost/1.64.0
module load mothur/1.39.5
mothur ul.batch
