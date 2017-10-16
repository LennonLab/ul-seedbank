#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=500gb,walltime=120:00:00
#PBS -M wisnoski@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/ul-seedbank
module load gcc
module load boost
module load mothur
mothur ul.batch
