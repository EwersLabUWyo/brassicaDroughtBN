#!/bin/bash
#
#SBATCH --account=plantanalytics
#SBATCH --time=168:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --constraint=hugemem
#SBATCH --partition=hugemem
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mstrong3@uwyo.edu
#SBATCH --output=net.out

module load intel
module load R/3.3.1

R < BL_largeRNApheno.R --no-save 
