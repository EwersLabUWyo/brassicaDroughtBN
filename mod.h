#!/bin/bash
#
#SBATCH --account=plantanalytics
#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mstrong3@uwyo.edu
#SBATCH --output=mod.out

module load intel
module load R/3.3.1

R < BHCdeClusteringExplorationFULL.R --no-save 
