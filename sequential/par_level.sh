#!/bin/bash
#SBATCH --job-name=chi_par_client
#SBATCH --output=logs/par_out.txt
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=Centaurus

./par_level_client "Tom Hanks" 3

