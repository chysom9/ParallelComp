#!/bin/bash
#SBATCH --job-name=chi_seq_client
#SBATCH --output=logs/seq_out.txt   
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=Centaurus

./level_client "Tom Hanks" 3








