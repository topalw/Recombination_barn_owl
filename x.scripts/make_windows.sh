#!/bin/bash
#SBATCH --partition cpu
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --job-name=windows
#SBATCH --account jgoudet_barn_owl

module load gcc bedtools2 r 

# run like make_windows.sh pyrho_total_file window_file

P=${1}
WIND=${2}

bedtools intersect -b ${WIND} -a ${P} -wo > "${P}.intersect.bed"

Rscript 2.make_windows.r "${P}.intersect.bed" ${WIND}

