#!/bin/bash
#SBATCH --job-name="TCV Te scan"
#SBATCH --output=logs/dream_%j.out
#SBATCH --error=logs/dream_%j.err
#SBATCH --time=40:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


export OMP_NUM_THREADS=4
cd /home/lowebl/DREAM/examples/warepinch
python3 basic.py

