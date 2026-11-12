#!/bin/bash
#SBATCH --job-name=DREAM_transport_sim   # 作业名称，便于在队列中识别
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
set -e

cd /data/zhzhou/DREAM/examples/transport

python generate.py
