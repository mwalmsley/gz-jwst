#!/bin/bash
#SBATCH --time=71:00:00                                # Time limit hrs:min:sec
#SBATCH --constraint=A100 
#SBATCH --mem=30G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=16

nvidia-smi

srun /share/nas2/walml/miniconda3/envs/zoobot39_cu118_dev/bin/python /share/nas2/walml/repos/zoobot_jwst/zoobot_jwst.py