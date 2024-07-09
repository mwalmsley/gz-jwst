#!/bin/bash
#SBATCH --time=71:15:00  
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu 4G
#SBATCH --job-name=download-jwst

pwd; hostname; date

nvidia-smi

PYTHON=/home/walml/envs/zoobot39_dev/bin/python

REPO_DIR=/project/def-bovy/walml/repos/gz-jwst
srun $PYTHON $REPO_DIR/TODO
