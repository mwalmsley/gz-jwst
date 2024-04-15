#!/bin/bash
#SBATCH --job-name=alldwn_%x.%A_%a.out
#SBATCH --output=alldwn_%x.%A_%a.out.log
#SBATCH --time=71:00:00                                # Time limit hrs:min:sec
#SBATCH --constraint=A100 
#SBATCH --mem=30G
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=16

srun /share/nas2/walml/miniconda3/envs/zoobot39_cu118_dev/bin/python /share/nas2/walml/repos/zoobot-foundation/foundation/experiments/fin