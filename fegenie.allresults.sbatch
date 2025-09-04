#!/bin/bash
#SBATCH --job-name=antismash_chunk_ao
#SBATCH --time=6-23:59:00                 # Time limit (D-HH:MM:SS)
#SBATCH --nodes=1
#SBATCH --ntasks=1                        # Number of tasks (1 task using multiple threads)
#SBATCH --cpus-per-task=16                # Number of CPU cores to use
#SBATCH --mem=16G                  # plenty vs. your ~2 GB usage, with headroom
#SBATCH --output=/work/dsamuelson/mmandolfo/job.%J.out
#SBATCH --error=/work/dsamuelson/mmandolfo/job.%J.err
#SBATCH --mail-user=mmandolfo@unmc.edu
#SBATCH --mail-type=END,FAIL