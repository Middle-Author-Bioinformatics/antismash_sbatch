#!/bin/bash
#SBATCH --job-name=antismash_chunk_ao
#SBATCH --time=6-23:59:00                 # Time limit (D-HH:MM:SS)
#SBATCH --nodes=1
#SBATCH --ntasks=1                        # Number of tasks (1 task using multiple threads)
#SBATCH --cpus-per-task=32                # Number of CPU cores to use
#SBATCH --mem=32G
#SBATCH --output=/work/dsamuelson/mmandolfo/job.%J.out
#SBATCH --error=/work/dsamuelson/mmandolfo/job.%J.err
#SBATCH --mail-user=mmandolfo@unmc.edu
#SBATCH --mail-type=END,FAIL


module load fegenie/1.2 # UPDATE IF NEEDED

hostname

echo "Downloading contig files..."
while read i; do
  wget $i
done < "all_contig_links.txt"

echo "Creating contig directory..."
mkdir ibdmd_contigs

echo "Moving contig files to directory..."
mv ./*fna.gz ibdmd_contigs

echo "Decompressing contig files..."
gzip -d ibdmd_contigs/*gz

echo "Running FeGenie on all contig files..."
FeGenie.py -bin_ext fna -bin_dir ibdmd_contigs/ -out all_fegenie_results -cpu 32 --all_results --heme --meta