#!/bin/bash
#SBATCH --job-name=antismash_chunk_aa
#SBATCH --time=6-23:59:00                 # Time limit (D-HH:MM:SS)
#SBATCH --nodes=1
#SBATCH --ntasks=1                        # Number of tasks (1 task using multiple threads)
#SBATCH --cpus-per-task=16                # Number of CPU cores to use
#SBATCH --mem=16G                  # plenty vs. your ~2 GB usage, with headroom
#SBATCH --output=/work/dsamuelson/mmandolfo/job.%J.out
#SBATCH --error=/work/dsamuelson/mmandolfo/job.%J.err
#SBATCH --mail-user=mmandolfo@unmc.edu
#SBATCH --mail-type=END,FAIL
CHUNK="chunk_aa"
RUNNAME="chunk_aa"

module load anaconda/25.3

conda activate antismash_meme

ANTISMASH_DATA_PATH=/work/HCC/BCRF/app_specific/antismash/8.0.2/

hostname

while read i; do
    wget $i
    contigs=$(basename $i)
    gzip -d ${contigs}
    contigs=${contigs%.gz}
    antismash --databases $ANTISMASH_DATA_PATH --genefinding-tool prodigal-m --output-dir antismash_${contigs%_contigs*} --output-basename ${contigs%_contigs*} --taxon bacteria -c "$SLURM_CPUS_PER_TASK" --cb-knownclusters --cc-mibig --tfbs --fimo --no-enable-genefunctions ${contigs}
    python /work/dsamuelson/mmandolfo/sbatch_code2/antismash2csv.py -o ${contigs%_contigs*}.csv antismash_${contigs%_contigs*}
    rm -rf ${contigs}
    rm -rf ${contigs}.gz
done < "$CHUNK"


