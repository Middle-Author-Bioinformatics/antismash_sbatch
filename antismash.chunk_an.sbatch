#!/bin/bash
#SBATCH --job-name=antismash_an           # Job name
#SBATCH --time=6-00:00:00                 # Time limit (D-HH:MM:SS)
#SBATCH --ntasks=1                        # Number of tasks (1 task using multiple threads)
#SBATCH --cpus-per-task=40                # Number of CPU cores to use
#SBATCH --mem-per-cpu=1024                # Memory per CPU (1 GB per core = 80 GB total)
#SBATCH --output=/work/dsamuelson/mmandolfo/job.%J.out
#SBATCH --error=/work/dsamuelson/mmandolfo/job.%J.err
#SBATCH --mail-user=mmandolfo@unmc.edu
#SBATCH --mail-type=END,FAIL

module load antismash/8.1

hostname

while read i; do
    wget $i
    contigs=$(basename $i)
    gzip -d ${contigs}
    contigs=${contigs%.gz}
    antismash --genefinding-tool prodigal-m --output-dir antismash_${contigs%_contigs*} --output-basename ${contigs%_contigs*} --taxon bacteria -c 40 --cb-knownclusters --cc-mibig --tfbs --fimo --no-enable-genefunctions ${contigs}
    python /work/dsamuelson/mmandolfo/sbatch_code2/antismash2csv.py -o ${contigs%_contigs*}.csv antismash_${contigs%_contigs*}
    rm ${contigs}
done < chunk_an

