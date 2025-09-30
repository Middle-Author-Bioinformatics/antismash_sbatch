#!/bin/bash
#SBATCH --job-name=antismash_chunk_aa
#SBATCH --time=6-23:59:00                 # Time limit (D-HH:MM:SS)
#SBATCH --nodes=1
#SBATCH --ntasks=1                       # Number of tasks (1 task using multiple threads)
#SBATCH --cpus-per-task=32                # Number of CPU cores to use
#SBATCH --mem=16G                  # plenty vs. your ~2 GB usage, with headroom
#SBATCH --output=/work/dsamuelson/mmandolfo/job.%J.out
#SBATCH --error=/work/dsamuelson/mmandolfo/job.%J.err
#SBATCH --mail-user=mmandolfo@unmc.edu
#SBATCH --mail-type=END,FAIL

hostname

module load anaconda/25.3

conda activate antismash_meme

ANTISMASH_DATA_PATH=/work/HCC/BCRF/app_specific/antismash/8.0.2/

# DIRECTORY CONTAINING ASSEMBELED METAGENOME CONTIGS
WORKING_DIR=/work/dsamuelson/mmandolfo/antismash_sbatch/
CONTIGS_DIR=/work/dsamuelson/mmandolfo/antismash_sbatch/contigs/
ANTISMASH_PARSED_DIR=/work/dsamuelson/mmandolfo/antismash_sbatch/antismash_parsed/
ANTISMASH_RAW_DIR=/work/dsamuelson/mmandolfo/antismash_sbatch/antismash_raw/
FEGENIEOUTDIR=/work/dsamuelson/mmandolfo/antismash_sbatch/fegenie_allresults/
BLASTOUTDIR=/work/dsamuelson/mmandolfo/antismash_sbatch/blast_outputs/
FINAL=/work/dsamuelson/mmandolfo/antismash_sbatch/final_outputs/

# update -bin_ext below if contig file extensions are different
echo "Running FeGenie on all contig files..."
FeGenie.py -bin_ext fna -bin_dir ${CONTIGS_DIR} -out ${FEGENIEOUTDIR} -cpu 32 --all_results --meta

cd ${CONTIGS_DIR} || exit

for i in *fna; do
    antismash --databases ${ANTISMASH_DATA_PATH} --genefinding-tool prodigal-m --output-dir antismash_${i%.*} --output-basename ${i%.*} --taxon bacteria -c "$SLURM_CPUS_PER_TASK" --cb-knownclusters --cc-mibig --tfbs --fimo --no-enable-genefunctions ${i%.*}
    python antismash2csv.py -o ${ANTISMASH_PARSED_DIR}/${i%.*}.csv ${ANTISMASH_RAW_DIR}/${i%.*}
    mv ${ANTISMASH_PARSED_DIR}/${i%.*}.mibig_hits.csv ${ANTISMASH_PARSED_DIR}/mibig/
    mv ${ANTISMASH_PARSED_DIR}/${i%.*}.tfbs.csv ${ANTISMASH_PARSED_DIR}/tfbs/

    conda activate gtotree
    gtt-genbank-to-AA-seqs -o ${ANTISMASH_PARSED_DIR}/${i%.*}/${i%.*}.faa -i ${ANTISMASH_PARSED_DIR}/${i%.*}/${i%.*}.gbk
    conda deactivate

    conda activate diamond
    diamond makedb --in ${ANTISMASH_PARSED_DIR}/${i%.*}/${i%.*}.faa --db ${ANTISMASH_PARSED_DIR}/${i%.*}/${i%.*}.dmnd
    diamond blastp --db ${i%.*}/${i%.*}.dmnd --query ${FEGENIEOUTDIR}/ORF_calls/$i-proteins.faa --outfmt 6 qseqid sseqid pident length qstart qend sstart send evalue bitscore --out ${BLASTOUTDIR}/${i%.*}.blast --max-target-seqs 1 --threads 16 --id 95
    rm ${i%.*}/${i%.*}.dmnd
    conda deactivate
done

cd ${WORKING_DIR} || exit

# OUTPUT SHOULD INCLUDE THREE FILES: antiSMASH-only.csv, antiSMASH-FeGenie.csv antiSMASH-FeGenie-rescued
python triage.py -a ${ANTISMASH_PARSED_DIR} -f ${FEGENIEOUTDIR} -b ${BLASTOUTDIR} -m ${ANTISMASH_PARSED_DIR}/mibig -o ${FINAL












}
