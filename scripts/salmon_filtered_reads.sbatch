#!/bin/bash
#SBATCH --job-name=salmon
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=40
#SBATCH --image=docker:mcfonsecalab/salmon:latest
#SBATCH --workdir=/mnt/beegfs/scratch/SALMEIDA/fcadete/TALES_RNAseq_analysis

mkdir salmon_on_hg38_output

srun -n1 -N1 --exclusive shifter salmon index -t ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all.fa.gz \
                                              -i ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all

for f in filtered_fastq/*1_filtered_paired.fq.gz; do

filename=$(basename "${f%.*}").gz

srun -n1 -N1 --exclusive shifter salmon quant -i ensembl_transcriptome/Homo_sapiens.GRCh38.cdna.all -l A -p 40 --gcBias \
                   -o salmon_on_hg38_output/${filename//_1_filtered_paired.fq.gz/} \
                   -1 filtered_fastq/$filename \
                   -2 filtered_fastq/${filename//_1_filtered/_2_filtered} &
done
wait


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

