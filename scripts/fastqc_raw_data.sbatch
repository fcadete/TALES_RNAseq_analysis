#!/bin/bash
#SBATCH --job-name=fastQC
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.1
#SBATCH --workdir=/mnt/beegfs/scratch/SALMEIDA/fcadete/TALES_RNAseq_analysis

mkdir fastqc_reports

for f in ~/TALES_RNAseq/C101HW18062594/raw_data/*fq.gz; do

filename=$(basename "${f%.*}")

srun --exclusive -n1 -N1 shifter bash -c "zcat $f > fastqc_reports/$filename; fastqc fastqc_reports/$filename; rm fastqc_reports/$filename" &

done

wait


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
