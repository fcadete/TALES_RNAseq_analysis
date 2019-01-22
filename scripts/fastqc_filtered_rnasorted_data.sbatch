#!/bin/bash
#SBATCH --job-name=fastQC
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:argrosso/htspreprocessing:0.1.1
#SBATCH --workdir=/mnt/beegfs/scratch/SALMEIDA/fcadete/TALES_RNAseq_analysis

for f in sortmerna_to_rRNA/*non_rRNA*filtered_paired.fq; do

srun --exclusive -n1 -N1 shifter fastqc $f &

done

wait


echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID
