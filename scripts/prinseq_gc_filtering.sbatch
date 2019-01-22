#!/bin/bash
#SBATCH --job-name=prinseq
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --image=docker:flowcraft/prinseq:0.20.4-1
#SBATCH --workdir=/mnt/beegfs/scratch/SALMEIDA/fcadete/TALES_RNAseq_analysis

mkdir prinseq_outputs

for f in filtered_fastq/*_1_filtered_paired.fq.gz; do

 filename=$(basename $f)

 filename2=${filename//1_filtered/2_filtered}

 srun -n1 -N1 --exclusive shifter bash -c "gunzip -c filtered_fastq/$filename > prinseq_outputs/${filename//.gz/}
                                           gunzip -c filtered_fastq/$filename2 > prinseq_outputs/${filename2//.gz}
                                           prinseq-lite.pl -fastq prinseq_outputs/${filename//.gz/} \
                                                           -fastq2 prinseq_outputs/${filename2//.gz} \
                                                           -out_bad null \
                                                           -min_gc 70
                                           rm prinseq_outputs/${filename//.gz/}
                                           rm prinseq_outputs/${filename2//.gz}" &
done
wait

sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID



