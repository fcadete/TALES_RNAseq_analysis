#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --time=4:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --image=docker:ummidock/trimmomatic:0.36
#SBATCH --workdir=/mnt/beegfs/scratch/SALMEIDA/fcadete/TALES_RNAseq_analysis

mkdir filtered_fastq

for f in ~/TALES_RNAseq/C101HW18062594/raw_data/*_1.fq.gz; do

filename=$(basename "${f%.*}").gz


srun --exclusive -n1 -N1 shifter bash -c "java -jar /NGStools/Trimmomatic-0.36/trimmomatic.jar PE -threads 16 \
      $f ${f//_1.fq/_2.fq} \
      filtered_fastq/${filename//_1.fq/_1_filtered_paired.fq} filtered_fastq/${filename//_1.fq/_1_filtered_unpaired.fq} \
      filtered_fastq/${filename//_1.fq/_2_filtered_paired.fq} filtered_fastq/${filename//_1.fq/_2_filtered_unpaired.fq} \
      ILLUMINACLIP:"/NGStools/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa":2:30:10 \
      ILLUMINACLIP:"/NGStools/Trimmomatic-0.36/adapters/TruSeq3-PE.fa":2:30:10 \
      SLIDINGWINDOW:5:20 MINLEN:50 &> filtered_fastq/${filename//_1.fq.gz/_trimmomatic_output.txt}" &

done

wait

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

