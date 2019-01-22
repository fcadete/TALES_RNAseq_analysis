#!/bin/bash
#SBATCH --job-name=rRNA-sortmerna
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=32
#SBATCH --workdir=/mnt/beegfs/scratch/SALMEIDA/fcadete/TALES_RNAseq_analysis

mkdir sortmerna_to_rRNA

SORTMERNAFOLDER=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/sortmerna-2.1b
RRNADATABASES=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/sortmerna-2.1b/rRNA_databases
RRNAINDEXES=/mnt/nfs/lobo/SALMEIDA-NFS/fcadete/TRF2_siRNA/sortmerna-2.1b/index

#srun -n1 -N1 --exclusive ./sortmerna-2.1b/indexdb_rna --ref $RRNADATABASES/silva-bac-16s-id90.fasta,$RRNAINDEXES/silva-bac-16s-db:\
#$RRNADATABASES/silva-bac-23s-id98.fasta,$RRNAINDEXES/silva-bac-23s-db:\
#$RRNADATABASES/silva-arc-16s-id95.fasta,$RRNAINDEXES/silva-arc-16s-db:\
#$RRNADATABASES/silva-arc-23s-id98.fasta,$RRNAINDEXES/silva-arc-23s-db:\
#$RRNADATABASES/silva-euk-18s-id95.fasta,$RRNAINDEXES/silva-euk-18s-db:\
#$RRNADATABASES/silva-euk-28s-id98.fasta,$RRNAINDEXES/silva-euk-28s:\
#$RRNADATABASES/rfam-5s-database-id98.fasta,$RRNAINDEXES/rfam-5s-db:\
#$RRNADATABASES/rfam-5.8s-database-id98.fasta,$RRNAINDEXES/rfam-5.8s-db

for f in filtered_fastq/*1_filtered_paired.fq.gz; do

filename=$(basename "${f%.*}")
filename=${filename//.fq/}

filetocheck="sortmerna_to_rRNA/${filename//_1_filtered/_non_rRNA_1_filtered}.fq"

if [ -s "$filetocheck" ]; then  continue; fi

srun -n1 -N1 --exclusive bash -c "gunzip -c $f > sortmerna_to_rRNA/${filename}.fq
                                  gunzip -c ${f//_1_filtered/_2_filtered} > sortmerna_to_rRNA/${filename//_1_filtered/_2_filtered}.fq
                                  $SORTMERNAFOLDER/scripts/merge-paired-reads.sh sortmerna_to_rRNA/${filename}.fq sortmerna_to_rRNA/${filename//_1_filtered/_2_filtered}.fq sortmerna_to_rRNA/${filename//_1_filtered/_merged_filtered}
                                  rm sortmerna_to_rRNA/${filename}.fq
                                  rm sortmerna_to_rRNA/${filename//_1_filtered/_2_filtered}.fq
                                 $SORTMERNAFOLDER/sortmerna --ref $RRNADATABASES/silva-bac-16s-id90.fasta,$RRNAINDEXES/silva-bac-16s-db:\
$RRNADATABASES/silva-bac-23s-id98.fasta,$RRNAINDEXES/silva-bac-23s-db:\
$RRNADATABASES/silva-arc-16s-id95.fasta,$RRNAINDEXES/silva-arc-16s-db:\
$RRNADATABASES/silva-arc-23s-id98.fasta,$RRNAINDEXES/silva-arc-23s-db:\
$RRNADATABASES/silva-euk-18s-id95.fasta,$RRNAINDEXES/silva-euk-18s-db:\
$RRNADATABASES/silva-euk-28s-id98.fasta,$RRNAINDEXES/silva-euk-28s:\
$RRNADATABASES/rfam-5s-database-id98.fasta,$RRNAINDEXES/rfam-5s-db:\
$RRNADATABASES/rfam-5.8s-database-id98.fasta,$RRNAINDEXES/rfam-5.8s-db \
--reads sortmerna_to_rRNA/${filename//_1_filtered/_merged_filtered} \
--aligned sortmerna_to_rRNA/${filename//_1_filtered/_rRNA_merged_filtered} \
--other sortmerna_to_rRNA/${filename//_1_filtered/_non_rRNA_merged_filtered} \
--fastx \
--log \
--paired_in \
-a 32
                                  $SORTMERNAFOLDER/scripts/unmerge-paired-reads.sh sortmerna_to_rRNA/${filename//_1_filtered/_non_rRNA_merged_filtered}.fastq \
                                                                                    sortmerna_to_rRNA/${filename//_1_filtered/_non_rRNA_1_filtered}.fq \
                                                                                    sortmerna_to_rRNA/${filename//_1_filtered/_non_rRNA_2_filtered}.fq
                                   rm sortmerna_to_rRNA/${filename//_1_filtered/_merged_filtered}
                                   rm sortmerna_to_rRNA/${filename//_1_filtered/_non_rRNA_merged_filtered}.fastq
                                   rm sortmerna_to_rRNA/${filename//_1_filtered/_rRNA_merged_filtered}.fastq" &
done
wait

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID


