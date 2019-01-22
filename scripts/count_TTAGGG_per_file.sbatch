#!/bin/bash
#SBATCH --job-name=Rep_Count
#SBATCH --ntasks=30
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1

mkdir repeat_counts

for file in ~/TALES_RNAseq/C101HW18062594/raw_data/*fq.gz
do
   filename=$(basename "${file%.*}")
   srun -n1 -N1 --exclusive bash -c "zgrep -o -n 'TTAGGG' $file | cut -d : -f 1 | uniq -c | sed -e 's/  \+//' > repeat_counts/${filename}_TERRA.repeat_count" &
   srun -n1 -N1 --exclusive bash -c "zgrep -o -n 'CCCTAA' $file | cut -d : -f 1 | uniq -c | sed -e 's/  \+//' > repeat_counts/${filename}_ARRET.repeat_count" &
done
wait

echo "Statistics for job $SLURM_JOB_ID:"
sacct --format="JOBID,NodeList,NNodes,Start,End,Elapsed,AllocCPUs,CPUTime,AveDiskRead,AveDiskWrite,MaxRSS,MaxVMSize,exitcode,derivedexitcode" -j $SLURM_JOB_ID

