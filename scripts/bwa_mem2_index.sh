#!/bin/bash
#SBATCH --job-name=bwa_mem2_index
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amanda.katzer@nih.gov
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=70GB
#SBATCH --output=/data/katzeram/hydractinia_qtl/genome/out_%j.log


cd /data/katzeram/hydractinia_qtl/genome/
module load bwa-mem2
bwa-mem2 index -p hydractiniaT2T GCF_029227915.1_HSymV2.1_genomic.fna
