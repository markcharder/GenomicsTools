#!/bin/bash -l

#SBATCH --job-name=gcompare
#SBATCH --account=y95
#SBATCH --time=24:00:00
#SBATCH --nodes=12
#SBATCH --export=none
#SBATCH --output=mapping_variants.out
#SBATCH --error=mapping_variants.err

module swap PrgEnv-cray PrgEnv-gnu
module load java/8u66
module load bwa
module load python/2.7.10
module load samtools
module load mpibash
module load mummer
module load bedtools

aprun -n 12 -N 1 -d 24 \
./mapping_variants.sh 24

