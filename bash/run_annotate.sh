#!/bin/bash -l

#SBATCH --job-name=gcompare
#SBATCH --account=y95
#SBATCH --time=24:00:00
#SBATCH --nodes=12
#SBATCH --export=ALL
#SBATCH --output=annotate.out
#SBATCH --error=annotate.err

module swap PrgEnv-cray PrgEnv-gnu
module load java/8u66
module load bwa
module load python/2.7.10
module load samtools
module load mpibash
module load mummer
module load bedtools
export PATH=/home/mderbyshire/localperl/bin:$PATH
export PERL5LIB=/home/mderbyshire/localperl/lib

aprun -n 12 -N 1 -d 24 \
./annotate.sh 24
