#!/usr/bin/env mpibash

enable -f mpibash.so mpi_init

threads=$1

mpi_init
mpi_comm_rank rank
mpi_comm_size size

echo "Script three, mapping and variants, began at `date +"%c"` for process $rank of $size."

readsf=($(ls /scratch/y95/mderbyshire/pan_genome/data/P_*R1.fastq))
readsr=($(ls /scratch/y95/mderbyshire/pan_genome/data/P_*R2.fastq))
readsu=($(ls /scratch/y95/mderbyshire/pan_genome/data/U_*.u.fastq))
reference="/scratch/y95/mderbyshire/pan_genome/data/Ssclerotiorum_v2_sorted.fasta"

for ((i=0; i < ${#readsf[@]}; i++)); do
  stub=${readsf[$i]}
  stub=${stub##*/}
  stub=${stub%.*}.$$
  newarray+=(${readsf[$i]},${readsr[$i]},${readsu[$i]})
  stubarray+=($stub)
done

#function check_assembly_finished(){
#  count=1
#  if ls *ass/*contigs.fasta 1> /dev/null 2>&1; then
#    denovoarray=($(ls *ass/*contigs.fasta))
#  else
#    denovoarray=0
#  fi
#  while [ ${#denovoarray[@]} -lt ${#readsf[@]} ]; do
#    if [ $count == 1 ]; then
#      echo ""
#      echo "Wating for assembly to finish."
#      echo ""
#    fi
#    ((count+=1))
#    sleep 0.5
#  done
#}

if [[ $rank == 0 ]]; then
  echo "Running mapping and variant calling for:"
  echo "${stubarray[@]}"
  echo "Running ${#stubarray[@]} on $size * $threads cpus."
fi

./mapping.pl \
--reads ${newarray[$rank]} \
--reference $reference \
--prefix ${stubarray[$rank]}.map \
--threads $threads

if [[ $? != 0 ]]; then
  echo "Problem with mapping for ${stubarray[$rank]}."
else
  echo "Mapping complete for ${stubarray[$rank]}."
fi

mpi_barrier

bamarray=($(ls /scratch/y95/mderbyshire/pan_genome/analysis/P_*/*rg.bam))

./assembly.pl \
--bam ${bamarray[$rank]} \
--prefix ${stubarray[$rank]}.add \
--mitochondrial /scratch/y95/mderbyshire/pan_genome/data/Sclerotinia_sclerotiorum_mtDNA.fasta \
--threads $threads

mpi_barrier

#check_assembly_finished

denovoarray=($(ls /scratch/y95/mderbyshire/pan_genome/analysis/*ass/*contigs.fasta))
bamarray=($(ls /scratch/y95/mderbyshire/pan_genome/analysis/P_*/*rg.bam))

./callVariants.pl \
--reference $reference \
--bam ${bamarray[$rank]} \
--nucmer --denovo ${denovoarray[$rank]} \
--filter-snps --depth 0 --prefix ${stubarray[$rank]}.var --combine \
--threads $threads

if [[ $? != 0 ]]; then
  echo "Problem with calling variants for ${stubarray[$rank]}."
else
  echo "Variant calls complete for ${stubarray[$rank]}."
fi

mpi_barrier
mpi_finalize

echo "Script three, mapping and variants, finished at `date +"%c"`"
