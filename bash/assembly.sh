#!/usr/bin/env mpibash

enable -f mpibash.so mpi_init

threads=$1

mpi_init
mpi_comm_rank rank
mpi_comm_size size

echo "Script one, assembly, began at `date +"%c"` for process $rank of $size."

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

if [[ $rank == 0 ]]; then
  echo "Running assembly for:"
  echo "${stubarray[@]}"
  echo "Running ${#stubarray[@]} on $size * $threads cpus."
fi

./assembly.pl \
--reads ${newarray[$rank]} \
--prefix ${stubarray[$rank]}.ass \
--threads $threads

if [[ $? != 0 ]]; then
  echo "Problem with assembly for ${stubarray[$rank]}."
else
  echo "Assembly complete for ${stubarray[$rank]}."
fi

mpi_barrier
mpi_finalize
echo "Script one, assembly, ended at `date +"%c"`"
