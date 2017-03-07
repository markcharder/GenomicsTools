#!/usr/bin/env mpibash

enable -f mpibash.so mpi_init

threads=$1

mpi_init
mpi_comm_rank rank
mpi_comm_size size

echo "Script two, annotate, started at `date +"%c"`"
denovoarray=($(ls /scratch/y95/mderbyshire/pan_genome/analysis/*ass/*contigs.fasta))

readsf=($(ls /scratch/y95/mderbyshire/pan_genome/data/P_*R1.fastq))

for ((i=0; i < ${#readsf[@]}; i++)); do
  stub=${readsf[$i]}
  stub=${stub##*/}
  stub=${stub%.*}.$$
  newarray+=(${readsf[$i]},${readsr[$i]},${readsu[$i]})
  stubarray+=($stub)
done

if [[ $rank == 0 ]]; then
  echo "Running annotation for:"
  echo "${denovoarray[@]}"
  echo "Running ${#denovoarray[@]} on $size * $threads cpus."
fi

./annotation.pl \
--type all \
--contigs ${denovoarray[$rank]} \
--species Sclerotinia_sclerotiorum \
--prot /scratch/y95/mderbyshire/pan_genome/data/Ssclerotiorum_v2.pep.fasta \
--threads $threads \
--prefix ${stubarray[$rank]}.ann

if [[ $? != 0 ]]; then
  echo "Problem with calling annotation for ${stubarray[$rank]}."
else
  echo "Annotation complete for ${stubarray[$rank]}."
fi
mpi_barrier
mpi_finalize

echo "Script two, annotate, finished at `data +"%c"`" 
