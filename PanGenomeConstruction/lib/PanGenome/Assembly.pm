package PanGenome::Assembly;
use strict;
use warnings;
use Cwd;
use File::Basename;

sub exitStatus {

  my $message = shift @_;
  my $exit = shift @_;

  if ($exit != 0){
    die ($message);
    }

  }


sub getAdditionalReads {

  my $prefix = shift @_;
  my $bam = shift @_;
  my $count = shift @_;
  my $baseName = basename($bam);
  my $dir = getcwd;

  my $exit = system("bam2fastq -f -o $prefix/$prefix#.fastq --no-aligned --unaligned --no-filter $bam >> $prefix/bam2fq.$$.logfile.txt 2>&1");

  exitStatus("Problem running bam2fastq.\n", $exit);  

  return ("$dir/$prefix/${prefix}_1.fastq","$dir/$prefix/${prefix}_2.fastq","$dir/$prefix/${prefix}_M.fastq");

  }

sub assembleA5 {

  my $dir = getcwd;
  my $prefix = shift @_;
  my @readsFiles = @_;
  my $baseName = basename($readsFiles[0]);

  if ((@readsFiles == 2) or (@readsFiles == 3)) {
    open my $fh, ">>", "$prefix/$prefix.lib";
    print $fh "[LIB]\np1=$readsFiles[0]\np2=$readsFiles[1]\n";

    if (@readsFiles == 3) {
      print $fh "up=$readsFiles[2]\n";
      }

    close $fh;
    }

  if (@readsFiles == 1) {
    open my $fh, ">>", "$prefix/$prefix.lib";
    print $fh "[LIB]\nup=$readsFiles[0]\n";
    }

  chdir "$prefix";

  my $exit = system("a5_pipeline.pl $prefix.lib $prefix.out >> A5.$$.logfile.txt 2>&1");
  
  chdir $dir;

  exitStatus("Problem running A5 pipeline for $prefix process $$.\n", $exit);
  return "$prefix/$prefix.out.contigs.fasta";

  }

sub filterMitochondrial {

  my $prefix = shift @_;
  my $mitochondrial = shift @_;
  my $contigs = shift @_;
  system("nucmer $mitochondrial $contigs --prefix $prefix/mitochondrial >> $prefix/filter_contigs.$$.logfile.txt 2>&1");
  system("show-coords $prefix/mitochondrial.delta > $prefix/mitochondrial.coords 2>> $prefix/filter_contigs.$$.logfile.txt");
  system("awk '{if(\$0~/^=/){s=1};if(s==1){print \$13}}' $prefix/mitochondrial.coords | sort -u > $prefix/mitochondrialScaffolds.txt 2>> $prefix/filter_contigs.$$.logfile.txt");
  return "$prefix/mitochondrialScaffolds.txt";

  }

sub collectAdditional {

  my ($prefix, @additionals, $reference) = @_;
  system("cat @additionals $reference > $prefix/concatenated.fa 2>> $prefix/cdhit.$$.logfile.txt");
  system("cdhit-est run_$$/concatenated.fa >> $prefix/cdhit.$$.logfile.txt");

  }
1;
