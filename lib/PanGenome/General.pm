package PanGenome::General;
use strict;
use warnings;

sub parseConfig {

  my ($UTILPATH,
  $SAMTOOLSPATH, 
  $REPEATMODELERPATH, 
  $QUARRYPATH, 
  $REPEATSCOUTPATH,
  $TRFPATH,
  $NSEGPATH,
  $RMBLASTPATH,
  $RECONPATH,
  $REPEATMASKERPATH,
  $BAM2FASTQPATH,
  $A5PATH,
  $STAMPYPATH,
  $BWAPATH,
  $MUMMERPATH,
  $GATKPATH,
  $PICARDPATH,
  $CODINGQUARRYPATH,
  $EXONERATEPATH,
  $MAKERPATH);

  open my $fh, "<", "config.txt";
  while (my $line = <$fh>){

    if ($line =~ /^#/){
      next;
      }

    chomp $line;
    my @fields = split /=/, $line;

    if ($line =~ /^QUARRY_PATH/){
      $QUARRYPATH = $fields[1];
      $ENV{QUARRY_PATH} = $QUARRYPATH;
      }

    if ($line =~ /^REPEATMODELER_PATH/){
      $REPEATMODELERPATH = $fields[1];
      }

    if ($line =~ /^BIN_PATH/){
      $UTILPATH = $fields[1];
      }

    if ($line =~ /^MAKER_PATH/){
      $MAKERPATH = $fields[1];
      }

    if ($line =~ /^CODINGQUARRY_PATH/){
      $CODINGQUARRYPATH = $fields[1];
      }

    if ($line =~ /^RMBLAST_PATH/){
      $RMBLASTPATH = $fields[1];
      }

    if ($line =~ /^NSEG_PATH/){
      $NSEGPATH = $fields[1];
      }

    if ($line =~ /^TRF_PATH/){
      $TRFPATH = $fields[1];
      }

    if ($line =~ /^PICARD_PATH/){
      $PICARDPATH = $fields[1];
      }

    if ($line =~ /^REPEATMASKER_PATH/){
      $REPEATMASKERPATH = $fields[1];
      }

    if ($line =~ /^RECON_PATH/){
      $RECONPATH = $fields[1];
      }

    if ($line =~ /^REPEATSCOUT_PATH/){
      $REPEATSCOUTPATH = $fields[1];
      }

    if ($line =~ /BAM2FASTQ_PATH/){
      $BAM2FASTQPATH = $fields[1];
      }

    if ($line =~ /^GATK_PATH/){
      $GATKPATH = $fields[1];
      }

    if ($line =~ /^A5_PATH/){
      $A5PATH = $fields[1];
      }

    if ($line =~ /^STAMPY_PATH/){
      $STAMPYPATH = $fields[1];
      }

    if ($line =~ /^BWA_PATH/){
      $BWAPATH = $fields[1];
      }

    if ($line =~ /^MUMMER_PATH/){
      $MUMMERPATH = $fields[1];
      }

    if ($line =~ /^EXONERATE_PATH/){
      $EXONERATEPATH = $fields[1];
      }

    if ($line =~ /^SAMTOOLS_PATH/){
      $SAMTOOLSPATH = $fields[1];
      }

    foreach my $path ($CODINGQUARRYPATH, 
                      $SAMTOOLSPATH, 
                      $REPEATMODELERPATH, 
                      $REPEATSCOUTPATH, 
                      $TRFPATH, 
                      $NSEGPATH, 
                      $RMBLASTPATH, 
                      $RECONPATH, 
                      $REPEATMASKERPATH, 
                      $BAM2FASTQPATH, 
                      $A5PATH, 
                      $STAMPYPATH, 
                      $BWAPATH, 
                      $MUMMERPATH, 
                      $GATKPATH, 
                      $PICARDPATH,
                      $EXONERATEPATH,
		      $MAKERPATH,
		      $UTILPATH){
      if ($path){
        $ENV{PATH} = "$ENV{PATH}:$path";
        }
      }
    }

  return ($GATKPATH, $PICARDPATH, $REPEATMASKERPATH, $EXONERATEPATH);
  close $fh;

  }

sub outputDir {

  my $output = shift @_;
  if ( ! -d $output) {
    system("mkdir $output");
    }
  }

sub parseOptions {

  my $readsFiles = shift @_;
  my @fields = split /,/, $readsFiles;

  return @fields;

  }


sub verboseMessage {

  my $message = shift @_;
  my $verbose = shift @_;

  if ($verbose) {
    print "$message\n";
    }

  }

sub cleanUp {

  system ("rm -r @_");

  }

sub getFastaSeqs {

  my $list = shift @_;
  my $fasta = shift @_;
  my $inverse = shift @_;
  my $out = shift @_;
  my %names;
  my $switch = 0;

  open my $fh, "<", $list or die ("\nError with filtering mitochondrial.\n");
  while (my $line = <$fh>) {
    chomp $line;
    $names{$line} = "";
    }
  close $list;

  open $fh, "<", $fasta or die ("\nError with filtering mitochondrial.\n");
  while (my $line = <$fh>) {
    chomp $line;

    if ($line =~ />/) {
      $switch = 0;
      $line =~ s/>//g;
      if ($inverse eq 'inverse'){
        if (!defined($names{$line})) {
          $switch = 1;
          }
        }
      if ($inverse eq 'n'){
        if (defined($names{$line})) {
          $switch = 1;
          }
        }
      if ($switch == 1) {
        open my $fho, ">>", $out;
        print $fho ">$line\n";
        close $fho;
        }
      }
    else {
      if ($switch == 1) {
        open my $fho, ">>", $out;
        print $fho "$line\n";
        close $fho;
        }
      }
    }
    close $fh;
  }

1;
