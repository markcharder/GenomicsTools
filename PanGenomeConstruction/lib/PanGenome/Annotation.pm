package PanGenome::Annotation;
use strict;
use warnings;
use Cwd;
use File::Basename;

# Subroutines for annotating a genome sequence using conding quarry, maker and repeat masker/modeler.

sub exitStatus {
	# Routine for dying with error upon failed exit status.
	my $exit	= shift @_;
	my $message	= shift @_;
	if ($exit != 0){
		die ($message);
	}
}

sub annotateCodingQuarry {
	# Annotate genome using conding quarry without RNASeq but trained on a related species with RNASeq.
	# This is a very specific usage case. To run coding quarry differently, it can be run independenty
	# and the result can be supplied to 'bin/annotate' as a parameter.
	my ($species, $threads, $contigs, $prefix) = @_;
	my $dir		= getcwd;
	chdir $prefix;
	my $exit	= system("CodingQuarry -p $threads -f $contigs -s $species >> CodingQuarry.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem running CodingQuarry.\n");
	$exit		= system("mv out CodingQuarry.$$ >> CodingQuarry.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem running CodingQuarry.\n");
	chdir $dir;
	return ("$prefix/CodingQuarry.$$/PredictedPass.gff3");
}


sub annotateRepeats {
	# Annotate repeats using repeat modeler and repeat masker. Again, can be run separately and results can be supplied 
	# as parameters to 'bin/annotate'.
	my $dir		= getcwd;
	my ($contigs, $threads, $REPEATMASKERPATH, $prefix) = @_;
	my $contigsBase	= basename($contigs);
	chdir $prefix;

	my $exit	= system("BuildDatabase -name $prefix.$$.db -engine ncbi $contigs >> RepeatMasker.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem running BuildDatabase of Repeat Modeler.\n");

	$exit		= system("RepeatModeler -database $prefix.$$.db >> RepeatMasker.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem running Repeat Modeler.\n");
	chdir $dir;

	if (! -e "$REPEATMASKERPATH/Libraries/RepeatMaskerLib.fasta") {
	$exit		= system("$REPEATMASKERPATH/util/buildRMLibFromEMBL.pl $REPEATMASKERPATH/Libraries/RepeatMaskerLib.embl > $REPEATMASKERPATH/Libraries/RepeatMaskerLib.fasta 2>> RepeatMasker.$$.logfile.txt");
	exitStatus($exit, "Problem running Repeat Masker conversion of Lib file.\n");
	}

	$exit		= system("cat $REPEATMASKERPATH/Libraries/RepeatMaskerLib.fasta $prefix/RM_*/consensi.fa.classified > $prefix/RepeatMaskerCombinedLib.$$.fasta 2>> $prefix/RepeatMasker.$$.logfile.txt");
	exitStatus($exit, "Problem finding Repeat Modeler and Repeat Masker libraries to concatenate.\n");

	chdir $prefix;
	system("mkdir RepeatMasker");
	chdir "RepeatMasker";
	$exit		= system("RepeatMasker -xsmall -gff -s -dir ./ -lib ../RepeatMaskerCombinedLib.$$.fasta $contigs >> ../RepeatMasker.$$.logfile.txt");
	chdir $dir;
	exitStatus($exit, "Problem running Repeat Masker.\n");

	my $baseName	= basename($contigs);
	return "$prefix/RepeatMasker/$baseName.out";
}


sub reformatRepeatMasker {
	# Reformats repeat masker output to gff3 for use with maker.
	my ($REPEATMASKERPATH, $output, $prefix) = @_;
	my $exit	= system("perl $REPEATMASKERPATH/util/rmOutToGFF3.pl $output > $prefix/$prefix.$$.rm.reformat.gff 2>> $prefix/RepeatMasker.$$.logfile.txt");
	exitStatus($exit, "Problem converting RepeatMasker output to gff3.\n");
	$exit		= system("awk 'BEGIN{OFS=\"\\t\";c=1}{if(\$0!~/^#/){print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\"ID=\"c\";\"\$9\" \"\$10\" \"\$11;c++}else{print \$0}}' $prefix/$prefix.$$.rm.reformat.gff > $prefix/temp 2>> $prefix/RepeatMasker.$$.logfile.txt");
	exitStatus($exit, "Problem converting RepeatMasker output to gff3.\n");
	$exit		= system("mv $prefix/temp $prefix/$prefix.$$.rm.reformat.gff >> $prefix/RepeatMasker.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem converting RepeatMasker output to gff3.\n");
	return "$prefix/$prefix.$$.rm.reformat.gff";
}

sub createMakerConfigs {
	# Automates creation of maker config files.
	my ($genome, $EXONERATEPATH, $REPEATMASKERPATH, $prot, $repeats, $abinitio, $prefix) = @_;
	open OUT, ">>$prefix/maker_bopts.ctl";
	print OUT "#-----BLAST and Exonerate Statistics Thresholds
blast_type=ncbi+
pcov_blastn=0.8
pid_blastn=0.85
eval_blastn=1e-10
bit_blastn=40
depth_blastn=0
pcov_blastx=0.5
pid_blastx=0.8
eval_blastx=1e-06
bit_blastx=30
depth_blastx=0
pcov_tblastx=0.8
pid_tblastx=0.85
eval_tblastx=1e-10
bit_tblastx=40
depth_tblastx=0
pcov_rm_blastx=0.5
pid_rm_blastx=0.4
eval_rm_blastx=1e-06
bit_rm_blastx=30
eva_pcov_blastn=0.8
eva_pid_blastn=0.85
eva_eval_blastn=1e-10
eva_bit_blastn=40
ep_score_limit=20
en_score_limit=20";
	close OUT;
	open OUT, ">>$prefix/maker_opts.ctl";
	print OUT "#-----Genome
genome=$genome
organism_type=eukaryotic
#-----Protein Homology Evidence
protein=$prot
#-----Repeat Masking
rm_gff=$repeats
#-----Gene Prediction
pred_gff=$abinitio";
	close OUT;
	open OUT, ">>$prefix/maker_exe.ctl";
	print OUT "#-----Location of Executables Used by MAKER/EVALUATOR
exonerate=$EXONERATEPATH/exonerate
RepeatMasker=$REPEATMASKERPATH/RepeatMasker";
	close OUT;
	return ("$prefix/maker_opts.ctl", "$prefix/maker_bopts.ctl", "$prefix/maker_exe.ctl");
}


sub makerAnnotate {
	# Runs maker using previous config files. The proteins supplied to exonerate for use with maker can be
	# high quality protein annotations from closely related species or strains.
	my ($prefix, $contigs, $opt, $bopt, $exe) = @_;
	my $dir		= getcwd();
	my $baseName	= basename($contigs);
	chdir $prefix;
	system("maker -fix_nucleotides $opt $bopt $exe >> Maker.$$.logfile.txt 2>&1");
	chdir $dir;
	system("mv $prefix/$prefix.maker.output $prefix/$prefix.$$.maker.output");
	return "$prefix/$prefix.$$.maker.output";
}

1;

