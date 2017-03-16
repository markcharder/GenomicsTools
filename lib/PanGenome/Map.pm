package PanGenome::Map;
use strict;
use warnings;
use File::Basename;

# Subroutines for mapping reads and formatting bam files. Uses Stampy and bwa.

sub exitStatus {
	# Die with error message on failed exit status.
	my $message	= shift @_;
	my $exit	= shift @_;
	if ($exit != 0){
		die ($message);
	}
}

sub filterMultiMap {
	# Filters out reads from a bam file that have a mapping quality less than 30.
	my $prefix = shift @_;
	my $bam = shift @_;

	my $exit = system("samtools view -q 30 -b $bam -o $prefix/$prefix.singleSites.bam >> $prefix/samtools.$$.logfile.txt 2>&1");
	exitStatus("Problem removing multiple mapping reads with samtools\n", $exit);
	return "$prefix/$prefix.singleSites.bam";
}

sub alignBwa {
	# Aligns reads with bwa.
	my $u_out;
	my $p_out;
	my $prefix = shift @_;
	my $referenceGenome = shift @_;
	my $threads = shift @_;
	my @files = @_;
	my $refbase = basename($referenceGenome);

	my $exit = system("bwa index $referenceGenome -p $prefix/$refbase >> $prefix/bwa.$$.logfile.txt 2>&1");
	exitStatus("Problem indexing with bwa.\n", $exit);
	$exit = system("ln -s $referenceGenome $prefix");
	exitStatus("Problem indexing with bwa.\nMaybe sym-links not permitted on drive?\n", $exit);

	foreach my $file (@files) {
		my $fileBase = basename($file);
		my $fileDir = dirname($file);
		$exit = system("bwa aln -q10 -t $threads $prefix/$refbase $file  -f $prefix/$fileBase.sai >> $prefix/bwa.$$.logfile.txt 2>&1");
		exitStatus("Problem  creating sam index with bwa.\n", $exit);
		$exit = system("ln -s $file $prefix");
		exitStatus("Problem creating sam index with bwa.\nMaybe sym-links not permitted on drive?\n", $exit);
	}


	if ((@files == 2) or (@files == 3)) {
		my $fileBaseF = basename($files[0]);
		my $fileBaseR = basename($files[1]);
		my $fileBaseO = basename($files[0]);
		my $fileDir = dirname($files[0]);
		$p_out = "$prefix/$prefix.$$.sam";

		$exit = system("bwa sampe $prefix/$refbase $prefix/$fileBaseF.sai $prefix/$fileBaseR.sai $files[0] $files[1] -f $prefix/$prefix.$$.sam >> $prefix/bwa.$$.logfile.txt 2>&1");
		exitStatus("Problem  mapping to reference with bwa.\n", $exit);
    
	}

	if (@files == 3) {
		my $fileBase = basename($files[2]);
		my $fileBaseO = basename($files[2]);
		my $fileDir = dirname($files[0]);

		$exit = system("bwa samse $prefix/$refbase $prefix/$fileBase.sai $files[2] -f $prefix/$prefix.$$.u.sam >> $prefix/bwa.$$.logfile.txt 2>&1");
		$u_out = "$prefix/$prefix.$$.u.sam";
		exitStatus("Problem  mapping to reference with bwa.\n", $exit);
	}

	if (@files == 1) {
		my $fileBase = basename($files[0]);
		my $fileBaseO = basename($files[0]);
		my $fileDir = dirname($files[0]);
		$u_out = "$prefix/$prefix.$$.u.sam";

		$exit = system("bwa samse $prefix/$refbase $prefix/$fileBase.sai $files[0] -f $prefix/$prefix.$$.u.sam >> $prefix/bwa.$$.logfile.txt 2>&1");
		exitStatus("Problem  mapping to reference with bwa.\n", $exit);
	}

	return($p_out,$u_out);
}


sub samToBam {
	# Converts sam files to bam files.
	my $prefix = shift @_;
	my $referenceGenome = shift @_;
	my $threads = shift @_;
	my $sam =  shift @_; 

	my $baseName = basename($sam);
	$baseName =~ s/(.*)\..*/$1/;

	my $exit = system("samtools view -@ $threads -bT $referenceGenome $sam > $prefix/$baseName.bam 2>> $prefix/samtools.$$.logfile.txt");
	exitStatus("Error converting sams to bams.\n", $exit);

	return "$prefix/$prefix.$$.bam";
}


sub mergeBams {
	# Merges bam files.
	my $prefix = shift @_;
	my @bams =  @_;
	my $baseName = basename($bams[0]);

	my $exit = system("samtools merge $prefix/$prefix.merged.bam @bams >> $prefix/samtools.$$.logfile.txt 2>&1");
	return "$prefix/$prefix.merged.bam";
	exitStatus("Problem merging bams.\n", $exit);
}


sub sortBam {
	# Sorts bam files.
	my $prefix = shift @_;
	my $bam = shift @_;
	my $fileBase = basename($bam);

	my $exit = system ("samtools sort $bam $prefix/$prefix.$$.sort >> $prefix/samtools.$$.logfile.txt 2>&1");
	exitStatus("Problem sorting bams.\n", $exit);

	return "$prefix/$prefix.$$.sort.bam";
}


sub alignStampy {
	# Aligns reads based on bwa mappings using Stampy.
	my $prefix = shift @_;
	my $referenceGenome = shift @_;
	my $threads = shift @_;
	my $bamFile = shift @_;
	my $fileBase = basename($referenceGenome);
	my $bamBase = basename($bamFile);

	my $exit = system("stampy.py -G $prefix/$prefix.$$ $referenceGenome; stampy.py -g $prefix/$prefix.$$ -H $prefix/$prefix.$$ >> $prefix/stampy.$$.logfile.txt 2>&1");
	$exit = system("stampy.py -g $prefix/$prefix.$$ -h $prefix/$prefix.$$ -t $threads --bamkeepgoodreads -M $bamFile -o $prefix/$prefix.$$.stampy.sam >> $prefix/stampy.$$.logfile.txt 2>&1");
	exitStatus("Problem mapping with stampy.\n", $exit);

	return "$prefix/$prefix.$$.stampy.sam";
}


sub fixBam {
	# Fixes badly formatted bam files in accordance with GATK recommendations:
	# http://gatkforums.broadinstitute.org/gatk/discussion/2909/how-to-fix-a-badly-formatted-bam
	my $prefixo = shift @_;
	my $PICARDPATH = shift @_;
	my $reference = shift @_;
	my $bam = shift @_;
	my $baseName = basename($bam);
	my @fileParts = split /\./, $reference;
	my $prefix = $fileParts[0];
	my $exit = system("java -jar $PICARDPATH/picard.jar CleanSam INPUT=$bam OUTPUT=$prefixo/$prefixo.$$.filt.bam >> $prefixo/picard.$$.logfile.txt 2>&1");
	exitStatus("Problem with fixing sams for picard.\n", $exit);

	if (! -e "$prefix.dict"){
		$exit = system("java -jar $PICARDPATH/picard.jar CreateSequenceDictionary R=$reference O=$prefix.dict >> $prefixo/picard.$$.logfile.txt 2>&1");
		exitStatus("Problem creating picard sequence dictionary.\n", $exit);
	}

	$exit = system("java -jar $PICARDPATH/picard.jar SortSam I=$prefixo/$prefixo.$$.filt.bam O=$prefixo/$prefixo.$$.sort.bam SORT_ORDER=coordinate >> $prefixo/picard.$$.logfile.txt 2>&1");
	exitStatus("Problem sorting sam with picard.\n", $exit);

	$exit = system("java -jar $PICARDPATH/picard.jar MarkDuplicates I=$prefixo/$prefixo.$$.sort.bam O=$prefixo/$prefixo.$$.md.bam M=$prefixo/$prefixo.$$.metrics >> $prefixo/picard.$$.logfile.txt 2>&1");
	exitStatus("Problem marking duplicates with picard.\n", $exit);

	$exit = system("java -jar $PICARDPATH/picard.jar AddOrReplaceReadGroups I=$prefixo/$prefixo.$$.md.bam O=$prefixo/$prefixo.$$.rg.bam RGID=$prefixo RGLB=lib.$prefixo RGPL=Illumina RGPU=Unit.$prefixo RGSM=sample.$prefixo >> $prefixo/picard.$$.logfile.txt 2>&1");
	exitStatus("Problem adding read groups with picard.\n", $exit);

	$exit = system("samtools index $prefixo/$prefixo.$$.rg.bam >> $prefixo/samtools.$$.logfile.txt 2>&1");
	exitStatus("Problem adding read groups with picard.\n", $exit);

	return ("$prefixo/$prefixo.$$.rg.bam",
		"$prefixo/$prefixo.$$.md.bam",
		"$prefixo/$prefixo.$$.sort.bam",
		"$prefixo/$prefixo.$$.filt.bam");
}

1;
