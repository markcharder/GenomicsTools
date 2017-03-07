package PanGenome::CallVariants;
use strict;
use warnings;
use File::Basename;

# Subroutines that act as wrappers for various software used in variant calling. Includes methods for calling SNPs with
# GATK and nucmer. Follows GATK best practices for filtering vcf files.

sub exitStatus {
	# Print message if failed exit status.
	my $exit	= shift @_;
	my $message	= shift @_;
	if ($exit != 0){
		die ($message);
	}
}

sub callVariants {
	# Call variants using GATK's HaplotypeCaller.
	# Output will be filtered according to best practices: https://software.broadinstitute.org/gatk/best-practices/
	# Requires both a dictionary and fasta index for the reference sequence within the same directory.
	my $prefix		= shift @_;
	my $GATKPATH		= shift @_;
	my $referenceGenome	= shift @_;
	my $bam			= shift @_;
	my $baseName		= basename($bam);

	my $exit 		= system("java -jar $GATKPATH/GenomeAnalysisTK.jar -R $referenceGenome -T HaplotypeCaller -I $bam --genotyping_mode DISCOVERY -stand_call_conf 30 -o $prefix/$prefix.$$.vcf >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with HaplotypeCaller.\n");

	$exit			= system("java -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $referenceGenome -V $prefix/$prefix.$$.vcf -selectType SNP -o $prefix/$prefix.$$.raw.vcf >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with SelectVariants.\n");

	$exit			= system("java -jar $GATKPATH/GenomeAnalysisTK.jar -T VariantFiltration -R $referenceGenome -V $prefix/$prefix.$$.raw.vcf --filterExpression \"QD < 2.0 || AF < 1.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0\" --filterName \"snp_filter\" -o $prefix/$prefix.$$.filt.vcf >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with VariantFiltration.\n");

	$exit			= system("java -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $referenceGenome -V $prefix/$prefix.$$.vcf -selectType INDEL -o $prefix/$prefix.$$.indels.raw.vcf >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with SelectVariants.\n");

	$exit = system("java -jar $GATKPATH/GenomeAnalysisTK.jar -T VariantFiltration -R $referenceGenome -V $prefix/$prefix.$$.indels.raw.vcf --filterExpression \"QD < 2.0 || AF < 1.0 || FS > 200.0 || ReadPosRankSum < -20.0\" --filterName \"indel_filter\" -o $prefix/$prefix.$$.indels.filt.vcf >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with VariantFiltration.\n");

	$exit = system("java -jar $GATKPATH/GenomeAnalysisTK.jar -T CombineVariants -R $referenceGenome --variant:a $prefix/$prefix.$$.indels.filt.vcf --variant:b $prefix/$prefix.$$.filt.vcf -o $prefix/$prefix.$$.highqual.vcf -genotypeMergeOptions PRIORITIZE -priority a,b >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with CombineVariants.\n");

	$exit = system(" java -jar $GATKPATH/GenomeAnalysisTK.jar -T BaseRecalibrator -R $referenceGenome -I $bam -knownSites $prefix/$prefix.$$.highqual.vcf -o $prefix/recal.table.$$ >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with BaseRecalibrator.\n");

	$exit = system("java -jar $GATKPATH/GenomeAnalysisTK.jar -T PrintReads -R $referenceGenome -I $bam -BQSR $prefix/recal.table.$$ -o $prefix/$prefix.$$.recal.bam >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with PrintReads.\n");

	$exit = system("java -jar $GATKPATH/GenomeAnalysisTK.jar -R $referenceGenome -T HaplotypeCaller -I $prefix/$prefix.$$.recal.bam --genotyping_mode DISCOVERY -stand_call_conf 30 -o $prefix/$prefix.$$.recal.vcf >> $prefix/gatk.$$.logfile.txt 2>&1");

	$exit = system("java -jar $GATKPATH/GenomeAnalysisTK.jar -R $referenceGenome -T HaplotypeCaller -I $prefix/$prefix.$$.recal.bam --emitRefConfidence GVCF -o $prefix/$prefix.$$.g.vcf >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with HaplotypeCaller.\n");

	return(	"$prefix/$prefix.$$.vcf", 
        	"$prefix/$prefix.$$.vcf.idx", 
		"$prefix/$prefix.$$.raw.vcf.idx", 
		"$prefix/$prefix.$$.filt.vcf.idx", 
		"$prefix/$prefix.$$.indels.raw.vcf.idx", 
		"$prefix/$prefix.$$.indels.filt.vcf.idx", 
		"$prefix/$prefix.$$.highqual.vcf.idx", 
		"$prefix/$prefix.$$.raw.vcf", 
		"$prefix/$prefix.$$.filt.vcf", 
		"$prefix/$prefix.$$.indels.raw.vcf", 
		"$prefix/$prefix.$$.indels.filt.vcf", 
		"$prefix/$prefix.$$.highqual.vcf", 
	 	"$prefix/recal.table.$$", 
		"$prefix/$prefix.$$.g.vcf",
		"$prefix/$prefix.$$.recal.vcf");
}

sub hardFilterVCF {
	# Once the vcf has been produced from the previous step, it can be further hard filtered based on quality.
	# This may be necessary if there is no reference snpDB entry for the organism.
	my $prefix = shift @_;
	my $referenceGenome = shift @_;
	my $GATKPATH = shift @_;
	my $vcf = shift @_;
	my $baseName = basename($vcf);

	my $exit = system("java -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $referenceGenome -V $vcf -selectType SNP -o $prefix/$prefix.hf.$$.raw.vcf >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with SelectVariants.\n");

	$exit = system("java -jar $GATKPATH/GenomeAnalysisTK.jar -T VariantFiltration -R $referenceGenome -V $prefix/$prefix.hf.$$.raw.vcf --filterExpression \"DP < 10.0 || AF < 1.0 || MQ < 30.0\" --filterName \"snp_filter\" -o $prefix/$prefix.hf.$$.filt.vcf >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with VariantFiltration.\n");

	$exit = system("java -jar $GATKPATH/GenomeAnalysisTK.jar -T SelectVariants -R $referenceGenome -V $vcf -selectType INDEL -o $prefix/$prefix.hf.$$.indels.raw.vcf >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with SelectVariants.\n");

	$exit = system("java -jar $GATKPATH/GenomeAnalysisTK.jar -T VariantFiltration -R $referenceGenome -V $prefix/$prefix.hf.$$.indels.raw.vcf --filterExpression \"DP < 10.0 || AF < 1.0 || MQ < 30.0\" --filterName \"indel_filter\" -o $prefix/$prefix.hf.$$.indels.filt.vcf >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with VariantFiltration.\n");

	$exit = system("java -jar $GATKPATH/GenomeAnalysisTK.jar -T CombineVariants -R $referenceGenome --variant:a $prefix/$prefix.hf.$$.indels.filt.vcf --variant:b $prefix/$prefix.hf.$$.filt.vcf -o $prefix/$prefix.hf.$$.highqual.vcf -genotypeMergeOptions PRIORITIZE -priority a,b >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with CombineVariants.\n");

	return( "$prefix/$prefix.hf.$$.raw.vcf",
		"$prefix/$prefix.hf.$$.filt.vcf",
		"$prefix/$prefix.hf.$$.indels.raw.vcf",
		"$prefix/$prefix.hf.$$.indels.filt.vcf",
        	"$prefix/$prefix.hf.$$.raw.vcf.idx",
		"$prefix/$prefix.hf.$$.filt.vcf.idx",
		"$prefix/$prefix.hf.$$.indels.raw.vcf.idx",
		"$prefix/$prefix.hf.$$.indels.filt.vcf.idx",
		"$prefix/$prefix.hf.$$.highqual.vcf");
}

sub genotypeGVCFs {
	# Combines gvcf files into a multi-genotyped GVCF using GATK's GenotypeGVCFs.
	my $prefix		= shift @_;
	my $GATKPATH		= shift @_;
	my $referenceGenome	= shift @_;
	my @GVCFs		= @_;
	foreach my $item (@GVCFs) {
		$item =~ s/$item/--variant $item/g;
	}

	my $exit	= system("java -jar ${GATKPATH}/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $referenceGenome @GVCFs -o $prefix/$prefix.$$.merged.vcf >> $prefix/gatk.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem with GenotypeGVCFs.\n");

	return ("$prefix/$prefix.$$.merged.vcf");
}

sub callWithNucmer {
	# Calls SNPs with nucmer's show-snps utility.
	my $prefix	= shift @_;
	my $reference	= shift @_;
	my $alt		= shift @_;
	my $baseName	= basename($alt);

	my $exit	= system("nucmer --prefix $prefix/$prefix $reference $alt >> $prefix/nucmer.$$.logfile.txt 2>&1");
	exitStatus($exit, "Problem running nucmer.\n");
	$exit		= system("delta-filter -q -r $prefix/$prefix.delta > $prefix/$prefix.filt.delta 2>> $prefix/nucmer.$$.logfile.txt");
	exitStatus($exit, "Problem running delta-filter.\n");
	$exit		= system("grep \">\" $prefix/$prefix.filt.delta | awk \'\{print \$2\}\' | sort -V -u > $prefix/$prefix.present 2>> $prefix/nucmer.$$.logfile.txt");
	exitStatus($exit, "Problem problem getting mapped contigs list.\n");
	$exit		= system("show-snps -Clr $prefix/$prefix.filt.delta > $prefix/$prefix.snps 2>> $prefix/nucmer.$$.logfile.txt");
	exitStatus($exit, "Problem running show-snps.\n");
	return("$prefix/$prefix.snps", "$prefix/$prefix.present"); 
}

sub coverageBedtools {
	# Gets coverage for each position in a reference genome based on a bam file of read mappings.
	my ($prefix, $bam, $reference)	= @_;
	my $exit	= system("bedtools genomecov -d -ibam $bam -g $reference > $prefix/$prefix.$$.coverage 2>> $prefix/bedtools.$$.logfile.txt");
	exitStatus($exit, "Problem getting coverage from bedtools.\n");
	return("$prefix/$prefix.$$.coverage");
}

1;
