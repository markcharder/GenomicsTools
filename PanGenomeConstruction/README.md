# PanGenomeConstruction

## Perl package used for 'pan genome' (if you like that word) assembly and analysis.

### Contents:
-	Requirements
-	Installation
-	Scripts
-	Notes on running on slurm cluster

#### Requirements:

-	Maker2				
-	Coding Quarry			v2.0
-	Repeat Modeler			
-	Repeat Masker
-	Tandem Repeat Finder
-	Samtools			v0.1.19
-	Genome Analysis Toolkit		v3.7
-	Picard Tools			v2.9.0
-	Bed Tools			v2.25.0
-	Bioperl
-	Mummer
-	A5 (assembly pipeline)
-	Stampy
-	BWA
-	Exonerate
-	Repeat Scout
-	Bam2Fastq
-	NSEG
-	RMBLAST
-	Repeat Scout
-	Recon

Versions listed have been tested but are not necessarily the only versions the package will work with.

#### Installation:

Scripts and modules are generally wrappers for other softwares. 
Paths to these softwares can either be specified in config.txt or added to $PATH.
GATK and Picard paths must be specified in the config file, bin/config.txt. 
In order to use CodingQuarry, the QUARRY_PATH variable must also be set in bin/config.txt.

Not all of these software packages are easy to install or configure. Unfortunately this package
does not come with an automatic install script.
Other than this, the package should run straight out of the box.

	
#### Scripts:

##### assembly

This script assembles short reads de novo using the A5 pipeline for automated parameter optimisation 
(http://bioinformatics.oxfordjournals.org/content/31/4/587.full.pdf).

There are three typical usage cases for this script.

	1)	assembly --reads /path/to/reads.R1.fq,/path/to/reads.R2.fq
	2)	assembly --bam /path/to/bam_file.bam
	3)	assembly --mitochondrial /path/to/mitochondrial_scaffolds.fa --contigs /path/to/contigs.fa

1) This will assemble reads from the fastq files provided and is simply a wrapper for A5.
It will automate generation of a library file (.lib) for use with the A5 script.

2) This will gather unmapped reads from a bam file and assemble them using A5. This can be used to produce
contigs that are not present in the reference sequence that the reads were mapped to originally.

3) This will determine contigs that are likely to be mitochondrial based on comparison to known mitochondrial sequence using nucmer.

Software used in assembly:

	Software	Version		Link
	--------	-------		----

	A5		20160825	http://gsl.hudsonalpha.org/information/software/bam2fastq
	Bam2fastq	1.1.0		https://github.com/jts/bam2fastq
	Nucmer		3.1		http://mummer.sourceforge.net/

###### annotation

This script annotates assembled sequences using Coding Quarry, Repeat Masker and Repeat Modeler, and Maker2.

There is one typical usage case for this script.

	1)	annotation --type all --contigs /path/to/contigs.fa --prot /path/to/proteins.fa --species species_name.

1) This will run Coding Quarry,  Repeat Modeler, Repeat Masker and Maker2 on the supplied contigs. The '--prot' option supplies a protein fasta file to maker that will be mapped to the genome using exonerate to provide protein evidence. This is any known high quality set of proteins. This could include for example reviewed proteins from SwissProt. The '--species' option is a species that Coding Quarry has been trained on. To use Coding Quarry or Maker with RNA sequencing evidence, this script should be modified.

Software used in annotation:

	Software	Version		Link
	--------	-------		----
	
	Coding Quarry	2.0		https://sourceforge.net/projects/codingquarry/
	Repeat Modeler	1.0.4		http://www.repeatmasker.org/RepeatModeler.html
	Repeat Masker	4.0.6		http://www.repeatmasker.org/
	Maker2		2.31.8		http://www.yandell-lab.org/software/maker.html

##### mapping

This script maps short reads to a reference sequence with the stampy software. It subsequently fixes the bam files using picard tools.

There is one typical usage case for this script.

	1)	mapping --reference /path/to/reference.fa --reads /path/to/reads.R1.fq,/path/to/reads.R2.fq

1) This will map the paired end reads to the genome using stampy in conjunction with BWA. It will then reformat the resulting bam file according to GATK specifications, adding read groups etc. It will also create a bam index.

Software used in mapping:

	Software	Version		Link
	--------	-------		----

	BWA		0.7.5a-r405	http://bio-bwa.sourceforge.net/
	Samtools	0.1.19-96b5f2294a	http://samtools.sourceforge.net/
	Picardtools	2.0.1		https://broadinstitute.github.io/picard/
	Stampy		1.0.29		http://www.well.ox.ac.uk/stampy

##### callVariants

This script uses supplied bam files to call variants using gatk. It follows gatk best practices and will hard filter vcf files.

There are three typical usage cases for this script.

	1)	callVariants --reference /path/to/reference.fa --bam /path/to/bam_file.bam
	2)	callVariants --reference /path/to/reference.fa --GVCFs /path/to/gvcf1.g.vcf,/path/to/gvcf2.g.vcf,/path/to/gvcf3.g.vcf
	3)	callVariants --nucmer --filter-snps --depth 0 --bam /path/to/bam_file.bam --reference /path/to/reference.fa --denovo /path/to/alternate_genome.fa

1) This will map reads to a reference genome and call SNPs between the reference and the reads. It will hard filter the variants based on depth, quality and allele frequency. It assumes a haploid, homokaryotic genome, so allele frequency can only be 1. It will produce both a vcf and gvcf file for the input bam.

2) This will combine a series of GVCF files into a single vcf, applying the same hard filter as in 1).

3) This will call variants in regions where assembled contigs align but there is low or no read mapping (specified by '--depth'). This is to include highly variable regions where reads do not map but that are nonetheless a contiguous part of the alternate assembly and not an insertion or deletion relative to the reference.

Software used in callVariants:

	Software	Version		Link
	--------	-------		----
	GATK		3.6-0-g89b7209	https://software.broadinstitute.org/gatk/download/
	Picardtools	2.0.1		https://broadinstitute.github.io/picard/
	Samtools	0.1.19-96b5f2294a		http://samtools.sourceforge.net/
	Bedtools	v2.17.0		http://bedtools.readthedocs.io/en/latest/content/installation.html
	Nucmer		3.1		http://mummer.sourceforge.net/

#### Notes on running on a slurm cluster

In its original implementation, this package was run on a slurm cluster using a reference genome and illumina reads for alternate isolates.
The whole process, encompassing all scripts, was run from beginning to end in parallel using mpibash with 24 internal threads for each instance.

The scripts used to accomplish this are in the bash directory.

To run the whole process takes some time - up to a week for a small genome of only 40 Mb. As there is a walltime limit on the supercomputer that
was used in the original implementation, the process was split across several jobs.
