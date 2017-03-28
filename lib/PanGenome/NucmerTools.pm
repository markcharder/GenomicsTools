package PanGenome::NucmerTools;
use strict;
use warnings;

sub new {

	my ($class, %args)	= @_;
	return bless {%args}, $class;

}

sub verbose_message{

	my ($self, $message, $verbose)	= @_;
	warn "\n$message\n\n" if $verbose;

}

sub exit_message {

	my ($self, $message_a, $message_b, $verbose)	= @_;

	if($self->{exit} != 0){
		die("\n$message_a\n\n");
	}
	else{
		$self->verbose_message("\n$message_b\n\n", $verbose);
	}

}

sub nucmer_align_filter {

	my ($self, $verbose, $nucmer_args, $filter_args)	= @_;
	my $exit	= system("nucmer  $nucmer_args  $self->{reference} $self->{alt} -p $self->{prefix} >> $self->{prefix}.nucmer.logfile 2>&1");
	$self->{exit}	= $exit;
	bless $self;

	$self->exit_message(	"Problem with Aligning $self->{alt} to $self->{reference} using Nucmer. Check $self->{prefix}.nucmer.logfile for details",
				"Successfully aligned $self->{alt} to $self->{reference} with nucmer with options $nucmer_args", $verbose		);

	$exit		= system("delta-filter  $filter_args  $self->{prefix}.delta > $self->{prefix}.filt.delta 2>> $self->{prefix}.nucmer.logfile");
	$self->{exit}	= $exit;
	bless $self;

	$self->exit_message(	"Problem with filtering $self->{prefix}.delta", 
				"Successfully filtered $self->{prefix}.delta with options $filter_args", $verbose	);

	$self->{delta}	= "$self->{prefix}.filt.delta";
	bless $self;

}

sub show_snps {

	my ($self, $verbose, $args)	= @_;
	my $exit		= system("show-snps $args $self->{delta} > $self->{prefix}.snps 2>>$self->{prefix}.nucmer.logfile");
	$self->{exit}	= $exit;
	bless $self;

	$self->exit_message(	"Problem with running show-snps for $self->{prefix}.delta", 
				"Successfully ran show-snps $self->{prefix}.delta with options $args", $verbose	);

	$self->{snp_file}	= "$self->{prefix}.snps";
	bless $self;

}

sub snps_to_table {

	my ($self, $verbose)	= @_;
	my %hash_headers;
	my @snps;
	my %hash_info;
	my @return_snps;

	open my $FH, "<", $self->{snp_file} or die ("\nCan't open snp file\n\n");
	while (my $line = <$FH>){

        	chomp $line;

        	unless ($line !~ /\s+\d/){

			my @fields = split(/\s+/, $line);
			my @return_fields = split(/\s+/, $line);
			if (!$hash_headers{$fields[17]}){
				my @array	= (1, $fields[9]);
				$hash_headers{$fields[17]}	= \@array;
			}

			my ($rcon, undef, undef)        = split("", $fields[12]);
			my ($acon, undef, undef)        = split("", $fields[13]);
			$hash_info{"$fields[17]-$fields[1]"} = "$rcon\t$acon" if not $hash_info{"$fields[17]-$fields[1]"};
			@fields = ($fields[17], $fields[1], ".", $fields[2], $fields[3], ".", ".", "ORIGIN=denovo\tGT:DP\t1/1:1");
			push @snps, \@fields;
			push @return_snps, \@return_fields;
		}
	}

	($self->{snps}, $self->{headers}, $self->{info}, $self->{all_snp_info})	= (\@snps, \%hash_headers, \%hash_info, \@return_snps);
	bless $self;
}

sub table_to_vcf{

	my ($self, $verbose)	= @_;

	my $snps		= $self->{snps};
	my $hash_info		= $self->{info};
	my $hash_mid_end	= $self->{mid_end};
	my @snps		= @$snps;
	my %drhash		= %$hash_info;

	my ($pname, $cpos, $ppos, $count, $palt, $rinstart)     = ("", "", "", 1, "", "");

	foreach my $val (@snps){
	 	my @array		= @$val;
		my $chrpos		= "$array[0]-$array[1]";
		my ($rcon, $acon)	= split("\t", $drhash{$chrpos}) if $drhash{$chrpos};
		if ($array[3] eq "\."){
			$drhash{$chrpos}  = $drhash{$chrpos} . $array[4] if $drhash{$chrpos};
		}
		elsif ($array[4] eq "\."){
			if (($array[1] != ($ppos + 1))or($array[0] ne $pname)){
				$rinstart       = $array[1];
			}
			$drhash{"$array[0]-$rinstart"}	= $drhash{"$array[0]-$rinstart"} . $array[3] if $drhash{"$array[0]-$rinstart"};
		}
		else {
			$drhash{$chrpos}	= join ("\t", @array);
		}
		$ppos   = $array[1];
		$pname  = $array[0];
	}

	$self->verbose_message("Parsed snps from $self->{snp_file} and formatted for vcf output", $verbose);
	$self->{info}		= \%drhash;
	bless $self;
	my %returnhash;
	foreach my $val (@snps){
		my @array	= @$val;
		my $pos		= $array[1] - 1;
		if ($array[3] eq "\."){
			if ($array[1] ne "$ppos"){
				my ($rcon, $acon)	= split("\t", $drhash{"$array[0]-$array[1]"}) if $drhash{"$array[0]-$array[1]"};
				my @retarray		= ($array[0], $pos, ".", $rcon, $acon, @array[5..$#array]);
				$returnhash{"$array[0]-$pos"}	= \@retarray;
                	}
        	}
		elsif ($array[4] eq "\."){
			$cpos   = $ppos + 1;
			if ($array[1] != "$cpos"){
				if ($drhash{"$array[0]-$array[1]"}){
					my ($rcon, $acon)	= split("\t", $drhash{"$array[0]-$array[1]"});
					my @vals		= split("", $acon);
					$acon			= shift @vals;
					$rcon			= $rcon . join("", @vals);
					my @retarray		= ($array[0], $pos, ".", $rcon, $acon, @array[5..$#array]);
					$returnhash{"$array[0]-$pos"}	= \@retarray;
				}
			}
		}	
		else{
			$returnhash{"$array[0]-$array[1]"}	=  \@array;
		}
	$ppos   = $array[1];
	}
	$self->verbose_message("Created vcf fields from snps in $self->{snp_file} using 'ORIGIN=denovo' and 'GT:DP 1/1:1'\nChanging these options is not currently supported though it can be achieved by modifying this script", $verbose);
	$self->{vcf_fields}	= \%returnhash;
	bless $self;

}

sub print_vcf {

	my ($self, $verbose, $species)	= @_;
	my $headers			= $self->{headers};
	my %drheaders			= %$headers;
	my @time			= localtime();
	my $year			= $time[5] += 1900;
	my $month			= sprintf("%02d", $time[4]);
	my $day				= sprintf("%02d", $time[3]);

	open my $FH, ">", "$self->{prefix}.nucmer.vcf" or die("\nProblem writing to $self->{prefix}.nucmer.vcf\n\n");
	print $FH 	"##fileformat=VCFv4.2\n##fileDate=$year$month$day\n" .
			"##source=GenomicsToolsV1\n##reference=$self->{reference}\n";
	for my $key (keys %drheaders){
		my @array	= @{$drheaders{$key}};
		print $FH 	"##contig<ID=$key,length=$array[1],species=\"$species\">\n";
	}
	print $FH	"##INFO=<ID=ORIGIN,Number=1,Type=String,Description=\"How this snp was imputed\">\n" .
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\"\n" .
			"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n" .
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$self->{alt}>\n";
	my $vcf_fields	= $self->{vcf_fields};
	my %drvcf	= %$vcf_fields;

	for my $key (sort{no warnings; my ($aa) = $a =~ /(\d+)/; my ($bb) = $b =~ /(\d+)/; $aa <=> $bb || $drvcf{$a}->[1] <=> $drvcf{$b}->[1]} keys %drvcf){
		my @array	= @{$drvcf{$key}};
		unless (($array[3] =~ /R|Y|S|W|K|M|B|D|H|V|N|\./) or ($array[4] =~ /R|Y|S|W|K|M|B|D|H|V|N|\./)){
			print $FH join("\t", @array) . "\n";
		}
	}
	close $FH;
}

sub fix_vcf_nucmer {

	my ($self, $verbose, $vcf) = @_;
	my $hash	= $self->{vcf_fields};
	my %drhash	= %$hash;
	$self->verbose_message("Fixing nucmer vcf with supplied vcf", $verbose);
	my %return_hash;
	open my $FH, "<", $vcf or die ("\nCan't open $vcf\n\n");
	while (my $line = <$FH>){
		if ($line =~ /^#/){
			next;
		}
		chomp $line;
		my @fields	= split("\t", $line);
		my $chrpos	= "$fields[0]-$fields[1]";
		if ($drhash{$chrpos}){
			my @array	= @{$drhash{$chrpos}};
			if ($array[3] ne $fields[3]){
				if (($fields[6] eq "PASS") or ($fields[6] eq "\.")){
					$array[3] =~ s/$array[3]/$fields[3]/g;
					$array[4] =~ s/$array[4]/$fields[4]/g;
					$drhash{$chrpos}	= \@array;
				}
				else{
					delete $drhash{$chrpos};
				}
			}
		}
		my $genotype	= $fields[8];
		if(eof){
			foreach my $key (keys %drhash){
				my @array		= @{$drhash{$key}};
				my @array_genotype	= split("\t", $array[7]);
				if ($array_genotype[1] ne $genotype){
					$array_genotype[1] =~ s/\Q$array_genotype[1]/$genotype/g;
					my @genotype_fields	= split(":", $genotype);
					my @genotype_info;
					for my $i (@genotype_fields){
						push @genotype_info, 1;
					}
					$genotype_info[0]	= "1/1";
					my $newinfo		= join(":", @genotype_info);
					$array[7] 		=~ s/\Q$array[7]/ORIGIN=denovo\t$genotype\t$newinfo/g;
					my @newarray		= ($array[0], $array[1], $array[2], $array[3], $array[4], $array[5], $array[6], $array[7]);
					$drhash{$key}		= \@newarray;
				}
			}
		}
		my @line		= split (/\s+/, $line);
		$return_hash{$chrpos}	= \@line;
	($self->{vcf_fields}, $self->{vcf_file})	= (\%drhash, \%return_hash);
	bless $self;
	}
}

sub filter_low_covered {

	my ($self, $verbose, $cov, $dp)	= @_;
	my $vcf_file			= $self->{vcf_fields};
	my %drvcf			= %$vcf_file;
	my %coverage_hash;
	my $total			= 0;
	my @coverage			= ();
	my %return_hash;
	$self->verbose_message("Reading in bedtools coverage file", $verbose);
	open my $FH, "<", $cov or die ("\nBedtools coverage file not readable\n\n");
	while(my $line	= <$FH>){
		chomp $line;
		my @fields	= split (/\s+/, $line);
		if ($drvcf{"$fields[0]-$fields[1]"}){
			$coverage_hash{"$fields[0]-$fields[1]"} = $fields[2];
		}
		push @coverage, $fields[2];
		$total += $fields[2];
	}
	my $mean	= int($total / @coverage);
	close $FH;
	foreach my $key (keys %drvcf){
		my @fields	= @{$drvcf{$key}};
		if ($coverage_hash{"$fields[0]-$fields[1]"}){
			if ($coverage_hash{"$fields[0]-$fields[1]"} < int($dp * $mean)){
				$return_hash{"$fields[0]-$fields[1]"}	= \@fields;
			}
		}
	}
	my $coverage_threshold	= int($dp * $mean);
	$self->verbose_message("Removed de novo SNPs with at least $coverage_threshold x coverage from reads", $verbose);
	$self->{vcf_fields}	= \%return_hash;
	bless $self;
}
1;
