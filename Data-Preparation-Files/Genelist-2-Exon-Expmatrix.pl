# SpliceNet Input Generator #
# Hari Krishna Y 20th Nov 2013 #
# Panwen Wang 6th Jan 2014 #
use Data::Dumper;
use strict;

my $GENELIST=$ARGV[0];

# Get all Gene symbol, NM ID and UCSC ID mapping#
my %ucsc2refseq;
my %refseq2ucsc;
my %ucsc2gsym;
my %refseq2gsym;

open (G, $GENELIST);
my @genes=<G>;
foreach my $g (@genes) 
{
	chomp $g;
	#print "grep \"\\s$g\\s\" UCSC-GENE-SYM-RefSeQ";
	my $str = `grep "\\s$g\\s" UCSC-GENE-SYM-RefSeQ`;
	#print $str;
	foreach my $line (split /\n/, $str) 
	{
		chomp $line;
		next if (!$line);
		my ($ucsc,$gene,$refseq) = split /\t/, $line;
		next if (!$refseq);
        my ($rucsc)=split /\./,$ucsc;
		$refseq2ucsc{$refseq} = $rucsc;
		#$ucsc2refseq{$ucsc} = $refseq;
	}
}
close(G);

#print scalar (keys %refseq2ucsc) . "\n";

#open(FW1,">$GENELIST"."-Exon-Exp-Matrix-NC.txt");
open(FW2,">$GENELIST"."-Exon-Exp-Matrix-C.txt");
foreach my $refseq (keys %refseq2ucsc) 
{
	my $str = `grep "^$refseq\\s" Map-Isoform-2-Exon-2`;
	#
	# NM_001025603	chr1:151315818-151315919:-	2	NM_001025603
	# NM_001025603	chr1:151316157-151316359:-	2	NM_001025603
	# NM_001025603	chr1:151316673-151316755:-	2	NM_001025603,NM_001025604
	#
	chomp $str;
	next if (!$str);
	#print $str;
	
	foreach my $line (split /\n/, $str) 
	{
		chomp $line;
		next if (!$line);
		my ($ref, $exon, $count, $share) = split /\t/, $line;
		my ($chrom, $coord, $strand) = split /:/, $exon;
		my ($start, $end) = split /-/, $coord;
		#my $cmd = "awk '{ split(\$1, tmp, \":\"); if (tmp[1]==\"$chrom\" && tmp[3]==\"$strand\") { split(tmp[2], coord, \"-\"); dstart = $start-coord[1]; dend = $end-coord[2]; if (dstart*dend<0 || (dstart<=5 && dstart>=-5) && (dend<=5 && dend>=-5)) print FILENAME\"\\t\"\$1\"\\t\"\$4;} }' *.exon_quantification.txt\n";
		my $cmd = "awk '{ split(\$1, tmp, \":\"); if (tmp[1]==\"$chrom\" && tmp[3]==\"$strand\") { split(tmp[2], coord, \"-\"); if (($start<=coord[1] && $end>=coord[2]) || ($start-coord[1]<=5 && coord[2]-$end<=5) || ($start<=coord[1] && coord[2]-$end<=5) || ($start-coord[1]<=5 && $end>=coord[2])) print FILENAME\"\\t\"\$1\"\\t\"\$4;} }' *.exon_quantification.txt";
		my $str2 = `$cmd`;
		#print $cmd . "\n";
		#print $str2 . "\n";
		#
		# TCGA-A7-A13E-11A-61R-A12P-07.exon_quantification.txt    chr1:11874-12227:+      0.345215972174439
		# TCGA-A7-A13E-11A-61R-A12P-07-2.exon_quantification.txt    chr1:11874-12227:+      0.345215972174439
		#
		#print $str2;
		my %sharediso = ();
		foreach my $sh (split /,/, $share) { $sharediso{$refseq2ucsc{$sh}} = 1; }
		# $sharediso {'uc001aaa'} = 1
		#print Dumper(%sharediso);
		#my %nc_sample_ex;
		my %c_sample_ex;
		
		foreach my $exfile (<*.exon_quantification.txt>)
		{
			#$nc_sample_ex{$exfile} = 0;
			$c_sample_ex{$exfile} = 0;
		}
		#print $str2;
		foreach my $ln (split /\n/, $str2) 
		{
			chomp $ln;
			next if (!$ln);
			my @items = split /\t/, $ln;
			my ($prefix) = split /\./, $items[0];
			#$nc_sample_ex{$items[0]} = $items[2];
			my $isex = "$prefix.rsem.isoforms.normalized_results";
			#print $isex . "\n";
			my $isosum = 0;
			my $part = 0;
			#print Dumper(%nc_sample_ex);
			open (IS, $isex);
			while (my $l = <IS>) {
				chomp $l;
				next if (!$l);
				my @its = split /\t/, $l;
				#print $its[0];
				#print @its . "\n";
				my ($tempp)=split(/\./,$its[0]);
                $its[0]=$tempp;
				#print "$its[0]\n";
				$isosum += $its[1] if (exists $sharediso{$its[0]} or $its[0] eq $refseq2ucsc{$ref});
				#my $tempp=split(/\./,$its[0]);
				#$its[0]=$tempp;
				$part = $its[1] if ($its[0] eq $refseq2ucsc{$ref});
			}
			#print "part:" . $part . "\n";
			#print "sum :" . $isosum . "\n";
			close (IS);
			$c_sample_ex{$items[0]} = scalar(keys %sharediso) == 0 
				? $items[2] 
				: ($isosum == 0 ? 0 : $items[2] * ($part / $isosum));
		}
		#print FW1 $refseq . "\t" . $exon . "\t" . join ("\t", (values %nc_sample_ex)) . "\n";
		print FW2 $refseq . "\t" . $exon . "\t" . join ("\t", (values %c_sample_ex)) . "\n";
	}
}

#close(FW1);
close(FW2);
                
# END #






