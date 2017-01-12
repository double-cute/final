# This code is to rename the untared files to sort them #
#!/usr/bin/perl

open(FF,"../../../file_manifest.txt");
my %filehash;
while(<FF>)
{
my $temp=$_;
chomp($tmep);
my @a=split(/\t/,$temp);
$filehash{$a[6]}=$a[5];
}
close(FF);

system("rm -r Data/");
system("mkdir Data");
while((my $k,my $v)=each(%filehash))
{
$gene="rsem.isoforms.normalized_results";
$isoform="rsem.genes.normalized_results";
$exon="exon_quantification.txt";
	if(($k=~/$gene/) || ($k=~/$isoform/)||($k=~/$exon/))
	{
		$k=~/edu\..*\.[0-9]+\.?(.*)/;
		#print "$1\n";
		$exoa=$1;
		$v='Data/'.$v;
		if($k=~/$exon/)
		{
			$exoa=~/(.*)\.exon/;
			$exoa=~s/$1//;
			$exoa=~s/\.//;
		}
		$v=$v.".$exoa";
		open(FF,"$k");
		my @a = <FF>;
		open(FW,">$v");
		print FW @a;
		close(FW);
	}
}

# Get the Matched  Data #
system("rm -r Cancer_Matched");
system("rm -r Normal_Matched");
system("mkdir Cancer_Matched");
system("mkdir Normal_Matched");
@allfiles=<Data/*>;
@files=<Data/*-11A-*.rsem.isoforms.normalized_results>;
foreach $file (@files)
{
        #$file=~/TCGA-(.*-)11A(-.*)\.rsem/;

        my @a=split(/rsem/,$file);
        my $normat=$a[0];
        $normat=~/(TCGA.*-)11A-.*R(.*)/;
        my $t1=$1;
        my $t2=$2;
        my $fcount=0;
        foreach my $afile(@allfiles)
        {
                if($afile=~/$t1.*$t2/)
                {
                        $fcount++;
                }
        }
        if($fcount == 6)
        {
                #print "HIT\n";
                foreach my $afile(@allfiles)
                {
                        if($afile=~/$t1.*$t2/)
                        {
                                $afile=~s/Data\///;
				if($afile=~/-01A-/){ system("cp Data/$afile Cancer_Matched/");}
                                if($afile=~/-11A-/){ system("cp Data/$afile Normal_Matched/");}
                        }
                }

                $fcount=0;
        }
}


# Get Cancer Data #
system("rm -r Cancer/");
system("mkdir Cancer");
@files=<Data/*-01A-*rsem.isoforms.normalized_results>;
foreach my $file (@files)
{
	my @a=split(/rsem/,$file);
	$ca[0]=$a[0]."exon_quantification.txt";
	$ca[1]=$a[0]."rsem.genes.normalized_results";
	$ca[2]=$a[0]."rsem.isoforms.normalized_results";
	if((open(FF,$ca[0]))&&(open(FF,$ca[1]))&&(open(FF,$ca[2])))
	{
		foreach $caf (@ca)
		{
			$caf=~s/Data\///;
			$cmd="cp Data/$caf Cancer/";
			system($cmd);
		}
	}
}


# Get Normal Data #
system("rm -r Normal/");
system("mkdir Normal");
@files=<Data/*-11A-*rsem.isoforms.normalized_results>;
foreach my $file (@files)
{
	my @a=split(/rsem/,$file);
        $no[0]=$a[0]."exon_quantification.txt";
        $no[1]=$a[0]."rsem.genes.normalized_results";
        $no[2]=$a[0]."rsem.isoforms.normalized_results";
        if((open(FF,$no[0]))&&(open(FF,$no[1]))&&(open(FF,$no[2])))
        {
                foreach $nof (@no)
                {
                        $nof=~s/Data\///;
			#print "$nof\n";
                        $cmd="cp Data/$nof Normal/";
                        system($cmd);
                }
        }

}
       
