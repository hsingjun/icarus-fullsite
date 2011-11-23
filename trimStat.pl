#!usr/bin/perl
#use warnings;

#Project: JL01
#Function: find where each read is trimmed

#my $tobetrimedfl  = "C:\\hsingjun.Graduate\\datapool\\P3kRawReads";
#my $trimfl = "C:\\hsingjun.Graduate\\datapool\\P3kTrimAdaptor";

if ($#ARGV < 5)
{
	print "usage:perl program raw-fastq adaptor-list qual-file trimA-file trim-common-file trimuniq-files[...]\n";
	exit();
}

#sample usage:
# perl trimStat.pl /home/zhangxin/pro/JL01/data/P3b ./raw adaptor  P3b.qual P3b.qual.trimA P3b.qual.trimA.trimCom JL4/ JL5/ JL6/


my $rawfl = $ARGV[0];
#my $tobetrimedfl = $ARGV[1];#trimCom file
my $uniqadafl  = $ARGV[1];# id-value pair of uniq part of each adaptor
my $qualfl     = $ARGV[2];
my $trimAfl    = $ARGV[3];
my $trimcomfl  = $ARGV[4];


my @library_full_name = split(/\//,$rawfl);
my $library_name = $library_full_name[-1];

$trimcomfl =~ s/\/$//;
my $trimuniqfl =();
my @trimedSet = @ARGV[5..$#ARGV];
my $tobetrimedfl = $trimcomfl;

#split large files into small ones
my $curpath = $ENV{'PWD'};
$curpath =~ s/\/$//;
$curpath .= "\/work";
`mkdir $curpath`;

my @samll_trimcom_file_fullname = split(/\//,$tobetrimedfl);
my $samll_trimcom_file_folder = "$curpath\/$samll_trimcom_file_fullname[-1]"; # folder to put splitted samll trimCom file
`mkdir $samll_trimcom_file_folder`;
`split -l 260000 $tobetrimedfl $samll_trimcom_file_fullname[-1].`;
`mv $samll_trimcom_file_fullname[-1].* $samll_trimcom_file_folder`;
$tobetrimedfl = $samll_trimcom_file_folder;


foreach my $ff ( 0..$#trimedSet)
{
	$ff =~ s/\/$//;
	my @small_trimuniq_file_fullname = split(/\//,$trimedSet[$ff]);
	my $small_trimuniq_file_folder   = "$curpath\/$small_trimuniq_file_fullname[-1]";
	`mkdir $small_trimuniq_file_folder`;
	`split -l 260000 $trimedSet[$ff] $small_trimuniq_file_fullname[-1].`;
	`mv $small_trimuniq_file_fullname[-1].* $small_trimuniq_file_folder`;
	$trimedSet[$ff] = $small_trimuniq_file_folder;
}


my $readsofadaptor = fullAdaptor($rawfl);

my %uniqadaptor = ();
my %notrimmed = ();
getReadId($trimcomfl,\%notrimmed);


#get uniq adaptor for each library
if ($uniqadafl ne '')
{
	open (IN,"$uniqadafl")||die "$!";
	while(<IN>)
	{
		s/\s+$//;
		if (/$library_name/)
		{
			my @ada = split(/\s+/,$_);
			$uniqadaptor{$ada[1]} = $ada[0];
		}
	}
	close IN;
}

my @rawreadsfl = ();
getfiles($tobetrimedfl, \@rawreadsfl);		#get samll raw files
my %howtotrim = ();

$tobetrimedfl =~ s/\/$//;
my @container = split(/\//,$tobetrimedfl);
my $outfile = $container[$#container];
my $readpool = $outfile;
$outfile .= ".trimpos";
#open (OUT,">$outfile")||die "$!";
open (TM,">$readpool.read.pool")||die "$!";#output telling one reads is trimmed by which adaptor(s)
my %readcomposition = ();

foreach my $trimfl (@trimedSet) # for each trimmed LIBRARY
{
	print "\n";
	my @trimedreadsfl = ();
	getfiles($trimfl,\@trimedreadsfl); # get small trimmed files
	
	my $finalreadfl = $trimfl;
	$finalreadfl =~ s/\/$//;
	open (OO,">$finalreadfl.clean.mir.fq")||die "$!";

	#print TM "\t";
	#foreach my $kk (sort keys %uniqadaptor)
	#{
	#	print TM "$kk\t";
	#}
	#print TM "\n";
	
	my %trimedreads = ();
	my %rawreads  = ();
	my @fastq_class = ();

	my $posindicator = 0;
	my $scannedRawReads = 0;
	my $scannedTrimmedReads = 0;
	
	readFastq($trimedreadsfl[0],\%trimedreads); # read one small trimmed file
	print "reading $trimedreadsfl[0] ... finished\n";
	shift @trimedreadsfl;
	

	foreach my $raw(@rawreadsfl) # each small raw file
	{
		#$posindicator = 0;
		open (RW,"$raw")||die "$!";
		while(<RW>)
		{
			s/\s+$//;
			$posindicator ++;
			push @fastq_class, $_;
			if ($posindicator == 4)
			{
				if ($trimedreads{$fastq_class[0]} ne '') # The read is in .trime file
				{
					if ($trimedreads{$fastq_class[0]}[1] ne $fastq_class[1]) # it is trimed
					{
						#$notrimmed{$fastq_class[0]} = 0;
						#$fastq_class[1] =~ m/$trimedreads{$fastq_class[0]}/;
						#if (pos($_) > 0)
						#{
						$_ = $fastq_class[1];
						while(m/$trimedreads{$fastq_class[0]}[1]/g)
						{
							my $cutpos =  pos $_; # find cut position
							my @aread = split(//,$fastq_class[1]);
							my $firsthalf = join('',@aread[0..$cutpos-1]);
							$cutpos ++;
	
							my $secondhalf = join('',@aread[$cutpos-1..$#aread]);
							# find which adaptor is trimmed
							if ($uniqadafl ne '')
							{
								foreach my $aa (keys %uniqadaptor)
								{
									#print TM "$fastq_class[0]\t";
									if ( $secondhalf eq $uniqadaptor{$aa}) # if the trimmed part can match uniq-adaptor
									{
										$readcomposition{"$aa:$secondhalf"} +=1;
										push @{$howtotrim{$fastq_class[0]}}, "$aa:$secondhalf";#if this read is trimmed by adaptor:$aa
										
										$notrimmed{$fastq_class[0]} +=1;
										
										foreach my $rr (@{$trimedreads{$fastq_class[0]}})
										{
											print OO "$rr\n";
										}
										#print TM "$aa\t";
										#last;
									}
									#print TM "\n";
								}
							
							}
							#This statement is right, but not needed now.
							#It outputs the cut position which separating the adaptor from microRNA part.
							#print OUT "$fastq_class[0]\t$firsthalf-$secondhalf\t$trimedreads{$fastq_class[0]}\t$cutpos\n";
						}	
						$scannedTrimmedReads ++;
					}
					else # untrimmed
					{
						#$notrimmed{$fastq_class[0]} +=1; # !!!this is not right to use. It will collect any reads that are trimmed any number of nt, whether it is  adaptor or not
						my $trimp = length($fastq_class[1])+1;
						#This statement is right, but not needed now.
						#It outputs the cut position which separating the adaptor from microRNA part.
						#print OUT "$fastq_class[0]\t$fastq_class[1]\t$trimedreads{$fastq_class[0]}\t$trimp\n";
					}
					
					delete $trimedreads{$fastq_class[0]};
				
					if ( (my $hsize = keys %trimedreads ) == 0)#read next small trimmed file
					{
						if ($#trimedreadsfl ne -1)
						{
							readFastq($trimedreadsfl[0],\%trimedreads);
							print "reading $trimedreadsfl[0] ... finished\n";
							shift @trimedreadsfl;
						}
						else
						{
							print "!!!! trimed reads are all finished\n";
						}
					}
				}
				$posindicator = 0;
				@fastq_class = ();
				$scannedRawReads ++;
				}
			}
		close RW;
			
	}
	
	#close OUT; #$outfile
	close OO;
}

#WildCard filter
my $rawfqc = readCount($rawfl);
my $wildcard_filter = wildcardNum($rawfl);#number of reads containing more than 2 wildcard 'N'
#get number of reads failed in quality filter, adaptor trimming steps (quality filter, polyA filter)
my $qualfilterc  = readCount($qualfl);
my $trimAfilterc = readCount($trimAfl);
my $trimcomc     = readCount($trimcomfl);
my @trimuniqc    = ();
#foreach my $rr(@)
#readCount($trimuniqfl);



open(OUT,">$outfile.stat")||die "$!";
my $lessthan18 = 0;
#print OUT "scannedRawReads\t$scannedRawReads\nscannedTrimmedReads\t$scannedTrimmedReads\n";
print OUT "#reads in raw fastq file:\t$rawfqc\n";
print OUT "#reads containing full length adaptor:\t$readsofadaptor\n";
print OUT "#reads containing more than 2 wildcards:\t$wildcard_filter\n";
$lessthan18 = $rawfqc - $wildcard_filter - $qualfilterc;
print OUT "#reads after QC:\t$qualfilterc (filtered by length < 18: $lessthan18)\n";
$lessthan18 = $qualfilterc - $trimAfilterc;
print OUT "#reads after QC,trim polyA:\t$trimAfilterc (filtered by length < 18: $lessthan18)\n";
$lessthan18 = $trimAfilterc - $trimcomc;
print OUT "#reads after QC,trim polyA, trim adaprot common part:\t$trimcomc (filtered by length < 18: $lessthan18)\n";
#print OUT "#reads after QC,trim polyA, trim adaprot common part,unique part:\t$trimuniqc\n";
#close OUT;

my $notrim = 0;
my $trimonce = 0;
my $trimtwice = 0;
my $trimthree = 0;
my $libc = $#trimedSet +1;

foreach my $tt (sort keys %notrimmed)
{	
	if ($notrimmed{$tt} == 2 )
	{
		$trimtwice ++;
		print "$tt\n";
	}
	elsif ($notrimmed{$tt} == 1)
	{
		$trimonce++;
	}
	elsif ($notrimmed{$tt} eq 0)
	{
		$notrim++;
	}
	elsif( $notrimmed{$tt} == 3)
	{
		$trimthree ++;
	}
	else{}
}

print OUT "$notrim reads not trimmed\n$trimonce reads trimmed once\n$trimtwice reads trimmed twice\n$trimthree reads trimmed three times\n";

#open (OUT,">read_composition.out")||die "$!";
foreach my $aa (keys %readcomposition)
{
	print OUT "$aa\t$readcomposition{$aa}\n";
}
#close OUT;

my %intersect = ();
foreach my $hh (keys %howtotrim)
{
	print TM "$hh\t";
	my $inters = "";
	my $cc = 0;
	foreach my $dd (@{$howtotrim{$hh}})
	{
		$inters .= "$dd-";
		print TM "$dd\t";
		$cc ++;
	}
	if ($cc > 1 )
	{
		$intersect{$inters} ++;
	}
	print TM "\n";
}
close TM;

print OUT "reads that can be trimmed by more than 2 adaptors:\n";
foreach my $ii (sort keys %intersect)
{
	print OUT "$ii\t$intersect{$ii}\n";
}

close OUT;



#--------------------------------sub funtions----------------------------------------------------
sub readFastq
{
	#return hash for each read.
	#4 lines for a read, keep fastq format
	
	my $file = shift;
	my $array = shift;
	my @fastqclass = ();
	my $indicator = 0;
	open (IN,"$file")||die "$!";
	while(<IN>)
	{
		s/\s+$//;
		push @fastqclass, $_;
		$indicator ++;
		if ($indicator == 4)
		{
			$$array{$fastqclass[0]} = [@fastqclass];
			$indicator = 0;
			@fastqclass = ();
		}
	}
	close IN;
}


sub getfiles
{
	my $dir = shift;
	my $files = shift;
	
	my @array = glob("$dir\/*");
	foreach my $item (@array)
	{
		if (-f $item)
		{
			push @$files, $item;
		}
	}
	
}



sub wildcardNum
{
	my $files = shift;
	my $wildcardc = 0;


		open (IN,"$files")||die "$!";
		my $fastqindicator = 0;
		while(<IN>)
		{
			s/\s+$//;
			$fastqindicator ++;
			if (($fastqindicator %4) == 2)
			{
				my $read = $_;
				my $ncount = ($read =~ s/N/N/g);
				if ($ncount > 2)
				{
					$wildcardc++;
				}
				#$fastqindicator = 0;
			}

		}
		close IN;
	
	return $wildcardc ;
}



sub readCount
{
	my $fqfile = shift;#fastq file
	my $readc = 0;
	open (IN,"$fqfile")||die "$!";
	while(<IN>)
	{
		s/\s+$//;
		if (/^@/)
		{
			$readc ++;
		}
	}
	close IN;

	return $readc;

}

sub fullAdaptor
{
	my $fqfile = shift;
	my $adaptorcommon = 'TCGTATGCCGTCTTCTGCTTG';
	my $fastqindicator = 0;
	my $readsc = 0;
	
	open (IN,"$fqfile")||die "$!";
	while(<IN>)
	{
		s/\s+$//;
		$fastqindicator++;
		if (($fastqindicator %4) ==2)
		{
			#my $read = $_;
			if (/$adaptorcommon/)
			{
				$readsc ++;
			}
		}
	}
	close IN;
	return $readsc;
}




sub getReadId
{
	my $infile = shift;
	my $array = shift;

	open (IN,"$infile")||die "$!";
	while(<IN>)
	{
		s/\s+$//;
		if (/^@/)
		{
			$$array{$_} = 0;
		}
	}
	my $size = keys %$array;
	close IN;
	
}
















