#!usr/bin/perl 
my $infile = $ARGV[0];
open (IN,"$infile")||die "$!";
while(<IN>)
{
	s/\s+$//;
	my $line = $_;
	my $total=( $line=~ s/:/:/g);
	if ($total > 1) 
	{
		print "$_\n";
	}
}

close IN;
