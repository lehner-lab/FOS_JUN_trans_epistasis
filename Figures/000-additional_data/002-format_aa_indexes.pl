open FILE, "002-aaindex1.txt" or die;
open OUT, ">002-index1_formatted.txt" or die;
while (<FILE>){
	chomp;
	$lh = substr $_, 0, 1;
	substr $_, 0, 2, "";
	if ($lh eq  "H"){$id = $_}
	if ($lh eq  "D"){
		$_ =~ s/\"//g; 
		$_ =~ s/\'//g; 
		$desc = $_;
		$_ = <FILE>;
		chomp;
		$lh = substr $_, 0, 1;
		substr $_, 0, 2, "";
		if($lh eq " "){
			$_ =~ s/\"//g; 
			$_ =~ s/\'//g; 
			$desc .= $_;
		}
	}
	
	if ($lh eq "I"){
		
		$_ =~ s/[ ]{1,}/\//g;
		if($_ ne "/A/L/R/K/N/M/D/F/C/P/Q/S/E/T/G/W/H/Y/I/V"){
			print "$_\n";
		}
		
		$nl1 = <FILE>;
		chomp $nl1;
		$nl1 =~ s/[ ]{1,}/;/g;
		@nl1 = split /;/, $nl1;
		@nl1 = @nl1[(1 .. $#nl1)];
		
		$nl2 = <FILE>;
		chomp $nl2;
		$nl2 =~ s/[ ]{1,}/;/g;
		@nl2 = split /;/, $nl2;
		@nl2 = @nl2[(1 .. $#nl2)];
		
		
		@m = (@nl1,@nl2);
		$out = join("\t", @m);
		print OUT "$id\t$desc\t$out\n";
	}
}

close FILE;
close OUT;
