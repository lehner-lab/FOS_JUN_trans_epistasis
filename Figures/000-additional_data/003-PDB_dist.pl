$pdb_file = "003-1fos.pdb";
$chainA = "E";
$chainB = "F";

%pdb = ();
# contacts in trans
open FILE, $pdb_file or die "Can't open $pdb_file\n";
while (<FILE>){
	
	# for coordinates lines
	if (substr($_, 0, 4) eq "ATOM"){
		# if corresponds to the chain of interest
		$c = substr($_, 21, 1);
		if ($c eq $chainA || $c eq $chainB){
			$atom = substr($_, 13, 3);		$atom =~ s/\s//g;
			# if it is not a side chain atom, skip
			if($atom eq "N" || $atom eq "C" || $atom eq "CA" || $atom eq "O"){next;}
			# get position in chain
			$pos = substr($_, 22, 4);		$pos =~ s/\s//g;
			# get residue
			$aa = substr($_, 17, 3);
			$pos .= "\t" . "$aa";
			# get coordinates
			$x = substr($_, 30, 8);		$x =~ s/\s//g;
			$y = substr($_, 38, 8);		$y =~ s/\s//g;
			$z = substr($_, 46, 8);		$z =~ s/\s//g;
			# store in hash
			if ($atom ne ""){
				$pdb{$c}{$pos}{$atom} = [$x,$y,$z];
			}
		}
	}
}
close FILE;

# for each position pair, compute the shortest distance between all pairs of atom at the two positions
open OUT, ">003-distances_trans.txt";
@chain = sort keys %pdb;
foreach $aaA (keys %{$pdb{$chain[0]}}){
	@atomA = keys %{$pdb{$chain[0]}{$aaA}};
	foreach $aaB (keys $pdb{$chain[1]}){
		@atomB = keys %{$pdb{$chain[1]}{$aaB}};
		$closest = 100000;
		foreach $atomA (@atomA){
			@posA = @{$pdb{$chain[0]}{$aaA}{$atomA}};
			foreach $atomB (@atomB){
				@posB = @{$pdb{$chain[1]}{$aaB}{$atomB}};
				$distx = $posA[0] - $posB[0];
				$disty = $posA[1] - $posB[1];
				$distz = $posA[2] - $posB[2];
				$dist = sqrt(($distx*$distx)+($disty*$disty)+($distz*$distz));
				if ($dist < $closest){
					$closest = $dist;
				}
			}
		}
		print OUT "$aaA\t$aaB\t$closest\n";
	}
}
close OUT;






# same but in cis
%pdb = ();
open FILE, $pdb_file or die "Can't open $pdb_file\n";
while (<FILE>){
	if (substr($_, 0, 4) eq "ATOM"){
		$c = substr($_, 21, 1);
		if ($c eq $chainA){
			$atom = substr($_, 13, 3);		$atom =~ s/\s//g;
			if($atom eq "N" || $atom eq "C" || $atom eq "CA" || $atom eq "O"){next;}
			$pos = substr($_, 22, 4);		$pos =~ s/\s//g;
			$aa = substr($_, 17, 3);
			$pos .= "\t" . "$aa";
			$x = substr($_, 30, 8);		$x =~ s/\s//g;
			$y = substr($_, 38, 8);		$y =~ s/\s//g;
			$z = substr($_, 46, 8);		$z =~ s/\s//g;
			if ($atom ne ""){
				$pdb{$pos}{$atom} = [$x,$y,$z];
			}
		}
	}
}
close FILE;





open OUT, ">003-distances_cis.txt";
foreach $aaA (keys %pdb){
	@atomA = keys %{$pdb{$aaA}};
	foreach $aaB (keys %pdb){
		@atomB = keys %{$pdb{$aaB}};
		$closest = 100000;
		foreach $atomA (@atomA){
			@posA = @{$pdb{$aaA}{$atomA}};
			foreach $atomB (@atomB){
				@posB = @{$pdb{$aaB}{$atomB}};
				$distx = $posA[0] - $posB[0];
				$disty = $posA[1] - $posB[1];
				$distz = $posA[2] - $posB[2];
				$dist = sqrt(($distx*$distx)+($disty*$disty)+($distz*$distz));
				if ($dist < $closest){
					$closest = $dist;
				}
			}
		}
        print OUT "$aaA\t$aaB\t$closest\n";
	}
}
close OUT;
