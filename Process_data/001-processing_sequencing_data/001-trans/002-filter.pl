# First processing step:
# 1 - matching forward and reverse reads
# 2 - Filter reads based on average quality, the ones containing Ns in variable or barcode regions, the ones with more than one codon mutated on either side of the paired read and the ones ending with AT
# 3 - Flag reads where the mutated base has a quality score < 30
# 4 - then measure sequencing errors in the primer annealing region (done after filtering but ignoring the flag so that it matches the data we will be working with)

use List::Util qw(sum min);


# which sub-library?
$lib = $ARGV[0];
# threshold at which reads with average quality score are filtered out
$Qthresh = $ARGV[1];


# hash to convert Phred score to numerical values
@phred = split //, '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ';
@phred{@phred} = (0..$#phred);


# WT sequences
$wtSeqF = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA";
$wtSeqR = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT";
$constantWtSeqF = "AACCGGAGGAGGGAGCTG";
$constantWtSeqR = "GAAAAAGGAAGCTGGAGAGA";


# length of the different regions of each read
$bLengthF = 11;
$constantLengthF = length($constantWtSeqF);
$bLengthR = 9;
$constantLengthR = length($constantWtSeqR);


# input and output files
$forward_file = "000-sequences/C7GCHANXX_FosJun1_15s014637-1-1_Li_lane5index" . $lib . "_1_sequence.txt";
$reverse_file = "000-sequences/C7GCHANXX_FosJun1_15s014637-1-1_Li_lane5index" . $lib . "_2_sequence.txt";
$out = "002-filtered/Q". $Qthresh . "_lib" . $lib . ".txt";
open FILE, $forward_file or die "Can't open $forward_file\n";
open FILE2, $reverse_file or die "Can't open $reverse_file\n";
open OUT, ">$out" or die "Can't write in $out\n";


# counters for sample stats
$nTotalReads = 0;
$lowQual = 0;
$tooManyMutatedCodons = 0;
$endAT = 0;
$readsPassingFilter = 0;
$numberErrorsConstantF = 0;
$numberErrorsConstantR = 0;
$SNPlowScore = 0;

%QcountConstantErrorsF = ();
%QcountConstantErrorsR = ();
%QcountConstantCorrectF = ();
%QcountConstantCorrectR = ();


while (<FILE>){
	
	
	# get reverse read
	$r = <FILE2>;
	
	# keep only the sequence and the quality score
	for ($i = 2; $i <= 4; $i++){
		$f = <FILE>;
		$r = <FILE2>;
		if ($i == 2){
			$seqF = $f;
			chomp $seqF;
			$seqR = $r;
			chomp $seqR;
			$nTotalReads++;
		}
		if ($i == 4){
			$qualF = $f;
			chomp $qualF;
			$qualR = $r;
			chomp $qualR;
		}
	}
	
	
	
	###########################
	### First filtering step
	###########################
		
	# get barcodes sequence and quality scores. Do not consider the first base of the barcode, which is most of the time a N
	$barF = substr $seqF, 1, $bLengthF - 1;
	$barR = substr $seqR, 1, $bLengthR - 1;
	
	# get variable regions sequence and quality scores
	$varF = substr $seqF, $bLengthF + $constantLengthF, 96;
	$varR = substr $seqR, $bLengthR + $constantLengthR, 96;
	$QvarF = substr $qualF, $bLengthF + $constantLengthF, 96;
	$QvarR = substr $qualR, $bLengthR + $constantLengthR, 96;
	# convert Phred scores to numerical values
	@QvarF = @phred{(split //, $QvarF)};
	@QvarR = @phred{(split //, $QvarR)};
	# compute average Phred score for variable region
	$avQF = sum(@QvarF) / 96;
	$avQR = sum(@QvarR) / 96;
	
	
	# Filter paired reads with low quality or N in barcode or variable region
	if($avQF <= $Qthresh || $avQR <= $Qthresh || $barF =~ /N/ || $barR =~ /N/ || $varF =~ /N/ || $varR =~ /N/){
		$lowQual++;
		next;
	}
	
	
	
	###########################
	### Find mutation positions
	###########################
	
	# Initialize
	@diffs_nt_F = (); @diffs_nt_R = ();
	@diffs_codon_F = (); @diffs_codon_R = ();
	%diffs_codon_F = (); %diffs_codon_R = ();
	
	# Find nt differences in forward read
	@diffs_nt_F = find_diffs($varF,$wtSeqF);
	# Identify codon position
	@diffs_codon_F = map {int($_ / 3)} @diffs_nt_F;
	# Remove duplicated codon positions (e.g. when two mutations in the same codon)
	@diffs_codon_F{@diffs_codon_F} = @diffs_codon_F;
	@diffs_codon_F = keys %diffs_codon_F;
	# same for reverse read
	@diffs_nt_R = find_diffs($varR,$wtSeqR);
	@diffs_codon_R = map {int($_ / 3)} @diffs_nt_R;
	@diffs_codon_R{@diffs_codon_R} = @diffs_codon_R;
	@diffs_codon_R = keys %diffs_codon_R;
	
	###########################
	### Second filtering steps
	###########################
	
	# Filter out reads with more than one mutated codon in bait or prey
	if (@diffs_codon_F > 1 || @diffs_codon_R > 1){
		$tooManyMutatedCodons++;
		next;
	}
	
	# Filter out mutated codon that end by A or T (not designed)
	if (@diffs_codon_F == 1){
		$endF = substr $varF, $diffs_codon_F[0]*3+2, 1;
		if ($endF eq "A" || $endF eq "T"){
			$endAT++;
			next;
		}
	}
	if (@diffs_codon_R == 1){
		$endR = substr $varR, $diffs_codon_R[0]*3+2, 1;
		if ($endR eq "A" || $endR eq "T"){
			$endAT++;
			next;
		}
	}
	
	# Flag reads with mutation with high Qscore
	@QmutF = @QvarF[@diffs_nt_F];
	@QmutR = @QvarR[@diffs_nt_R];
	if ((grep {$_ < 30} @QmutF) > 0 || (grep {$_ < 30} @QmutR) > 0){
		$SNPlowScore++;
		$flagLowScore = 1;
	}else{
		$flagLowScore = 0;
	}
	
	
	$readsPassingFilter++;
	
	print OUT "$barF$barR\t$varF\t$varR\t@diffs_codon_F\t@diffs_codon_R\t@diffs_nt_F\t@diffs_nt_R\t$flagLowScore\n";
	
	
	
	#####################################
	### Check quality in constant region
	#####################################
	
	# get sequences and quality scores
	$constantF = substr $seqF, $bLengthF, $constantLengthF;
	$constantR = substr $seqR, $bLengthR, $constantLengthR;
	$QconstantF = substr $qualF, $bLengthF, $constantLengthF;
	$QconstantR = substr $qualR, $bLengthR, $constantLengthR;
	# convert Phred scores into numerical values
	@QconstantF = @phred{(split //, $QconstantF)};
	@QconstantR = @phred{(split //, $QconstantR)};
	# make a hash with position in constant region as keys and quality score as values
	%QconstantF = ();
	%QconstantR = ();
	@QconstantF{( 0 .. $constantLengthF-1)} = @QconstantF;
	@QconstantR{( 0 .. $constantLengthR-1)} = @QconstantR;
	
	# number and positions of errors
	@diffsConstantF = find_diffs($constantF,$constantWtSeqF);
	@diffsConstantR = find_diffs($constantR,$constantWtSeqR);
	$numberErrorsConstantF += scalar @diffsConstantF;
	$numberErrorsConstantR += scalar @diffsConstantR;
	
	# quality score of errors
	@constantErrorsF = delete @QconstantF{@diffsConstantF};
	@constantErrorsR = delete @QconstantR{@diffsConstantR};
	# quality score of correct positions
	@constantCorrectF = values %QconstantF;
	@constantCorrectR = values %QconstantR;
	
	# count number of bases with this quality score
	$QcountConstantErrorsF{$_}++ foreach (@constantErrorsF);
	$QcountConstantErrorsR{$_}++ foreach (@constantErrorsR);
	$QcountConstantCorrectF{$_}++ foreach (@constantCorrectF);
	$QcountConstantCorrectR{$_}++ foreach (@constantCorrectR);
	
}
close FILE;
close FILE2;
close OUT;


$propErrorsConstantF = $numberErrorsConstantF / $constantLengthF / $readsPassingFilter;
$propErrorsConstantR = $numberErrorsConstantR / $constantLengthR / $readsPassingFilter;

print "Filter threshold = $Qthresh
Total number of read = $nTotalReads
No reads filtered because of low quality = $lowQual
No reads filtered because too many codons mutated = $tooManyMutatedCodons
No reads filtered because codon ends with AT = $endAT
No reads passing filter = $readsPassingFilter
No reads flagged because mutations has low Qscore = $SNPlowScore
Proportion of error per base in constant region (forward read) = $propErrorsConstantF
Proportion of error per base in constant region (reverse read) = $propErrorsConstantR
\n";


print "Distribution of quality scores for errors in constant region (forward read)\n";
print "$_:\t$QcountConstantErrorsF{$_}\n" foreach (sort {$a <=> $b} keys %QcountConstantErrorsF);

print "\nDistribution of quality scores for errors in constant region (reverse read)\n";
print "$_:\t$QcountConstantErrorsR{$_}\n" foreach (sort {$a <=> $b} keys %QcountConstantErrorsR);

print "\nDistribution of quality scores for correct positions in constant region (forward read)\n";
print "$_:\t$QcountConstantCorrectF{$_}\n" foreach (sort {$a <=> $b} keys %QcountConstantCorrectF);

print "\nDistribution of quality scores for correct positions in constant region (reverse read)\n";
print "$_:\t$QcountConstantCorrectR{$_}\n" foreach (sort {$a <=> $b} keys %QcountConstantCorrectR);

print "\n\n\n\n";




use Inline C => << 'EOC';
void find_diffs(char* x, char* y) {										   
  int i;																	  
  Inline_Stack_Vars;														  
  Inline_Stack_Reset;														 
  for(i=0; x[i] && y[i]; ++i) {											   
	if(x[i] != y[i]) {														
	  Inline_Stack_Push(sv_2mortal(newSViv(i)));							  
	}																		 
  }																		   
  Inline_Stack_Done;														  
}																			 
EOC

