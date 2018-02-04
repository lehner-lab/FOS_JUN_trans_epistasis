use List::Util qw(sum min);
use String::Approx qw(amatch);


# which sub-library?
$lib = $ARGV[0];
# threshold at which reads with average quality score are filtered out
$Qthresh = $ARGV[1];


# hash to convert Phred score to numerical values
@phred = split //, '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJ';
@phred{@phred} = (0..$#phred);


# WT sequences
$wtSeq = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA";
$constantWtSeqF = "AACCGGAGGAGGGAGCTG";
$constantWtSeqR = "GAGTTCATCCTGGCAGC";


# length of the different regions of each read
$bLengthF = 11;
$constantLengthF = length($constantWtSeqF);
$bLengthR = 8;
$constantLengthR = length($constantWtSeqR);


# input and output files
$in = "002-Merged_seq/FOSintra_lib" . $lib . "_merged.txt.assembled.fastq";
$out = "003-filtered/Q" . $Qthresh . "_lib" . $lib . ".txt";
open FILE, $in or die "Can't open $in";
open OUT, ">$out" or die "Can't write in $out\n";


# counters for sample stats
$nTotalReads = 0;
$lowQual = 0;
$readsPassingFilter = 0;
$numberErrorsConstantF = 0;
$numberErrorsConstantR = 0;
$SNPlowScore = 0;

%QcountConstantErrorsF = ();
%QcountConstantErrorsR = ();
%QcountConstantCorrectF = ();
%QcountConstantCorrectR = ();

while ($l = <FILE>){
	
	
	# keep only the sequence and the quality score
	for ($i = 2; $i <= 4; $i++){
		$f = <FILE>;
		if ($i == 2){
			$seqF = $f;
			chomp $seqF;
			$nTotalReads++;
		}
		if ($i == 4){
			$qualF = $f;
			chomp $qualF;
		}
	}
	
	
	###########################
	### First filtering steps
	###########################
		
	# get barcodes sequence and quality scores. Do not consider the first base of the barcode, which is most of the time a N (last base in case of the reverse barcode)
	$barF = substr $seqF, 1, $bLengthF - 1;
	$barR = substr $seqF, -$bLengthR, $bLengthR - 1;
	
	# get variable regions sequence and quality scores
	$varF = substr $seqF, $bLengthF + $constantLengthF, 96;
	$QvarF = substr $qualF, $bLengthF + $constantLengthF, 96;
	# convert Phred scores to numerical values
	@QvarF = @phred{(split //, $QvarF)};
	# compute average Phred score for variable region
	$avQF = sum(@QvarF) / 96;
	
	
	# Filter paired reads with low quality or N in barcode or variable region
	if($avQF <= $Qthresh || $barF =~ /N/ || $barR =~ /N/ || $varF =~ /N/){
		$lowQual++;
		next;
	}
	
	
	
	###########################
	### Find mutation positions
	###########################
	
	# Initialize
	@diffs_nt_F = ();
	@diffs_codon_F = ();
	%diffs_codon_F = ();
	
	# Find nt differences in forward read
	@diffs_nt_F = find_diffs($varF,$wtSeq);
	# Identify codon position
	@diffs_codon_F = map {int($_ / 3)} @diffs_nt_F;
	# Remove duplicated codon positions (e.g. when two mutations in the same codon)
	@diffs_codon_F{@diffs_codon_F} = @diffs_codon_F;
	@diffs_codon_F = sort {$a <=> $b} keys %diffs_codon_F;
	
	
	# Flag reads with mutation with high Qscore
	@QmutF = @QvarF[@diffs_nt_F];
	if ((grep {$_ < 30} @QmutF) > 0){
		$SNPlowScore++;
		$flagLowScore = 1;
	}else{
		$flagLowScore = 0;
	}
	
	
	$readsPassingFilter++;
	
	print OUT "$barF$barR\t$varF\t@diffs_codon_F\t@diffs_nt_F\t$flagLowScore\n";
	
	
	
	#####################################
	### Check quality in constant region
	### However, note that paired-end reads do not overlap over the constant region
	#####################################
	
	# get sequences and quality scores
	$constantF = substr $seqF, $bLengthF, $constantLengthF;
	$constantR = substr $seqF, -$bLengthR-$constantLengthR, $constantLengthR;
	$QconstantF = substr $qualF, $bLengthF, $constantLengthF;
	$QconstantR = substr $qualF, -$bLengthR-$constantLengthR, $constantLengthR;
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
close OUT;

$propErrorsConstantF = $numberErrorsConstantF / $constantLengthF / $readsPassingFilter;
$propErrorsConstantR = $numberErrorsConstantR / $constantLengthR / $readsPassingFilter;

print "Filter threshold = $Qthresh
Total number of read = $nTotalReads
No reads filtered because of low quality = $lowQual
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

