$Qthresh = $ARGV[0];


open FILE, "../../000-data/000-genetic_code.txt" or die;
while (<FILE>){
    chomp;
    @a = split /\t/, $_;
    $code{$a[0]} = $a[1];
}
close FILE;

# amino acid sequences
$wt_F = "TDTLQAETDQLEDEKSALQTEIANLLKEKEKL";
# nucleotide sequences
$wt_F2 = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA";


open FILE1, "003-filtered/Q" . $Qthresh . "_lib1.txt" or die;
open FILE2, "003-filtered/Q" . $Qthresh . "_lib2.txt" or die;
open FILE3, "003-filtered/Q" . $Qthresh . "_lib3.txt" or die;
open FILE4, "003-filtered/Q" . $Qthresh . "_lib4.txt" or die;
open FILE5, "003-filtered/Q" . $Qthresh . "_lib5.txt" or die;
open FILE6, "003-filtered/Q" . $Qthresh . "_lib6.txt" or die;


%var = ();
while (<FILE1>){
	chomp;
	$bar = substr $_, 0, 18, ""; # including the tab, to get rid of it. Won't affect UMI counting
	$flag = substr $_, -2, 2, ""; # including the tab
	$var{$_}{1}{$bar} = 1;
}
close FILE1;

while (<FILE2>){
	chomp;
	$bar = substr $_, 0, 18, ""; # including the tab
	$flag = substr $_, -2, 2, ""; # including the tab
	$var{$_}{2}{$bar} = 1;
}
close FILE2;

while (<FILE3>){
	chomp;
	$bar = substr $_, 0, 18, ""; # including the tab
	$flag = substr $_, -2, 2, ""; # including the tab
	$var{$_}{3}{$bar} = 1;
}
close FILE3;

while (<FILE4>){
	chomp;
	$bar = substr $_, 0, 18, ""; # including the tab
	$flag = substr $_, -2, 2, ""; # including the tab
	$var{$_}{4}{$bar} = 1;
}
close FILE4;

while (<FILE5>){
	chomp;
	$bar = substr $_, 0, 18, ""; # including the tab
	$flag = substr $_, -2, 2, ""; # including the tab
	$var{$_}{5}{$bar} = 1;
}
close FILE5;

while (<FILE6>){
	chomp;
	$bar = substr $_, 0, 18, ""; # including the tab
	$flag = substr $_, -2, 2, ""; # including the tab
	$var{$_}{6}{$bar} = 1;
}
close FILE6;



open OUT1, ">../../000-data/001-raw_data/cis_Q" . $Qthresh . ".txt";
print OUT1 "F\tseq\ti1\ti2\ti3\to1\to2\to3\tnt_pos_F\twt_nt_F\tmut_nt_F\taa_pos_F\twt_aa_F\tmut_aa_F\tns_aa_pos\tns_aa_wt\tns_aa_mut\n";
foreach $var (keys %var){
    @counts = ();
    for ($i = 1; $i <= 6; $i++){
        if (exists $var{$var}{$i}){
            $count = keys %{$var{$var}{$i}};
            push @counts, $count;
        }else{
            push @counts, 0;
        }
    }
	
	
	@var = split /\t/, $var;
	# get aa and nt substitution
	if ($var[0] eq $wt_F2){
		$F = "wt";
		$nt_pos_F = NA;
		$wt_nt_F = NA;
		$mut_nt_F = NA;
		$aa_pos_F = NA;
		$wt_aa_F = NA;
		$mut_aa_F = NA;
		$ns_pos = NA;
		$ns_wt_aa = NA;
		$ns_mut_aa = NA;
	}else{
		@aa_pos_F = split / /, $var[1];
		@wt_aa_F = map {substr $wt_F, $_, 1} @aa_pos_F;
		@codon = map {substr $var[0], $_*3, 3} @aa_pos_F;
		@mut_aa_F = @code{@codon};
		@nt_pos_F = split / /, $var[2];
		@mut_nt_F = map {substr $var[0], $_, 1} @nt_pos_F;
		@wt_nt_F = map {substr $wt_F2, $_, 1} @nt_pos_F;
		
		$nt_pos_F = $var[2];
		$mut_nt_F = join(' ',@mut_nt_F);
		$wt_nt_F = join(' ',@wt_nt_F);
		$aa_pos_F = $var[1];
		$wt_aa_F = join(' ',@wt_aa_F);
		$mut_aa_F = join(' ',@mut_aa_F);
		
		# which substitutions are actually non-synonymous
		@ns_aa_subst_pos = ();
		@ns_aa_subst_wt = ();
		@ns_aa_subst_mut = ();
		for($i = 0; $i < @aa_pos_F; $i++){
			if($wt_aa_F[$i] ne $mut_aa_F[$i]){
				push @ns_aa_subst_wt, $wt_aa_F[$i];
				push @ns_aa_subst_mut, $mut_aa_F[$i];
				push @ns_aa_subst_pos, $aa_pos_F[$i];
			}
		}
		$ns_pos = join(' ', @ns_aa_subst_pos);
		$ns_wt_aa = join(' ', @ns_aa_subst_wt);
		$ns_mut_aa = join(' ', @ns_aa_subst_mut);
		
		$aa_seq = translate($var[0]);
		if($aa_seq eq $wt_F){
			$F = "S";
		}else{
			$F = "NS";
		}
	}
	
	$s = $var[0];
    $counts = join ("\t", @counts);
    print OUT1 "$F\t$s\t$counts\t$nt_pos_F\t$wt_nt_F\t$mut_nt_F\t$aa_pos_F\t$wt_aa_F\t$mut_aa_F\t$ns_pos\t$ns_wt_aa\t$ns_mut_aa\n";
}
close OUT1;
print "done\n";



sub translate{
	$s = $_[0];
	$length = length ($s);
	$nb_aa = $length / 3;
	@aa = ();
	for ($i = 0; $i < $nb_aa; $i++){
		$codon = substr $s, 0, 3, '';
		push @aa, $code{$codon};
	}
	$aa_seq = join('',@aa);
	return ($aa_seq);
}