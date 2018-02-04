use List::MoreUtils "uniq";

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
$wt_R = "IARLEEKVKTLKAQNSELASTANMLREQVAQL";
# nucleotide sequences
$wt_F2 = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA";
$wt_R2 = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT";


open FILE1, "002-filtered/Q" . $Qthresh . "_lib1.txt" or die;
open FILE2, "002-filtered/Q" . $Qthresh . "_lib2.txt" or die;
open FILE3, "002-filtered/Q" . $Qthresh . "_lib3.txt" or die;
open FILE4, "002-filtered/Q" . $Qthresh . "_lib4.txt" or die;
open FILE5, "002-filtered/Q" . $Qthresh . "_lib5.txt" or die;
open FILE6, "002-filtered/Q" . $Qthresh . "_lib6.txt" or die;


%var = ();
while (<FILE1>){
	chomp;
	$bar = substr $_, 0, 19, ""; # including the tab, to get rid of it. Won't affect UMI counting
	$flag = substr $_, -2, 2, ""; # including the tab
	$var{$_}{1}{$bar} = 1;
}
close FILE1;

while (<FILE2>){
	chomp;
	$bar = substr $_, 0, 19, ""; # including the tab
	$flag = substr $_, -2, 2, ""; # including the tab
	$var{$_}{2}{$bar} = 1;
}
close FILE2;

while (<FILE3>){
	chomp;
	$bar = substr $_, 0, 19, ""; # including the tab
	$flag = substr $_, -2, 2, ""; # including the tab
	$var{$_}{3}{$bar} = 1;
}
close FILE3;

while (<FILE4>){
	chomp;
	$bar = substr $_, 0, 19, ""; # including the tab
	$flag = substr $_, -2, 2, ""; # including the tab
	$var{$_}{4}{$bar} = 1;
}
close FILE4;

while (<FILE5>){
	chomp;
	$bar = substr $_, 0, 19, ""; # including the tab
	$flag = substr $_, -2, 2, ""; # including the tab
	$var{$_}{5}{$bar} = 1;
}
close FILE5;

while (<FILE6>){
	chomp;
	$bar = substr $_, 0, 19, ""; # including the tab
	$flag = substr $_, -2, 2, ""; # including the tab
	$var{$_}{6}{$bar} = 1;
}
close FILE6;



open OUT1, ">../../000-data/001-raw_data/trans_Q" . $Qthresh . ".txt";
print OUT1 "F\tJ\tseq\ti1\ti2\ti3\to1\to2\to3\tnt_pos_F\twt_nt_F\tmut_nt_F\taa_pos_F\twt_aa_F\tmut_aa_F\tnt_pos_J\twt_nt_J\tmut_nt_J\taa_pos_J\twt_aa_J\tmut_aa_J\n";
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
	}else{
		$wt_aa_F = substr $wt_F, $var[2], 1;
		$codon = substr $var[0], $var[2]*3, 3;
		$mut_aa_F = $code{$codon};
		@nt_pos_F = split / /, $var[4];
		@mut_nt_F = map {substr $var[0], $_, 1} @nt_pos_F;
		@wt_nt_F = map {substr $wt_F2, $_, 1} @nt_pos_F;
		if($wt_aa_F eq $mut_aa_F){
			$F = "S";
		}else{
			$F = "NS";
		}
		$nt_pos_F = $var[4];
		$mut_nt_F = join(' ',@mut_nt_F);
		$wt_nt_F = join(' ',@wt_nt_F);
		$aa_pos_F = $var[2];
	}
	if ($var[1] eq $wt_R2){
		$R = "wt";
		$nt_pos_R = NA;
		$wt_nt_R = NA;
		$mut_nt_R = NA;
		$aa_pos_R = NA;
		$wt_aa_R = NA;
		$mut_aa_R = NA;
	}else{
		$wt_aa_R = substr $wt_R, $var[3], 1;
		$codon = substr $var[1], $var[3]*3, 3;
		$mut_aa_R = $code{$codon};
		@nt_pos_R = split / /, $var[5];
		@mut_nt_R = map {substr $var[1], $_, 1} @nt_pos_R;
		@wt_nt_R = map {substr $wt_R2, $_, 1} @nt_pos_R;
		if($wt_aa_R eq $mut_aa_R){
			$R = "S";
		}else{
			$R = "NS";
		}
		$nt_pos_R = $var[5];
		$mut_nt_R = join(' ',@mut_nt_R);
		$wt_nt_R = join(' ',@wt_nt_R);
		$aa_pos_R = $var[3];
	}
	
	$s = $var[0] . $var[1];
    $counts = join ("\t", @counts);
    print OUT1 "$F\t$R\t$s\t$counts\t$nt_pos_F\t$wt_nt_F\t$mut_nt_F\t$aa_pos_F\t$wt_aa_F\t$mut_aa_F\t$nt_pos_R\t$wt_nt_R\t$mut_nt_R\t$aa_pos_R\t$wt_aa_R\t$mut_aa_R\n";
}
close OUT1;
print "done\n";

