#! /bin/perl

# Created 2014/12/08
# This is executed by a the master bash script. 
# It employs tabix -indexed files for speed.
# Use the followign command to return tabix info.
# $line = qx{tabix resources/snpData/snpIndex.1000gvcf.chr1.allInfo.bgz chr1:10505-10505};

my %SNPs; #hash holds the snps that we're analyzing.
my %chipSNPs; # hash holds the chip SNPs.

print "Base file name: $ARGV[0]\n";
$argstring = "$ARGV[0]";
@filename = split ("\\.", $argstring);
$fileTitle=$filename[0];
$dir = "$ARGV[3]";

# Puts original SNPs into a hash so we can sort by going over the list once.
print "\nGenerating a hash of the SNPs of interest.\n";
open (USERIN, "$ARGV[0]") || die();
while (<USERIN>) {
	@s = split (/\s+/);
    $snp = $s[0];
    if ( !exists $SNPs{$snp} ) {
    	$SNPs{$snp} = "$snp";
    }
} 
close USERIN; 


# Load the chip SNPs into memory.
open (CHIPIN, "$ARGV[1]") || die(print "Could not open $ARGV[1]. Script terminating.");
while(<CHIPIN>) {
	chomp;
	$snp=$_;
	if ( !exists $chipSNPs{$snp} ) {
   		$chipSNPs{$snp} = "$snp";
   	}
}

print "Generating bed file for original SNPs...";
open INALL, "gunzip -c $dir/resources/snpData/snpIndex.1000gvcf.allInfo.dat.gz |" || die(print "Failed opening $dir/resources/snpData/snpIndex.1000gvcf.allInfo.dat.gz\n");
open (OUT, ">BEDfiles/$filename[0].original.bed") || die();


# creates a .bed file of the original SNPs and enter them into a hash. SNPs at the HLA region are excluded.
# HLA region is defined as on chromosome 6: 29570005 and ending at chr6 33377658 for build37/hg19.
# Also, this loads the proper CHIP SNPs into memory. It filters out SNPs at HLA and w/ a MAF < 0.05. 
my %BEDHASH;
open (SNPDATA, ">dataFiles/$filename[0].snpData.dat");

my @tssArray=();
my @ldArray=();
my @snpArray=();
my @chipArray=();

while (<INALL>) {
	chomp;
	$line = $_;
	@s = split (/\s+/, $line);
	$chromo = $s[0];
	$pos = $s[1];
	$pos1 = $pos+1;
   	$id = $s[2];
   	$maf = $s[6];
	if ( exists $SNPs{$id} ) {
    	if (($chromo eq "chr6") && ($pos >= 29570005) && ($pos <= 33377658)) {
    		$excluded=1;
    	}
		else {
    		print OUT "$chromo\t$pos\t$pos1\t$id\t0\t.\n";
    		print SNPDATA "$line\n";
    		push (@snpArray, $line);
    		push (@ldArray, $s[7]);	
			push (@tssArray, $s[8]);
		}
	}
	if ( exists $chipSNPs{$id} ) {
		if ((($chromo eq "chr6") && ($pos >= 29570005) && ($pos <= 33377658)) || ($maf < 0.05)) {
    		$excluded=1;
    	}
    	else {
    		push (@chipArray, $line);
    	}
	}
}
close INALL;
close OUT;
untie %SNPs;
untie %chipSNPs;
print "...Done.\n";

@sortedTSS = sort { $a <=> $b } @tssArray;
@sortedLD = sort { $a <=> $b } @ldArray;

#sorts TSS distances to get quartiles
$tssNum = @tssArray;
$ldNum = @ldArray;
#note, these are casted to INT in case that the total number of numbers is not a multiple of 4.
$tssq1 = $sortedTSS[int($tssNum*(1/4))-1];
$tssq2 = $sortedTSS[int($tssNum*(2/4))-1];
$tssq3 = $sortedTSS[int($tssNum*(3/4))-1];
$tssq4 = $sortedTSS[int($tssNum*(4/4))-1];

$ldq1 = $sortedLD[int($ldNum*(1/4))-1];
$ldq2 = $sortedLD[int($ldNum*(2/4))-1];
$ldq3 = $sortedLD[int($ldNum*(3/4))-1];
$ldq4 = $sortedLD[int($ldNum*(4/4))-1];

print "TSS Quantiles:\n";
print "Q1: $tssq1\tQ2: $tssq2\tQ3: $tssq3\tQ4: $tssq4\n\n";
print "LD Quantiles:\n";
print "Q1: $ldq1\tQ2: $ldq2\tQ3: $ldq3\tQ4: $ldq4\n\n";

# bin original input SNPs.
my $bin11=0;
my $bin12=0;
my $bin13=0;
my $bin14=0;
my $bin21=0;
my $bin22=0;
my $bin23=0;
my $bin24=0;
my $bin31=0;
my $bin32=0;
my $bin33=0;
my $bin34=0;
my $bin41=0;
my $bin42=0;
my $bin43=0;
my $bin44=0;

foreach(@snpArray) {
	@b = split(/\s+/);
	$tssDistance = $b[8];
	$ldNum = $b[7];
	if (($tssDistance >= 0 ) && ($tssDistance <= $tssq1)) {
		if ($ldNum <= $ldq1) {
			$bin11++;
		}
		elsif (($ldNum > $ldq1) && ($ldNum <= $ldq2)) {
			$bin12++;
		}
		elsif (($ldNum > $ldq2) && ($ldNum <= $ldq3)) {
			$bin13++;
		}
		elsif ($ldNum > $ldq3) {
			$bin14++;
		}
	}
	elsif (($tssDistance > $tssq1 ) && ($tssDistance <= $tssq2)) {
		if ($ldNum <= $ldq1) {
			$bin21++;
		}
		elsif (($ldNum > $ldq1) && ($ldNum <= $ldq2)) {
			$bin22++;
		}
		elsif (($ldNum > $ldq2) && ($ldNum <= $ldq3)) {
			$bin23++;
		}
		elsif ($ldNum > $ldq3){
			$bin24++;
		}
	}
	elsif (($tssDistance > $tssq2 ) && ($tssDistance <= $tssq3)) {
		if ($ldNum <= $ldq1) {
			$bin31++;
		}
		elsif (($ldNum > $ldq1) && ($ldNum <= $ldq2)) {
			$bin32++;
		}
		elsif (($ldNum > $ldq2) && ($ldNum <= $ldq3)) {
			$bin33++;
		}
		elsif ($ldNum > $ldq3){
			$bin34++;
		}
	}
	elsif ($tssDistance > $tssq3 ) {
		if ($ldNum <= $ldq1) {
			$bin41++;
		}
		elsif (($ldNum > $ldq1) && ($ldNum <= $ldq2)) {
			$bin42++;
		}
		elsif (($ldNum > $ldq2) && ($ldNum <= $ldq3)) {
			$bin43++;
		}
		elsif ($ldNum > $ldq3) {
			$bin44++;
		}
	}
}

$b11 = "bin11 $bin11";
$b12 = "bin12 $bin12";
$b13 = "bin13 $bin13";
$b14 = "bin14 $bin14";
$b21 = "bin21 $bin21";
$b22 = "bin22 $bin22";
$b23 = "bin23 $bin23";
$b24 = "bin24 $bin24";
$b31 = "bin31 $bin31";
$b32 = "bin32 $bin32";
$b33 = "bin33 $bin33";
$b34 = "bin34 $bin34";
$b41 = "bin41 $bin41";
$b42 = "bin42 $bin42";
$b43 = "bin43 $bin43";
$b44 = "bin44 $bin44";

@binArray = ($b11, $b12, $b13, $b14, $b21, $b22, $b23, $b24, $b31, $b32, $b33, $b34, $b41, $b42, $b43, $b44);

# bin all chip SNPs.
@bin11=();
@bin12=();
@bin13=();
@bin14=();
@bin21=();
@bin22=();
@bin23=();
@bin24=();
@bin31=();
@bin32=();
@bin33=();
@bin34=();
@bin41=();
@bin42=();
@bin43=();
@bin44=();

$notAdded1 = 0;
$notAdded2 = 0;
$notAdded3 = 0;
$notAdded4 = 0;
$notAdded = 0;

foreach(@chipArray) {
	$line = $_;
	@s = split(/\s+/);
	$chromo = $s[0];
	$p1 = $s[1];
	$p2 = $p1+1;
	$id = $s[2];
	$tssDistance = $s[8];
	$ldNum = $s[7];
	$bedLine = "$line";
	if (($tssDistance >= 0 ) && ($tssDistance <= $tssq1)) {
		if ($ldNum <= $ldq1) {
			push (@bin11, $bedLine);
		}
		elsif (($ldNum > $ldq1) && ($ldNum <= $ldq2)) {
			push (@bin12, $bedLine);
		}
		elsif (($ldNum > $ldq2) && ($ldNum <= $ldq3)) {
			push (@bin13, $bedLine);
		}
		elsif ($ldNum > $ldq3) {
			push (@bin14, $bedLine);
		}
		else {
			$notadded1++;
		}
		
	}
	elsif (($tssDistance > $tssq1 ) && ($tssDistance <= $tssq2)) {
		if ($ldNum <= $ldq1) {
			push (@bin21, $bedLine);
		}
		elsif (($ldNum > $ldq1) && ($ldNum <= $ldq2)) {
			push (@bin22, $bedLine);
		}
		elsif (($ldNum > $ldq2) && ($ldNum <= $ldq3)) {
			push (@bin23, $bedLine);
		}
		elsif ($ldNum > $ldq3){
			push (@bin24, $bedLine);
		}
		else {
			$notadded2++;
		}
	}
	elsif (($tssDistance > $tssq2 ) && ($tssDistance <= $tssq3)) {
		if ($ldNum <= $ldq1) {
			push (@bin31, $bedLine);
		}
		elsif (($ldNum > $ldq1) && ($ldNum <= $ldq2)) {
			push (@bin32, $bedLine);
		}
		elsif (($ldNum > $ldq2) && ($ldNum <= $ldq3)) {
			push (@bin33, $bedLine);
		}
		elsif ($ldNum > $ldq3){
			push (@bin34, $bedLine);
		}
		else {
			$notadded3++;
		}
	}
	elsif ($tssDistance > $tssq3 ) {
		if ($ldNum <= $ldq1) {
			push (@bin41, $bedLine);
		}
		elsif (($ldNum > $ldq1) && ($ldNum <= $ldq2)) {
			push (@bin42, $bedLine);
		}
		elsif (($ldNum > $ldq2) && ($ldNum <= $ldq3)) {
			push (@bin43, $bedLine);
			$binLineNum = @bin43;
		}
		elsif ($ldNum > $ldq3) {
			push (@bin44, $bedLine);
		}
		else {
			$notadded4++;
		}
	}
	else {
		$notAdded++;
	}
}

$num11 = @bin11;
$num12 = @bin12;
$num13 = @bin13;
$num14 = @bin14;
$num21 = @bin21;
$num22 = @bin22;
$num23 = @bin23;
$num24 = @bin24;
$num31 = @bin31;
$num32 = @bin32;
$num33 = @bin33;
$num34 = @bin34;
$num41 = @bin41;
$num42 = @bin42;
$num43 = @bin43;
$num44 = @bin44;

print "\nOriginal SNP Distritutions:\n\t[,1]\t[,2]\t[,3]\t[,4]\n[1,]\t$bin11\t$bin12\t$bin13\t$bin14\n[2,]\t$bin21\t$bin22\t$bin23\t$bin24\n[3,]\t$bin31\t$bin32\t$bin33\t$bin34\n[4,]\t$bin41\t$bin42\t$bin43\t$bin44\n";
print "\nChipSNP Bin Distrututions:\n\t[,1]\t[,2]\t[,3]\t[,4]\n[1,]\t$num11\t$num12\t$num13\t$num14\n[2,]\t$num21\t$num22\t$num23\t$num24\n[3,]\t$num31\t$num32\t$num33\t$num34\n[4,]\t$num41\t$num42\t$num43\t$num44\n";

# generate random files
$fileCounter = 1;
$randFiles = $ARGV[2];
while ($fileCounter <= $randFiles) {
	open (OUTFILE, ">BEDfiles/$fileTitle"."-randomSNPs-"."$fileCounter.bed");
	foreach (@binArray){
		@b = split(/\s+/);
		$bin = $b[0];
		$num = $b[1];
		$randSNP = 0;
		while ($randSNP < $num) {
			@tempArray2 = @{$bin};
			$numInArray2 = @tempArray2;
			$randomElement = $tempArray2[rand @tempArray2];
			@l = split(/\s+/, $randomElement);
			$chromo = $l[0];
			$p1 = $l[1];
			$p2 = $p1+1;
			$id = $l[2];
			print OUTFILE "$chromo\t$p1\t$p2\t$id\t0\t.\n";
			$randSNP++;
		}
	}
	$fileCounter++;
	close OUTFILE;
}

$fileCounter--;
print "Generated $fileCounter random files.\n";

print "Done.\n\n"
