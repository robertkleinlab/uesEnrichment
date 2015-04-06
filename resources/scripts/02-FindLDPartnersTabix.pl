#!/usr/bin/perl

######################################################################################
#
# This script will make find the LD SNPs for every SNP in each file in the directory. 
# There should only be .bed files in the BEDfiles directory
#
######################################################################################

$dir = "$ARGV[0]";
@files = <BEDfiles/*.bed>;	#Holds all the files in the directory

foreach $file (@files) {
	open (VARIN, "$file") || die(print"$file failed.\n");
	@f = split(/\//, $file);
	open (LDBED, ">ldBedFiles/ld.$f[1]");
	while (<VARIN>) {
		my %ldHash=();
		chomp;
		@bedLine = split(/\s+/);
		$chromo = $bedLine[0];
	   	$p1 = $bedLine[1];
		$p2 = $bedLine[2];
		$id = $bedLine[3];
		$tabixCommand = "$chromo:"."$p1"."-"."$p1";
		
		print LDBED "$chromo\t$p1\t$p2\t$id\t.\t$id\n";
		$tabixIndex = qx{tabix $dir/resources/snpData/snpIndex.1000gvcf.$chromo.allInfo.bgz $tabixCommand};
		
		$ldHash{$id} = 1;
		
		$tabixLine = qx{tabix $dir/resources/eurLdTabix/EURonly-LD1000g.$chromo.v5.20130502.annotated.final.bgz $tabixCommand};
		@tabixArray = split(/\n/, $tabixLine); 
		foreach (@tabixArray) {
			chomp;
			$ldLine = $_;
			@s = split(/\s+/, $ldLine);
			$id1 = $s[2];
			$pos1 = $s[1];
			$pos11 = $pos1+1;
			$id2 = $s[4];
			$pos2 = $s[3];
			$pos22 = $pos2+1;
			$t1 = "$chromo:"."$pos1"."-"."$pos1";
			$t2 = "$chromo:"."$pos2"."-"."$pos2";
			if (($id1 eq $id) || ($id2 eq $id)) {
				if ( !exists $ldHash{$id1} ) {
					$ldHash{$id1} = 1;
					$t1line = qx{tabix $dir/resources/snpData/snpIndex.1000gvcf.$chromo.allInfo.bgz $t1};
					print LDBED "$chromo\t$pos1\t$pos11\t$id1\t.\t$id\n";
				}
				elsif ( !exists $ldHash{$id2} ) {
					$ldHash{$id2} = 1;
					$t2line = qx{tabix $dir/resources/snpData/snpIndex.1000gvcf.$chromo.allInfo.bgz $t2};
					print LDBED "$chromo\t$pos2\t$pos22\t$id2\t.\t$id\n";
				}
			}
		}
		untie %ldHash;
	}
	close VARIN;
	close LDBED;
}

