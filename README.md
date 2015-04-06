# uesEnrichment
- James Hayes
- PhD Candidate, Weill Cornell Graduate School of Medical Sciences, New York, NY
- Lab PI: Robert Klein, PhD, at Icahn School of Medicine at Mount Sinai, New York, NY
- Updated 2015/04/06

## ues (Uncovering Enrichment through Simulation) Algorithm Development
The UES (Uncovering Enrichment through Simulation) algorithm was written to help interpret results from genome-wide association studies (GWAS) using publicly available datasets. Our manuscript is currently under review and is available on the [bioRxvi](http://biorxiv.org) preprint server. The published version of the algorithm can be found at the [Klein lab website](http://research.mssm.edu/kleinlab/ues/).

## Executing the algorithm
UES is comprised of scripts written in bash (3.2.53(1)-release) and PERL  (v5.18.2) scripts. It relies on the following programs:
    • tabix  (v0.2.5)
    • bedtools (v2.20.1)
These programs need to be installed and added to your path in order for UES to operate properly.

UES uses pre-computed SNP that is included in the download for this package.  These file contain the following data: positional data, alleles, MAF, number of LD partners, and distance to TSS.  There is also a separate set of files containing the LD data for all the SNPs. All of these files are compressed, but they are still a large amount of data. 

The most recent version of the UES pipeline can be downloaded: from http://research.mssm.edu/kleinlab/software/ues

UES has 2 main bash scripts: making the matched-random SNPs (ues), and performing the intersection and calculating the enrichment score (uesSingleIntersect & uesBatchIntersect).  

## Running UES: 
Input: The user is required to provide a file that contains the list of SNPs for which they want to test enrichment in a column. No other data is required. UES will pull all of the required positional/maf/ld data on its own. The list of lymphoma SNPs  has been provided. When providing your own SNPs, do not put a header in your file. An example can be found here: ues/testData/lymphoma/lymphoma.snps.txt. These example SNPs have been previously reported and were manually selected from the NHGRI GWAS catalogue.

Options: UES requires 3 types of options: the original input SNPs, a chip code, and the number of random files to generate.

    -h             This flag will open up this manual document.
	-i <FILE>      File containing the column of original SNPs
	-n <integer>   Number of random files to generate

	Chip Codes (Capital Lettes)
     -A          Affymetrix Snp6 Chip (found in 1000Genomes)
     -I          Illumina Omni1-Quad and Human-1M (found in 1000Genomes)
     -C          Combined Affymetrix Snp6, Illumina Omni1-Quad, and Illumina Human-1M (found in 1000Genomes)
     -T          1000 Genomes
     -H          HapMapPhase3
     -U <FILE>   User provided list. As with the input SNPs, they should be in a column, with one SNP per line and no header.

Output:
UES creates the following folders within the directory that you ran the analysis: 
     BEDfiles/   Contains the bed (6 column format: chromosome, start position, stop position, rsid, score, strand) for the original input SNPs and matched random SNPs.
     dataFiles/  Contains the files that have the SNP data
     ldBedFiles/ Contains the bed (6 column format) for the original and matched SNP lists, along with all of the LD partners. Note, column 6 contains the original SNP for which the current SNP is an LD partner.
     analysis/   This folder is initially empty. The intersection scripts will populate this folder with the analysis results. 

The names of the generated output files are based on the name of the original input list of SNPs; it takes the string of characters before the first “.” in your input file. In this example glioma.snps.txt will produce output files that begin with glioma.

Executing the script: The script should be run from the folder where you have the proper input file. The output will be generated within that initial folder.

Example Command:
The SNP list for Hayes et al. is provided. In order to proceed, move to the following directory: ues/test-lymphoma/
The ues bash script is executable. Run the following command:
../ues -i lymphomaCllSnps.txt –C -n 10000

Note: this may take a while to run. For testing purposes, reduce -n to increase speed.

## UES Intersection Enrichment:
Two enrichment scripts are provided; one will perform the enrichment analysis on a single genomic interval track file, while the other will perform them on all of the genomic interval track files in a folder. BedTools requires these tracks to be in tab-delineated, 6-column bed format. Space delineated or files containing >6 columns will cause BedTools to fail.  This will work with either compressed or uncompressed files. Please ensure the files are uncompressed. The DNase tracks used in the published analysis can be found at: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseUniform/

These scripts must be run from the main folder of the ues run, containing the analysis/ and ldBedFiles/ folders.

### uesSingleIntersect:

Options:
      -f <FILE>  Full path to the genomic interval track file
      -o <NAME>  Output name

Generated Output:
      <NAME>.output             Folder containing the intersection files for the collection of random SNPs and the original SNPs.
      summaryStats.<NAME>.txt   This file contains the enrichment calculation. There are 5 columns in the file.
      Track:                    Name of the track against which the intersection was performed.
      OrigLoci                  The number of loci from the original input that co-localized with genomic marks.
      RandLoci>=Orig            The number of occurrences where a random SNP set had >= the number of loci that co-localized to a mark than the original input SNPs.
      RandAvg                   The average number of co-localizations of the random SNP sets.
      pValue                    Calculated pValue.

Example Command:  
From the same folder as before (ues/test-lymphoma) before executing the example command.
NOTE: Please provide the full path to the output file.  Using "../" may cause BedTools to fail.

	../uesSingleIntersect -f [DNaseFolder]wgEncodeAwgDnaseUwdukeGm12878UniPk.narrowPeak.gz -o lymphomaTest-singleIntersect

Output:
	Track     OrigLoci     RandLoci>=Orig     RandAvg     pValue
	wgEncodeAwgDnaseUwdukeGm12878UniPk.narrowPeak.gz     16     0     4.521     <0.0001


### uesBatchIntersect:

Options:
      -f <FILE>  Full path to the folder containing the genomic interval track files
      -o <NAME>  Output name

Generated Output:
The output is essentially identical to what was previously described above.  The resultant summaryStats file is sorted by pValue in ascending order.

Example Command:
From the same folder as before (ues/test-lymphoma) before executing the example command.
NOTE: Please provide the full path to the output file.  Using "../" may cause BedTools to fail.

	../uesBatchIntersect –f [DNaseFolder] -o lymphomaTest-batchIntersect

## UES Minerva Scripts for parallel computing. 
There are some scripts that are written specifically for use on the [Minerva cluster](hpc.mssm.edu) at MSSM. They allow for distribution over the cluster, billing to the proper project account, and submitting to a user-specified queue. You can modify these scripts for use on your own cluster.

### uesMinerva
This is virtually the same script as the main UES script. The only difference is that it has been written to distribue across Minerva to increase speed. There are 2 additional flags required beyond the flags for the normal algorithm run:
	- P <Account>	project account
	- l|a|p			specifies the queue where the jobs are submitted; low, alloc, and preimium, respectively
	
### uesMinervaIntersect
This script is run the same way as the uesSingleIntersect or uesBatchIntersect scripts. Run the script from the main folder of your UES run. It will send the output to the analysis/ folder. As with the uesMinerva script, there are 2 additional flags required:
    -P <Account>	project account
    -l|a|p			specifies the queue where the jobs are submitted; low, alloc, and preimium, respectively
	
### uesMinervaCleanUp
Since the uesMinervaIntersect script runs the intersection analysis asynchronously, the output must be sorted after the fact. Once all of the jobs have completed, run the uesMinervaCleanUp script to sort the output and provide the proper headings. This script takes one argument: the file summaryStats.YOUR-FILE-NAME.txt.temp. The file will be renamed without the .temp once sorting is complete.
