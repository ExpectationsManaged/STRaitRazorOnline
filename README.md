# STRaitRazorOnline

STRait Razor Online (SRO) serves as the user-interface (UI) for analyzing sequencing data with STRait Razor. Using the SRO, converts STRait Razor* sequence-based allele calls into genotype tables for and/or bar plots for downstream analysis. STRait Razor Online also allows the user to control the genotype calls using several thresholds (e.g., heterozygous balance, strand balance, etc.). Read strands are merged and reported as single strand (Default = Forward) in Genotype Tables.

*STRait Razor use is detailed in separate manual. Use of this workbook assumes a degree of STRait Razor knowledge prior to use.


# Installation

#Note to macOS users: Please install Xquartz from https://www.xquartz.org/ and then logout and log back in.

1)	Download or clone the repository.
2)	Unzip the STRaitRazorOnline.zip file. 
	
	a.	Linux users: rename str8rzr_linux to str8rzr
	b.	macOS users: rename str8rzr_osX to str8rzr

4)	If you don’t have R installed, install R from your mirror of choice https://www.r-project.org/
5)	If you don’t have RStudio installed, install RStudio Desktop from https://rstudio.com/products/rstudio/download/#download
6)	Once you have R and RStudio installed, a few packages need to be installed.
    #tidyverse; Biostrings; shinydashboard; rhandsontable; tcltk; plotly; stringi; data.table
    a.	Note: These packages are Imports not Suggests. So, make sure you get them all. 😊
6)	Once all the packages are installed, open ‘app_standalone.R’ script from the same directory.

# Landing Page
1.	Load the database files into the environment using the ‘Load DB’ Button.
2.	In sidebar, use the dropdown under ‘Select Kit’ based on the amplification primers* used.
    *Note: ForenSeq PMB is the default state for analysis.
3.	Select ‘Single’ or ‘Batch’ Tab for analysis of one “sample” or >= 2 “samples”.
    Note: the process of compiling the list of fastq files can occur recursively or not (default = TRUE). If you have subdirectories with the same sample name, results may be overwritten ☹. See ‘Settings’ tab on sidebar to change. 
4.	STRait Razor Online accepts both fastq/fastq.gz and allsequences.txt files from previous analyses to reduce analysis time of “reprocessing” fastq files. If processing allsequences.txt, please select ‘v2s’ or ‘v3’ format.

# Data Analysis
 **For more information on the data analysis pipeline see the PDF manual embedded in the app sidebar**
 
# Settings
While this page will look significantly different with the final batch of settings added, this implementation controls a few different functions.

a.	Number of cores for v3: Allocate system resources for fastq processing using str8rzr program

b.	Global Read Depth Thresh: Filter out haplotypes with fewer than X reads

c.	Recursive Directories: When processing a batch of files (either fastq or allsequences), process chosen directory or directory plus subdirectories

d.	AutoPass Loci: When processing single files, “passing” green loci are automatically moved to final data frame and “warning” orange loci are passed to ‘Data Analysis’ tab for interpretation (with some more optimization, will likely move this from Default = False  True)

e.	Save Bar Plots: To save a .png of all loci bar plot to sample folder additionally a conditional setting for separating STR and SNP loci into separate image files

# Troubleshooting
#
•	Issue: When launching electron app  , users are met with blank screen.

•	Solution: View --> Force Reload
#
•	Issue: When launching electron app  , users are met with error message “An error has occurred”.

•	Solution: Open Task Manager and close any instances of electron-quick-start or R for Windows front-end. After these processes have been closed, relaunch the application.
#
•	Issue: When selecting Settings --> AutoPass Loci, Pass Loci Count changes to zero.

•	Solution: View --> Force Reload

# Appendix: STRaitRazorAnalysisConfig
For this implementation, we will be focusing on a subset of the columns.

Kit: Amplification kit(s) used for target enrichment*
  *Note: ForenSeq PMB is the default state for analysis.
  
Locus: Full list of markers in each kit

Repeat_Size: Period of the repeat (e.g., CSF1PO: ATCT; 4 base repeat)

HBT: Heterozygote Balance Threshold on a per marker & per kit basis. This is used for assignment of second allele prior to data frame passing to UI

  Danger_math_ahead: The threshold is calculated by dividing the second largest, in terms of coverage or read depth, by the largest allele (e.g., the allele 12 has 932 reads associated with it and allele 17.3 has 950 reads. Thus, the heterozygote balance is 0.98). However, the heterozygote balance output in the genotype tables is calculated by dividing the largest (lexicographically) allele by the second largest allele (e.g., G/A. or  14/10 ).

RDT: Read Depth Threshold on a per marker & per kit basis filter prior to data frame passing to UI

SBT: Strand Balance Threshold on a per marker & per kit basis filter prior to data frame passing to UI

RAPT: Relative Allele Proportion Threshold on a per marker & per kit basis filter prior to data frame passing to UI

AutoRAPT: Locus flagging variable. In this beta, loci with a proportion of called alleles above this value will be flagged “green” (e.g., CSF1PO; 1000 reads aligned to the locus, 12: 450, 14: 450, AutoRAPT = 0.857; 900/1000 == 0.90; 0.9 > 0.857, CSF1PO = “green”)

Marker_Type: Categories of marker type (e.g., SNP, microhaplotype, STR) used for parsing markers


# Changelog
v 0.1.2: 07/20/2020
	-Launch Day!!!!!!!!!!!

v 0.1.3: 08/19/2020
	-minor bug fixes

v 0.1.4: 09/01/2020
	-added STRidER tab
	-STRidER input fields
	-added CV to RAP batch output
	-fixed bug affecting SNP loci with secondary SNP variant in matching sequence
	-added progress bar for batch mode or file processing

v 0.1.5: 09/17/2020
	-added data.table package to address STRidER appending issue
	-finalized STRidER functions for multi-sample processing
	-cleaned up in script comments regarding page titles
	-removed readxl package (not used in current implementation)
-added code to clean up UI on start-up
	-Launched Windows Electron application
	-added function subThresh for saving data frame of reads >= GlobalThreshold, but < locus threshold
	-added Indels (e.g., Amelogenin) to All Loci and Final Profile ggplots

v 0.1.51 (app_online.R only): 09/18/2020
	-Corrected bug for fastq processing related to relative path of configs in unzipped fastq pipeline (app_online.R only)
	
v 0.1.6 
	Calculation of heterozygote balance for isoalleles
	Toggle switch for expected vs. observed
	Hold 100% upon restart
	Showing locus name for loci missing data
	Report allele rather than numeric in table
	Removed closed beta tabs
	
v 0.1.7 (config update...no code): 01/05/2021
	Updated recommended PowerSeq config v2-->v2.1

v 0.1.8 (05/10/2021)
	-Changed Kit IDs
		ForenSeq --> ForenSeq DNA Signature
		PowerSeq --> PowerSeq 46GY
	-Added Kit ForenSeq MainstAY
	-Updated config files
		ForenSeq DNA Signature: ForenSeqv1.25 --> ForenSeqv1.26
		GlobalFilerNGSv2: GFNGSv2_v7 --> GFNGSv2_v7.1
		PowerSeq 46GY: PowerSeqv2.1 --> PowerSeqv3.1
	-Updated Sample Name bug in All Loci Tab to display sample name rather than locus
	
v 0.1.9 (config update...no code): 07/23/2021
	-Added Kit IDseek SNP85
	
v 0.2: 08/19/2021
	-Added Locus to haplotype for revComp filter to account for repeat region similarities

v 0.2.1: 10/04/2021	
	-Added hotfix for processing multi-format data alongside single-format (i.e., STRs with STR & SNP data in batch processing) credit: S.P. from VCU
