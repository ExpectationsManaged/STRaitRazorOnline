##############################################
#
#
#Version ID: 0.2.8
#
#
##############################################
# options(repos = BiocManager::repositories())
# options(repos = c(CRAN = "https://cran.rstudio.com/",
#                   BioCsoft = "https://bioconductor.org/packages/3.12/bioc",
#                   BioCann = "https://bioconductor.org/packages/3.12/data/annotation",
#                   BioCexp = "https://bioconductor.org/packages/3.12/data/experiment",
#                   BioCworkflows = "https://bioconductor.org/packages/3.12/workflows",
#                   BioCbooks = "https://bioconductor.org/packages/3.12/books")
#         )

library(data.table)
library(shinydashboard)
library(rhandsontable)
library(tidyverse)
library(shiny)
library(Biostrings)
library(stringi)




##############################################
#Core data processing functions and a global variable

Palette272 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4", "#f9f9f9", "#2a2526", "#e5e6e7", "#6b6a6b", "#a6a4a4")

split_path <- function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

STRaitRazorv3<- function(ActiveFASTQ, SampleName, numcores = 6, Kit = c("ForenSeq", "GlobalFilerNGSv2", "mtDNA", "PowerSeq"), KitPath){

  ActiveFASTQ <- paste0('"', ActiveFASTQ, '"')
  
  #Derive sample names
  fileType <- ifelse(substr(SampleName, nchar(SampleName) - 5, nchar(SampleName)) == ".fastq", ".fastq", ".fastq.gz")
  
  SampleName <- if(fileType == ".fastq"){
    sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(SampleName))
  }else{
    sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(SampleName))))
  }
  
  #Assign kit config
  KitPaths <- read_csv(KitPath)
  TempDF <- data.frame("Kitid" = Kit)
  TempKit <- left_join(TempDF, KitPaths, by = "Kitid")
  v3_Config <- paste0(TempKit[1, 2])
  
  Sys.chmod("./bin/str8rzr", mode="777")
  
  #Call STRaitRazor    
  if (.Platform$OS.type == "windows"){
    if (fileType == ".fastq"){
      STR8RZR_Call <- paste0("bin\\str8rzr -c  db\\configs\\", v3_Config, " -p ", numcores, " ", ActiveFASTQ)
      STRaitRazorIO <- read.table(text = shell(STR8RZR_Call, intern = TRUE), sep = "\t")
    }else{
      TEMP_STR8RZROutput <- shell(paste0('bin\\7z.exe x -so ', ActiveFASTQ, ' | bin\\str8rzr -c db\\configs\\', v3_Config,' -p ', numcores), intern = TRUE)
      STRaitRazorIO <- read.table(text = tail(TEMP_STR8RZROutput, length(TEMP_STR8RZROutput)), sep = "\t")
    }
  }
  
  if (.Platform$OS.type == "unix"){
    if (fileType == ".fastq"){
      STR8RZR_Call <- paste0("./bin/str8rzr -c db/configs/", v3_Config, " -p ", numcores, " ", ActiveFASTQ)
      STRaitRazorIO <- read.table(text = system(STR8RZR_Call, intern = TRUE), sep = "\t")
    }else{
      TEMP_STR8RZROutput <- system(paste0('zcat ', ActiveFASTQ, ' | ./bin/str8rzr -c db/configs/', v3_Config,' -p ', numcores), intern = TRUE)
      STRaitRazorIO <- read.table(text = tail(TEMP_STR8RZROutput, length(TEMP_STR8RZROutput)), sep = "\t")
    }
  }
  
  return(list(SampleName, STRaitRazorIO))
}

convert2pre3 <- function(STRaitRazorIO, Config) {
  
  #Rename column headers   
  names(STRaitRazorIO) <- c("LocusAllele", "AlleleLength", "Haplotype", "HaplotypeCount", "HaplotypeCountRC")    
  
  #Breakout alleles, ditch allele length, and convert Amelogenin alleles
  STRaitRazorIO <- separate(data = STRaitRazorIO, into = c("Locus", "Allele"), col = LocusAllele, sep = ":") %>% 
    select(-(AlleleLength))
  
  STRaitRazorIO <- STRaitRazorIO %>% 
    mutate(Real = (ifelse(str_detect(Allele, "-"), 0, 1))) %>% 
    filter(Real == 1) %>% 
    select(-(Real))
  
  STRaitRazorIO <- type_convert(STRaitRazorIO)
  
  #Enrich for STRs
  Type <- select(Config, Locus, Marker_Type)
  
  STRaitRazorIO <- left_join(STRaitRazorIO, Type, by = "Locus") %>% 
    select(Marker_Type, everything())
  
  STRaitRazorIO$Haplotype <- as.character(STRaitRazorIO$Haplotype)
  
  #Create reverse complement haplotypes (don't @ me I do realize there's a RC function)
  STRaitRazorIO$HaplotypeRC <- chartr("AGTC", "TCAG", Biostrings::reverse(STRaitRazorIO$Haplotype))
  
  #Remove variable
  rm(Type)
  
  return(STRaitRazorIO)
}

idcomphaps <- function(STRaitRazorIO, Config) {
  
  #Rename column headers  
  names(STRaitRazorIO) <- c("LocusAllele", "AlleleLength", "Haplotype", "HaplotypeRC", "HaplotypeCount")    
  
  #Breakout alleles, ditch allele length, and convert Amelogenin alleles
  STRaitRazorIO <- separate(data = STRaitRazorIO, into = c("Locus", "Allele"), col = LocusAllele, sep = ":") %>% 
    select(-(AlleleLength))
  
  X_Allele <- (STRaitRazorIO$Allele[STRaitRazorIO$Allele == "X"] <- 0)
  Y_Allele <- (STRaitRazorIO$Allele[STRaitRazorIO$Allele == "Y"] <- 6)
  
  #Enrich for STRs
  Type <- select(Config, Locus, Marker_Type)
  
  STRaitRazorIO <- left_join(STRaitRazorIO, Type, by = "Locus")
  
  #Merge locus and haplotype in case someone wants to do repeat region only
  STRaitRazorIO <- STRaitRazorIO %>% 
    unite(LH, Locus, Haplotype, remove = FALSE) %>% 
    unite(LHRC, Locus, HaplotypeRC, remove = FALSE)
  
  #Create tables for each haplotype
  LHRC_RD <- select(STRaitRazorIO, LHRC, HaplotypeCount) %>% 
    dplyr::rename(`HaplotypeCountRC` = `HaplotypeCount`, `LH` = `LHRC`)
  
  LH <- select(STRaitRazorIO, LH, HaplotypeCount)
  
  #Join two tables and remove temporary variables
  ljh <- left_join(LH, LHRC_RD, by = "LH")
  ljh[is.na(ljh)] <- 0
  STRaitRazorIO <- left_join(STRaitRazorIO, ljh, by = c("LH", "HaplotypeCount")) %>% 
    select(Marker_Type, Locus, Allele, Haplotype, HaplotypeRC, HaplotypeCount, HaplotypeCountRC)
  
  #Remove variable
  rm(Type)
  
  #Convert data types
  STRaitRazorIO <- type_convert(STRaitRazorIO)  
  
  return(STRaitRazorIO)
}

STRaitRazorFilter <- function(allsequences, SampleName, Kit, GlobalRDT = 2, mkallseq_version = "v3", DBPath){
  
  #Identify paths and read in database variables
  DatabasePaths <- read_csv(DBPath)
  Config <- read_csv(file.path("db", "db", DatabasePaths[[2, 2]])) %>% 
    filter(Kitid == Kit)
  HaplotypeDatabase <- read_csv(file.path("db", "db", DatabasePaths[[3, 2]])) %>% 
    filter(Kitid == Kit)
  ResultsDir <- paste0(file.path("data", DatabasePaths[[4, 2]]), "/", Kit)
  OutputDir <- paste0(DatabasePaths[[5, 2]])
  
  #Create folder for results
  if (file.exists(OutputDir)){
    ""
  }else{
    dir.create(file.path(ResultsDir, OutputDir), showWarnings = FALSE, recursive = TRUE)
  }
  
  #Load STRait Razor result
  STRaitRazorIO <- allsequences
  attach(STRaitRazorIO)
  
  #Enrich for target strand and arrange data
  if (mkallseq_version == "v2s"){
    STRaitRazorIO <- idcomphaps(STRaitRazorIO, Config)
  }else{
    STRaitRazorIO <- convert2pre3(STRaitRazorIO, Config)
  }
  
  #Isolate desired strand and columns
  Motif <- select(Config, Locus, InternalMotif1, InternalMotif2, InternalMotif3, Repeat_Size)
  STRaitRazorIO <- left_join(STRaitRazorIO, Motif, by = "Locus")
  
  STRaitRazorIO <- STRaitRazorIO %>% 
    mutate(Keep = (ifelse(str_detect(STRaitRazorIO$Haplotype, STRaitRazorIO$InternalMotif1), 1, 0) + ifelse(str_detect(STRaitRazorIO$Haplotype, STRaitRazorIO$InternalMotif2), 1, 0) + ifelse(str_detect(STRaitRazorIO$Haplotype, STRaitRazorIO$InternalMotif3), 1, 0))) %>% 
    mutate(KeepRC = (ifelse(str_detect(STRaitRazorIO$HaplotypeRC, STRaitRazorIO$InternalMotif1), 1, 0) + ifelse(str_detect(STRaitRazorIO$HaplotypeRC, STRaitRazorIO$InternalMotif2), 1, 0) + ifelse(str_detect(STRaitRazorIO$HaplotypeRC, STRaitRazorIO$InternalMotif3), 1, 0)))
  
  STRaitRazorIO <- STRaitRazorIO %>% 
    mutate(FinalHaplotype = ifelse(STRaitRazorIO$Keep >= STRaitRazorIO$KeepRC, STRaitRazorIO$Haplotype, STRaitRazorIO$HaplotypeRC)) %>% 
    mutate(FinalHaplotype2 = paste0(Locus, "_", FinalHaplotype))
  
  STRaitRazorIO <- STRaitRazorIO %>% 
    mutate(ContainsComp = ifelse((duplicated(STRaitRazorIO$FinalHaplotype2) | duplicated(STRaitRazorIO$FinalHaplotype2, fromLast = TRUE)), 1, 0))
  
  STRaitRazorIO <- STRaitRazorIO %>% 
    mutate(DesiredStrand = ifelse(ContainsComp == 0, 1, ifelse(Keep >= KeepRC, 1, 0))) %>% 
    filter(DesiredStrand == 1) %>% 
    select(Locus, Allele, FinalHaplotype, HaplotypeCount, HaplotypeCountRC, Repeat_Size) %>% 
    dplyr::rename(`Haplotype` = `FinalHaplotype`)
  
  #Calculate locus depth, strand balance
  STRaitRazorIO <- STRaitRazorIO %>% 
    mutate(SB = ifelse(HaplotypeCount == 0 | HaplotypeCountRC == 0, 0, HaplotypeCount/HaplotypeCountRC)) %>% 
    mutate(SBC = ifelse(SB > 1, 1/SB, SB)) %>%
    mutate(HaplotypeSum = (HaplotypeCount + HaplotypeCountRC))
  
  #Temporary merge locus and haplotype for DB
  HaplotypeDatabase <- HaplotypeDatabase %>% 
    unite(LH, Locus, Haplotype, remove = FALSE)
  STRaitRazorIO <- STRaitRazorIO %>% 
    unite(LH, Locus, Haplotype, remove = FALSE)
  
  #Merge with database
  HaplotypeDatabase <- select(HaplotypeDatabase, LH, Nomenclature, Allele_Char)
  STRaitRazorIO <- left_join(STRaitRazorIO, HaplotypeDatabase, by = "LH") %>% 
    select(-(LH))
  
  #Consolidate alleles for reporting
  STRaitRazorIO <- STRaitRazorIO %>% 
    mutate(FinalAllele = ifelse(is.na(Allele_Char), as.character(Allele), Allele_Char)) %>% 
    select(-(Allele_Char))
  
  #Combine Locus and CE state into nomenclature
  STRaitRazorIO$Nomenclature <- ifelse(is.na(STRaitRazorIO$Nomenclature), "", paste0(STRaitRazorIO$Locus," [CE ", as.character(STRaitRazorIO$FinalAllele),"]", STRaitRazorIO$Nomenclature)) 
  rm(HaplotypeDatabase)
  
  #Extract summary data
  LocusSummary <- summarize(dplyr::group_by(STRaitRazorIO, Locus), UniqueHaps = n(), TotalReads = sum(HaplotypeSum))
  AlleleSummary <- summarize(dplyr::group_by(STRaitRazorIO, Locus, Allele), UniqueHaps = n(), TotalReads = sum(HaplotypeSum))    
  
  #Apply global threshold
  STRaitRazorIO <- STRaitRazorIO %>% 
    filter(HaplotypeSum >= GlobalRDT)
  
  #Check if all haplotypes have been removed    
  if(nrow(STRaitRazorIO) == 0){
    
    #Sample subdirectories
    if (file.exists(file.path(ResultsDir, OutputDir, SampleName))){
      ""
    }else{
      dir.create(file.path(ResultsDir, OutputDir, SampleName), showWarnings = FALSE, recursive = TRUE)
    }
    OutputPath <- file.path(ResultsDir, OutputDir, SampleName)
    
    # #Save tables
    # write.csv(as.data.frame(STRaitRazorIO), paste0(OutputPath, "/", SampleName, "_DataTable.csv"), row.names = FALSE)
    # write.csv(as.data.frame(AlleleSummary), paste0(OutputPath, "/", SampleName, "_AlleleSummary.csv"), row.names = FALSE)
    # write.csv(as.data.frame(LocusSummary), paste0(OutputPath, "/", SampleName, "_LocusSummary.csv"), row.names = FALSE)
    
    #Return variables
    return(list(STRaitRazorIO, AlleleSummary, LocusSummary, SampleName, OutputPath))
    
  }else{ 
    #Relative Allele Proportion
    STRaitRazorIO <- left_join(STRaitRazorIO, LocusSummary, by = "Locus")
    
    STRaitRazorIO <- STRaitRazorIO %>% 
      mutate(RAP = (HaplotypeSum/TotalReads)) %>% 
      select(-(UniqueHaps), -(TotalReads), -(HaplotypeCount), -(HaplotypeCountRC), -(Repeat_Size))
    
    #Set thresholds
    Thresholds <- select(Config, Locus, HBT, RDT, SBT, RAPT)
    STRaitRazorIO <- left_join(STRaitRazorIO, Thresholds, by = "Locus")
    rm(Thresholds)
    
    #Apply thresholds and rank haplotype abundance
    STRaitRazorIO <- filter(STRaitRazorIO, HaplotypeSum >= RDT, SBC >= SBT, RAP >= RAPT) %>% 
      select(-(RDT), -(SBT), -(RAPT), -(SBC)) %>% 
      dplyr::group_by(Locus) %>% 
      mutate(LocusRank = rank(-HaplotypeSum, ties.method = "first"))
    
    STRaitRazorIO <- STRaitRazorIO[with(STRaitRazorIO, order(Locus, LocusRank)),]
    
    #Error handling for low-template DNA samples
    if(nrow(STRaitRazorIO) == 0){
      
      #Sample subdirectories
      if (file.exists(file.path(ResultsDir, OutputDir, SampleName))){
        ""
      }else{
        dir.create(file.path(ResultsDir, OutputDir, SampleName), showWarnings = FALSE, recursive = TRUE)
      }
      
      OutputPath <- file.path(ResultsDir, OutputDir, SampleName)
      
      # #Save tables
      # write.csv(as.data.frame(STRaitRazorIO), paste0(OutputPath, "/", SampleName, "_DataTable.csv"), row.names = FALSE)
      # write.csv(as.data.frame(AlleleSummary), paste0(OutputPath, "/", SampleName, "_AlleleSummary.csv"), row.names = FALSE)
      # write.csv(as.data.frame(LocusSummary), paste0(OutputPath, "/", SampleName, "_LocusSummary.csv"), row.names = FALSE)
      
      #Return variables
      return(list(STRaitRazorIO, AlleleSummary, LocusSummary, SampleName, OutputPath))
      
    }else{

      #Create abundance ratio
      STRaitRazorIO <- left_join(STRaitRazorIO, select(mutate(left_join(dplyr::rename(aggregate(HaplotypeSum ~ Locus, STRaitRazorIO, function(x) sort(x, decreasing = T)[1]), First = HaplotypeSum), dplyr::rename(aggregate(HaplotypeSum ~ Locus, STRaitRazorIO, function(x) sort(x, decreasing = T)[2]), Second = HaplotypeSum), by = "Locus"), AR = Second/First), Locus, AR), by = "Locus")
      
      #Single haplotype adjustment
      STRaitRazorIO$AR <- ifelse(is.na(STRaitRazorIO$AR), 0, STRaitRazorIO$AR)
      
      #Auto call alleles
      STRaitRazorIO$AutoCall <- ifelse(STRaitRazorIO$AR >= STRaitRazorIO$HBT & STRaitRazorIO$LocusRank <= 2, "a", ifelse(STRaitRazorIO$AR < STRaitRazorIO$HBT, ifelse(STRaitRazorIO$LocusRank == 1, "a", ""), ""))
      
      #Clean up heterozygote balance
      STRaitRazorIO$AR <- ifelse(STRaitRazorIO$AutoCall == "a" & STRaitRazorIO$LocusRank <= 2 & STRaitRazorIO$AR >= as.numeric(STRaitRazorIO$HBT), STRaitRazorIO$AR, "")
      
      #Spread data frame to determine balance direction
      TempBalance <- filter(select(STRaitRazorIO, Locus, FinalAllele, LocusRank, AutoCall), AutoCall == "a") %>% 
        spread(LocusRank, FinalAllele)
      
      #This accounts for homozygous only profiles (i.e., low-template DNA)
      if(length(TempBalance) == 3){
        TempBalance$AB <- 1
      }else{
        TempBalance$AB <- ifelse(is.na(TempBalance$`2`), "", ifelse(is.na(as.numeric(TempBalance$`1`)), ifelse(stri_cmp_lt(TempBalance$`1`, TempBalance$`2`), 1, 0), ifelse(as.numeric(TempBalance$`1`) < as.numeric(TempBalance$`2`), 1, 0)))
      }
      
      STRaitRazorIO <- left_join(STRaitRazorIO, select(TempBalance, Locus, AB), by = "Locus")
      STRaitRazorIO$AB <- ifelse(STRaitRazorIO$AR > 0, STRaitRazorIO$AB, "")
      
      #Calculate heterozygote balance
      STRaitRazorIO <- STRaitRazorIO %>% 
        type_convert() %>% 
        mutate(HB = ifelse(AB == 0, (1/as.numeric(AR)), ifelse(AB == 1, AR, ""))) %>% 
        select(-(AB), -(HBT))
      
      #Convert final allele to character to account for SNP-less SNP data
      STRaitRazorIO$FinalAllele <- as.character(STRaitRazorIO$FinalAllele)
      
      #Sample subdirectories
      if (file.exists(file.path(ResultsDir, OutputDir, SampleName))){
        ""
      }else{
        dir.create(file.path(ResultsDir, OutputDir, SampleName), showWarnings = FALSE, recursive = TRUE)
      }
      
      OutputPath <- file.path(ResultsDir, OutputDir, SampleName)
      
      # #Save tables
      # write.csv(as.data.frame(STRaitRazorIO), paste0(OutputPath, "/", SampleName, "_DataTable.csv"), row.names = FALSE)
      # write.csv(as.data.frame(AlleleSummary), paste0(OutputPath, "/", SampleName, "_AlleleSummary.csv"), row.names = FALSE)
      # write.csv(as.data.frame(LocusSummary), paste0(OutputPath, "/", SampleName, "_LocusSummary.csv"), row.names = FALSE)
      
      #Return variables
      return(list(STRaitRazorIO, AlleleSummary, LocusSummary, SampleName, OutputPath))                          
      
    }
  }}

STRaitRazorFigures <- function(STRaitRazorIO, AlleleSummary, LocusSummary, GlobalRDT = 2, SampleName, OutputPath){
  #Generate color palette
  colorCount = max(STRaitRazorIO$LocusRank)
  randomColor = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  #Plot sequence-based bar plot
  if(colorCount > 272){
    Electrofakogram <- ggplot(STRaitRazorIO, aes(Allele, HaplotypeSum, fill = as.factor(LocusRank))) + geom_bar(stat = "identity", show.legend = FALSE) + facet_wrap(~Locus, scales = "free_x") + theme_classic() + scale_fill_manual(values = sample(randomColor,colorCount)) + theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_text(vjust = 0.5, hjust = 0.5)) + labs(title = paste0(SampleName, " Sequence-based Barplot"), y = "Read Depth", x = "Alleles") + geom_text(aes(x = Allele, y = (0-(max(HaplotypeSum)/10)), label = round(Allele, 1), angle = 90, hjust = 0.5, vjust = 0.5), size = 2) + geom_hline(yintercept = 0) + scale_y_continuous(expand = c(0.15,0))
    
  }else{
    Electrofakogram <- ggplot(STRaitRazorIO, aes(Allele, HaplotypeSum, fill = as.factor(LocusRank))) + geom_bar(stat = "identity", show.legend = FALSE) + facet_wrap(~Locus, scales = "free_x") + theme_classic() + scale_fill_manual(values = Palette272) + theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_text(vjust = 0.5, hjust = 0.5)) + labs(title = paste0(SampleName, " Sequence-based Barplot"), y = "Read Depth", x = "Alleles") + geom_text(aes(x = Allele, y = (0-(max(HaplotypeSum)/10)), label = round(Allele, 1), angle = 90, hjust = 0.5, vjust = 0.5), size = 2) + geom_hline(yintercept = 0) + scale_y_continuous(expand = c(0.15,0))      
  }
  
  # ggsave(paste0(OutputPath, "/", SampleName, "/", SampleName, "_SB_Electrofakogram.jpg"), plot = Electrofakogram, width = 24, height = 24, dpi = 320)
  
  #Plot locus summary
  LocusSummaryPlot <- ggplot(LocusSummary) + geom_bar(aes(reorder(Locus, TotalReads), TotalReads), stat = "identity", show.legend = FALSE) + geom_point(aes(reorder(Locus, TotalReads), UniqueHaps)) + theme_classic() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(hjust = 0.5, vjust = 0.5)) + labs(title = paste0(SampleName," Locus Summary"), x = "Locus", y = "Total Reads") + coord_flip() + scale_y_continuous(expand = c(0, 0))
  
  # ggsave(paste0(OutputPath, "/", SampleName, "/", SampleName, "_LocusSummaryPlot.jpg"), plot = LocusSummaryPlot, width = 9, height = 24, dpi = 320)
  
  #Plot length-based bar plot
  ElectroFakogram_LB <- ggplot(filter(AlleleSummary, TotalReads >= GlobalRDT), aes(Allele, TotalReads)) + geom_bar(stat = "identity") + facet_wrap(~Locus, scales = "free_x") + theme_classic() + theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_text(vjust = 0.5, hjust = 0.5)) + labs(title = paste0(SampleName, " Length-based Barplot"), y = "Read Depth", x = "Alleles") + geom_text(aes(x = Allele, y = (0-(max(TotalReads)/10)), label = round(Allele, 1), angle = 90, hjust = 0.5, vjust = 0.5), size = 2) + geom_hline(yintercept = 0) + scale_y_continuous(expand = c(0.15,0))
  
  # ggsave(paste0(OutputPath, "/", SampleName, "/", SampleName, "_LB_ElectroFakogram.jpg"), plot = ElectroFakogram_LB, width = 24, height = 24, dpi = 320)
}

STRaitRazorFigures_byType <- function(STRaitRazorIO, AlleleSummary, LocusSummary, GlobalRDT = 2, SampleName, OutputPath){
  
  #Generate color palette
  colorCount = max(STRaitRazorIO$LocusRank)
  randomColor = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  
  #Save STR figures 
  if(colorCount > 272){
    Electrofakogram <- ggplot(filter(STRaitRazorIO, Marker_Type == "STR"), aes(Allele, HaplotypeSum, fill = as.factor(LocusRank))) + geom_bar(stat = "identity", show.legend = FALSE) + facet_wrap(~Locus, scales = "free_x") + theme_classic() + scale_fill_manual(values = sample(randomColor,colorCount)) + theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_text(vjust = 0.5, hjust = 0.5)) + labs(title = paste0(SampleName, " Sequence-based STR Barplot"), y = "Read Depth", x = "Alleles") + geom_text(aes(x = Allele, y = (0-(max(HaplotypeSum)/10)), label = round(Allele, 1), angle = 90, hjust = 0.5, vjust = 0.5), size = 2) + geom_hline(yintercept = 0) + scale_y_continuous(expand = c(0.15,0))
    
  }else{
    Electrofakogram <- ggplot(filter(STRaitRazorIO, Marker_Type == "STR"), aes(Allele, HaplotypeSum, fill = as.factor(LocusRank))) + geom_bar(stat = "identity", show.legend = FALSE) + facet_wrap(~Locus, scales = "free_x") + theme_classic() + scale_fill_manual(values = Palette272) + theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_text(vjust = 0.5, hjust = 0.5)) + labs(title = paste0(SampleName, " Sequence-based STR Barplot"), y = "Read Depth", x = "Alleles") + geom_text(aes(x = Allele, y = (0-(max(HaplotypeSum)/10)), label = round(Allele, 1), angle = 90, hjust = 0.5, vjust = 0.5), size = 2) + geom_hline(yintercept = 0) + scale_y_continuous(expand = c(0.15,0))      
  }
  
  # ggsave(paste0(OutputPath, "/", SampleName, "/", SampleName, "_SB_STR_Electrofakogram.jpg"), plot = Electrofakogram, width = 12, height = 12, dpi = 320)
  
  LocusSummaryPlot <- ggplot(LocusSummary) + geom_bar(aes(reorder(Locus, TotalReads), TotalReads), stat = "identity", show.legend = FALSE) + geom_point(aes(reorder(Locus, TotalReads), UniqueHaps)) + theme_classic() + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(hjust = 0.5, vjust = 0.5), axis.text.y = element_text(size = 3)) + labs(title = paste0(SampleName," Locus Summary"), x = "Locus", y = "Total Reads") + coord_flip() + scale_y_continuous(expand = c(0, 0))
  
  # ggsave(paste0(OutputPath, "/", SampleName, "/", SampleName, "_LocusSummaryPlot.jpg"), plot = LocusSummaryPlot, width = 9, height = 12, dpi = 320)
  
  ElectroFakogram_LB <- ggplot(filter(AlleleSummary, TotalReads >= GlobalRDT, Marker_Type == "STR"), aes(Allele, TotalReads)) + geom_bar(stat = "identity") + facet_wrap(~Locus, scales = "free_x") + theme_classic() + theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_text(vjust = 0.5, hjust = 0.5)) + labs(title = paste0(SampleName, " Length-based STR Barplot"), y = "Read Depth", x = "Alleles") + geom_text(aes(x = Allele, y = (0-(max(TotalReads)/10)), label = round(Allele, 1), angle = 90, hjust = 0.5, vjust = 0.5), size = 2) + geom_hline(yintercept = 0) + scale_y_continuous(expand = c(0.15,0))
  
  # ggsave(paste0(OutputPath, "/", SampleName, "/", SampleName, "_LB_STR_ElectroFakogram.jpg"), plot = ElectroFakogram_LB, width = 12, height = 12, dpi = 320)
  
  #Save SNP figures 
  
  #Sequence-based bar plot
  if(colorCount > 272){
    Electrofakogram <- ggplot(filter(STRaitRazorIO, Marker_Type == "SNP"), aes(Allele, HaplotypeSum, fill = as.factor(LocusRank))) + geom_bar(stat = "identity", show.legend = FALSE) + facet_wrap(~Locus) + theme_classic() + scale_fill_manual(values = sample(randomColor,colorCount)) + theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_text(vjust = 0.5, hjust = 0.5)) + labs(title = paste0(SampleName, " Sequence-based SNP Barplot"), y = "Read Depth", x = "Alleles") + geom_text(aes(x = Allele, y = (0-(max(HaplotypeSum)/10)), label = round(Allele, 1), angle = 90, hjust = 0.5, vjust = 0.5), size = 2) + geom_hline(yintercept = 0) + scale_y_continuous(expand = c(0.15,0), limits = c(-20, 20))
    
  }else{
    Electrofakogram <- ggplot(filter(STRaitRazorIO, Marker_Type == "SNP"), aes(Allele, HaplotypeSum, fill = as.factor(LocusRank))) + geom_bar(stat = "identity", show.legend = FALSE) + facet_wrap(~Locus) + theme_classic() + scale_fill_manual(values = Palette272) + theme(axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_text(vjust = 0.5, hjust = 0.5)) + labs(title = paste0(SampleName, " Sequence-based SNP Barplot"), y = "Read Depth", x = "Alleles") + geom_text(aes(x = Allele, y = (0-(max(HaplotypeSum)/10)), label = round(Allele, 1), angle = 90, hjust = 0.5, vjust = 0.5), size = 2) + geom_hline(yintercept = 0) + scale_y_continuous(expand = c(0.15,0)) + scale_x_continuous(limits = c(-20, 20))      
  }
  
  # ggsave(paste0(OutputPath, "/", SampleName, "/", SampleName, "_SB_SNP_Electrofakogram.jpg"), plot = Electrofakogram, width = 24, height = 24, dpi = 320)
  
}

create_PL_Haplotypes_dbs <- function(AlleleSummary, STRaitRazorIO, LocusSummary, SampleName){
  
  #Create db files
  PL_AlleleDepth <- AlleleSummary %>% 
    mutate(ID = paste0(Locus, "_", Allele)) %>% 
    select(ID, TotalReads)
  
  PL_LocusDepth <- LocusSummary %>% 
    select(Locus, TotalReads)  
  
  PL_HaplotypeDepth <- STRaitRazorIO %>%  
    mutate(ID = paste0(Locus, "_", Haplotype)) %>%  
    select(ID, FinalAllele, Nomenclature, HaplotypeSum)
  
  PL_HaplotypeRAP <- STRaitRazorIO %>%  
    mutate(ID = paste0(Locus, "_", Haplotype)) %>%  
    select(ID, FinalAllele, Nomenclature, RAP)  
  
  PL_HaplotypeHB <- filter(STRaitRazorIO, LocusRank == "1") %>% 
    select(Locus, HB)
  
  PL_HaplotypeSB <- STRaitRazorIO %>%  
    mutate(ID = paste0(Locus, "_", Haplotype)) %>%  
    select(ID, FinalAllele, Nomenclature, SB)    
  
  #Tag the column with the name
  colnames(PL_AlleleDepth)[ncol(PL_AlleleDepth)] <- SampleName
  colnames(PL_LocusDepth)[ncol(PL_LocusDepth)] <- SampleName
  colnames(PL_HaplotypeDepth)[ncol(PL_HaplotypeDepth)] <- SampleName
  colnames(PL_HaplotypeRAP)[ncol(PL_HaplotypeRAP)] <- SampleName
  colnames(PL_HaplotypeHB)[ncol(PL_HaplotypeHB)] <- SampleName
  colnames(PL_HaplotypeSB)[ncol(PL_HaplotypeSB)] <- SampleName
  
  #Return STRaitRazorIO   
  return(list(PL_AlleleDepth, PL_LocusDepth, PL_HaplotypeDepth, PL_HaplotypeRAP, PL_HaplotypeHB, PL_HaplotypeSB))
  
}

append_PL_Haplotypes_dbs <- function(AlleleSummary, STRaitRazorIO, LocusSummary, PL_AlleleDepth, PL_LocusDepth, PL_HaplotypeDepth, PL_HaplotypeRAP, PL_HaplotypeHB, PL_HaplotypeSB, SampleName){
  
  #Create temporary db for appending
  temp_PL_AlleleDepth <- AlleleSummary %>% 
    mutate(ID = paste0(Locus, "_", Allele)) %>% 
    select(ID, TotalReads)
  
  temp_PL_LocusDepth <- LocusSummary %>% 
    select(Locus, TotalReads)  
  
  temp_PL_HaplotypeDepth <- STRaitRazorIO %>% 
    mutate(ID = paste0(Locus, "_", Haplotype)) %>% 
    select(ID, FinalAllele, Nomenclature, HaplotypeSum)
  
  temp_PL_HaplotypeRAP <- STRaitRazorIO %>% 
    mutate(ID = paste0(Locus, "_", Haplotype)) %>%  
    select(ID, FinalAllele, Nomenclature, RAP)  
  
  temp_PL_HaplotypeHB <- filter(STRaitRazorIO, LocusRank == "1") %>% 
    select(Locus, HB)
  
  temp_PL_HaplotypeSB <- STRaitRazorIO %>%  
    mutate(ID = paste0(Locus, "_", Haplotype)) %>% 
    select(ID, FinalAllele, Nomenclature, Nomenclature, SB)  
  
  
  #Find new alleles
  new_PL_AlleleDepth <- anti_join(select(temp_PL_AlleleDepth, ID), PL_AlleleDepth, by = "ID")
  PL_AlleleDepth <- bind_rows(PL_AlleleDepth, new_PL_AlleleDepth)
  
  new_PL_LocusDepth <- anti_join(select(temp_PL_LocusDepth, Locus), PL_LocusDepth, by = "Locus")
  PL_LocusDepth <- bind_rows(PL_LocusDepth, new_PL_LocusDepth)  
  
  new_PL_HaplotypeDepth <- anti_join(select(temp_PL_HaplotypeDepth, ID, FinalAllele, Nomenclature), PL_HaplotypeDepth, by = "ID")
  PL_HaplotypeDepth <- bind_rows(PL_HaplotypeDepth, new_PL_HaplotypeDepth)  
  
  new_PL_HaplotypeRAP <- anti_join(select(temp_PL_HaplotypeRAP, ID, FinalAllele, Nomenclature), PL_HaplotypeRAP, by = "ID")
  PL_HaplotypeRAP <- bind_rows(PL_HaplotypeRAP, new_PL_HaplotypeRAP)  
  
  new_PL_HaplotypeHB <- anti_join(select(temp_PL_HaplotypeHB, Locus), PL_HaplotypeHB, by = "Locus")
  PL_HaplotypeHB <- bind_rows(PL_HaplotypeHB, new_PL_HaplotypeHB)
  
  new_PL_HaplotypeSB <- anti_join(select(temp_PL_HaplotypeSB, ID, FinalAllele, Nomenclature), PL_HaplotypeSB, by = "ID")
  PL_HaplotypeSB <- bind_rows(PL_HaplotypeSB, new_PL_HaplotypeSB)
  
  
  #Combine new results with compilation
  PL_AlleleDepth <- left_join(PL_AlleleDepth, select(temp_PL_AlleleDepth, ID, TotalReads), by = "ID")
  PL_AlleleDepth[is.na(PL_AlleleDepth)] <- ""
  
  PL_LocusDepth <- left_join(PL_LocusDepth, select(temp_PL_LocusDepth, Locus, TotalReads), by = "Locus")
  PL_LocusDepth[is.na(PL_LocusDepth)] <- ""  
  
  PL_HaplotypeDepth <- left_join(PL_HaplotypeDepth, select(temp_PL_HaplotypeDepth, ID, HaplotypeSum), by = "ID")
  PL_HaplotypeDepth[is.na(PL_HaplotypeDepth)] <- ""
  
  PL_HaplotypeRAP <- left_join(PL_HaplotypeRAP, select(temp_PL_HaplotypeRAP, ID, RAP), by = "ID")
  PL_HaplotypeRAP[is.na(PL_HaplotypeRAP)] <- ""  
  
  PL_HaplotypeHB <- left_join(PL_HaplotypeHB, select(temp_PL_HaplotypeHB, Locus, HB), by = "Locus")
  PL_HaplotypeHB[is.na(PL_HaplotypeHB)] <- ""  
  
  PL_HaplotypeSB <- left_join(PL_HaplotypeSB, select(temp_PL_HaplotypeSB, ID, SB), by = "ID")
  PL_HaplotypeSB[is.na(PL_HaplotypeSB)] <- ""  
  
  
  
  #Tag the column with the name
  colnames(PL_AlleleDepth)[ncol(PL_AlleleDepth)] <- SampleName
  colnames(PL_LocusDepth)[ncol(PL_LocusDepth)] <- SampleName
  colnames(PL_HaplotypeDepth)[ncol(PL_HaplotypeDepth)] <- SampleName
  colnames(PL_HaplotypeRAP)[ncol(PL_HaplotypeRAP)] <- SampleName
  colnames(PL_HaplotypeHB)[ncol(PL_HaplotypeHB)] <- SampleName
  colnames(PL_HaplotypeSB)[ncol(PL_HaplotypeSB)] <- SampleName
  
  #Remove variables
  rm(temp_PL_AlleleDepth, new_PL_AlleleDepth, temp_PL_LocusDepth, new_PL_LocusDepth, temp_PL_HaplotypeDepth, new_PL_HaplotypeDepth, temp_PL_HaplotypeRAP, new_PL_HaplotypeRAP, temp_PL_HaplotypeHB, new_PL_HaplotypeHB, temp_PL_HaplotypeSB, new_PL_HaplotypeSB)
  
  #Return STRaitRazorIO   
  return(list(PL_AlleleDepth, PL_LocusDepth, PL_HaplotypeDepth, PL_HaplotypeRAP, PL_HaplotypeHB, PL_HaplotypeSB))
}
#########STRidER########
# saveSeqs <- function(df){
#   write_tsv(select(df, -Locus), paste0("data/AnalysisOutput/STRidER/STRidER_Data_File_", unique(df$Locus),".tsv"))
#   return(df)
# }

# STRidER_formatting <- function(longformData, Kit, DBPath){
#   
#   
#   #Identify paths and read in database variables
#   DatabasePaths <- read_csv(file.path("db", "db", DBPath))
#   Config <- read_csv(file.path("db", "db", DatabasePaths[[2, 2]])) %>% 
#     filter(Kitid == Kit)
#   
#   longGTDF <- read_tsv(longformData)
#   
#   longGTDF <- left_join(longGTDF, select(Config, Locus, STRidER_Up))
#   
#   
#   longGTDF$Allele <- ifelse(longGTDF$Locus == "Amelogenin", ifelse(longGTDF$Allele == 0, "X", ifelse(longGTDF$Allele == 1, "Y", ifelse(longGTDF$Allele == 6, "Y", longGTDF$Allele))), longGTDF$Allele)
#   
#   
#   #Filter for the called alleles from the STRidER designated loci for sequence saving
#   STRidERSeqUp <- longGTDF %>% 
#     filter(AutoCall =="a", STRidER_Up == 1) %>% 
#     select(SampleName, Locus, Allele, Haplotype) %>% 
#     mutate(fastaHeader = paste0(">", SampleName, "_", Locus, "_", Allele)) %>% 
#     select(Locus, fastaHeader, Haplotype) %>% 
#     arrange(Locus, fastaHeader)
#   
#   #Save sequnece files grouped by locus
#   STRidERSeqUp %>% 
#     group_by(Locus) %>% 
#     do(saveSeqs(.))
#   
#   
#   #Filter uncalled alleles and create id's for joins
#   longGTDF <- longGTDF %>% 
#     filter(AutoCall =="a", STRidER_Up == 1) %>% 
#     select(SampleName, Locus, Allele) %>% 
#     group_by(SampleName, Locus) %>% 
#     mutate(AlleleRank = rank(Allele, ties.method = "first"), ID = paste0(SampleName,"_", Locus), ID2 = paste0(Locus, "_", AlleleRank)) %>% 
#     arrange(SampleName, Locus, AlleleRank) %>% 
#     ungroup()
#   
#   
#   #Filter for STRidER Loci
#   Config <- Config %>% 
#     filter(STRidER_Up == 1)
#   
#   
#   #Pull list of loci typed
#   Loci <- unique(Config$Locus)
#   
#   
#   #Duplicate loci for genotype table
#   Loci <- append(Loci, Loci) %>% 
#     sort() %>%
#     as.data.frame()
#   
#   names(Loci) <- c("Locus")
#   
#   
#   #Create id's for joining
#   Loci <- Loci %>%
#     mutate(Order = 1:nrow(Loci)) %>% 
#     group_by(Locus) %>% 
#     mutate(Count = row_number()) %>% 
#     mutate(ID2 = paste0(Locus, "_", Count)) %>%
#     ungroup()
#   
#   longGTDF <- left_join(longGTDF, select(Loci, ID2, Order), by = "ID2")
#   
#   longGTDF <- longGTDF %>% 
#     mutate(ID3 = paste0(SampleName, "_", ID2), SisterAlleleID = ID3)
#   
#   
#   #Make a list of all locus-sample combinations
#   SampleList <- as.data.frame(rep(unique(longGTDF$SampleName),nrow(Loci))) 
#   
#   names(SampleList) <- "SampleName"
#   
#   SampleList <- SampleList %>% 
#     arrange(SampleName) %>% 
#     group_by(SampleName) %>% 
#     mutate(Order = row_number()) %>% 
#     ungroup() 
#   
#   SampleList <- left_join(SampleList, select(Loci, Order, Locus), by = "Order") %>% 
#     group_by(SampleName, Locus) %>% 
#     mutate(Count = row_number(), ID3 = paste0(SampleName, "_", Locus, "_", Count)) %>% 
#     ungroup()
#   
#   antiAlleles  <- anti_join(SampleList, longGTDF, by = "ID3") %>% 
#     mutate(SisterAlleleID = if_else(Count == 1, "", paste0(SampleName, "_", Locus, "_1")))
#   
#   antiAlleles <- left_join(antiAlleles, select(longGTDF, SisterAlleleID, Allele), by = "SisterAlleleID")
#   
#   longGTDF_XX<- bind_rows(longGTDF, select(antiAlleles, SampleName, Locus, Allele, Order, ID3)) %>% arrange(SampleName, Order)
#   
#   
#   #Create and merge header and wide form
#   headerButNotHeader <- Loci %>% 
#     mutate(SampleName = "SampleName") %>% 
#     select(SampleName, Locus, Order) %>% 
#     pivot_wider(names_from = Order, values_from = Locus)
#   
#   STRidER_Wide <- longGTDF_XX %>% 
#     select(SampleName, Order, Allele) %>% 
#     pivot_wider(names_from = Order, values_from = Allele)
#   
#   STRidER_Wide <- rbindlist(list(headerButNotHeader, STRidER_Wide), fill = TRUE)
#   
#   colnames(STRidER_Wide)[1] <- "X1"
#   
#   return(STRidER_Wide)
#   
# }


##############################################
#UI script
ui <- dashboardPage(
  
  dashboardHeader(
    title = "STRait Razor Analysis v0.2.7"
  ),

  ##############################################
  #Dashboard sidebar
  dashboardSidebar(     
    
    sidebarMenu(
      id = "tabs",
      
      menuItem("Home", tabName = "home", icon = icon("home"), selected = TRUE),

      menuItem("Kit/Locus Selection", tabName = "pickkit", icon = icon("fingerprint"), startExpanded = TRUE,
        menuSubItem(uiOutput("kits"), icon = NULL),
        menuSubItem(uiOutput("loci"), icon = NULL),
        menuSubItem(checkboxInput("cLocus", "Manual Locus Selection"), icon = NULL)
      ),

      menuItem("Data Analysis", tabName = "SLDA", icon = icon("dna")),
      
      menuItem("Thresholds", tabName = "thresholds", icon = icon("sliders-h"),       
        menuSubItem(sliderInput("RAP_Threshold", "Relative Allele Proportion Threshold", min = 0, max = 1, value = 0), icon = NULL),
        menuSubItem(sliderInput("SB_Threshold", "Strand Balance Threshold", min = 0, max = 1, value = 0), icon = NULL),
        menuSubItem(numericInput("HD_Threshold", "Haplotype Depth Threshold", value = 0, min = 1, max = 99999999), icon = NULL)      
      ),
      
      menuItem("All Loci", tabName = "barplot_all", icon = icon("chart-bar")),
      
      menuItem("Final Profile", tabName = "final_df", icon = icon("save")),
      
      # menuItem("STRidER", tabName = "strider", icon = icon("file-upload")),
      
      menuItem("Manual", tabName = "manual", icon = icon("file-pdf")),
      
      menuItem("Contact", tabName = "email", icon = icon("envelope")),
            
      menuItem("Settings", tabName = "settings", icon = icon("ellipsis-v")),
      
      tags$hr(),
      
      valueBoxOutput("progressBox", width = "100%"),
      
      valueBoxOutput("primaryBoxCount", width = "80%"),
      
      valueBoxOutput("warningBoxCount", width = "80%")
      
    )
  ),
  
  ##############################################
  #Dashboard page body 
  dashboardBody(
    
    ##############################################
    #Apply CSS formatting   
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "STRaitRazoR.css")
    ),    
    
    
    ##############################################
    #Page/tab contents   
    tabItems(
      
      
      ##############################################
      #Home page       
      tabItem(tabName = "home",
              
              fluidPage(
                id = 'landingpage',
                fluidRow(
                  column(
                    width = 12,
                    align = "center",
                    tags$a(
                      href = "https://www.unthsc.edu/graduate-school-of-biomedical-sciences/laboratory-faculty-and-staff/strait-razor/",
                      target = "_blank",
                      tags$img(
                        src = "STRaitRazor_30pxPaddingDNA.png",
                        width = "259",
                        height = "300",
                        align = "center",
                        style = "margin: 10px 10px"
                      )
                    )
                  )
                ),
                
                tags$hr(class = "thiccHr"),
                
                fluidRow(
                  column(
                    offset = 3, 
                    width = 12,
                    align = "left",
                    
                    tags$p(
                      tags$strong("Important update v0.2.7 (07/12/2023): "),
                      "Length-based (CE) allele nomenclature for DYS627 reduced 6 repeats (e.g., 26 --> 20) in accordance with Draft ISFG Recommendations."
                    )
                    
                  )
                ),
                
                tags$hr(class = "thiccHr"),
                fluidRow(
                  id = "dbpath",
                  column(
                    offset = 4, 
                    width = 12,
                    align = "center",
                    box(
                      id = "pathchoosing",
                      title = h3("First Select Database Path", style = "font-family: Times New Roman"), 
                      width = 4,
                      actionButton("loaddb", "Load Database")
                    )
                  )
                ),
                
                tags$hr(class = "thiccHr"),
                
                fluidRow(
                  column(
                    offset = 3, 
                    width = 12,
                    align = "center",
                  
                    tabBox(
                      title = h3("Then Upload Sample", style = "font-family: Times New Roman"),
                      width = 6,
                      id = "fileselect",
                      tabPanel("Single", 
                               div(style = "padding: 0 0 0 0",
                                   radioButtons("versionS", "Version of allsequences.txt", c("v2s", "v3"), selected = "v3", inline = TRUE)),
                               tags$hr(),
                               div(style = "padding: 0 0 0 0",
                                   textInput("sampleID", "Sample Name*", value = ""),
                                   "*Only for allsequences.txt. Will cause crash without file extension."),
                               tags$hr(),
                               div(style = "display: inline-block; padding: 0 10px 0 0", 
                                   fileInput("fastq", "FASTQ File")), 
                               div(style = "display: inline-block; padding: 0 0 0 10px", 
                                   fileInput("allseq", "allsequences.txt File"))
                      )#,
                      # tabPanel("Batch",                         
                      #          div(style = "padding: 0 0 0 0",
                      #              radioButtons("versionB", "Version of allsequences.txt", c("v2s", "v3"), selected = "v3", inline = TRUE)),
                      #          tags$hr(),
                      #          div(style = "display: inline-block; padding: 0 10px 0 0",
                      #              actionButton("fastqs", "FASTQs Dir")),
                      #          div(style = "display: inline-block; padding: 0 0 0 10px",
                      #              actionButton("allseqs", "allsequences Dir"))
                      # )
                    )
                  )
                ),
                
                tags$div(
                  id = "footers",
          ##############################################
          #Apply primary footer
                  
                  tags$footer(
                    
                    tags$a(
                      href = "https://www.unthsc.edu/graduate-school-of-biomedical-sciences/laboratory-faculty-and-staff/",
                      target = "_blank",
                      tags$img(
                        src = "UNTCHI_RDU.png",
                        width = "75",
                        height = "75",
                        align = "left",
                        style = "margin: 10px 10px"
                      )
                    ), 
                    
                    tags$a(
                      href = "https://www.untchi.org/",
                      target = "_blank",
                      tags$img(
                        src = "CHI_2018.png",
                        width = "75",
                        height = "75",
                        align = "right",
                        style = "margin: 10px 10px"
                      )
                    ),
                    
                    tags$div(
                      id = "social-media-icons",
                      
                      tags$br(),
                      
                      tags$a(
                        href ="https://twitter.com/UNTCHI_RDU",
                        target = "_blank",
                        tags$img(src = "twitter.png",
                                 title = "UNTCHI RDU on Twitter",
                                 width = "30",
                                 height = "30",
                                 class = "displayed",
                                 style = "margin: 10px 10px"
                        )
                      ),
                      
                      tags$a(
                        href ="https://www.instagram.com/untchi_rdu/",
                        target = "_blank",
                        tags$img(src = "instagram.png",
                                 title = "UNTCHI RDU on Instagram",
                                 width = "30",
                                 height = "30",
                                 class = "displayed",
                                 style = "margin: 10px 10px"
                        )
                      ),
                      
                      tags$a(
                        href = "https://www.facebook.com/UNTCenterForHumanIdentification",
                        target = "_blank",
                        tags$img(src = "facebook.png",
                                 title = "UNTCHI RDU on Facebook",
                                 width = "30",
                                 height = "30",
                                 class = "displayed",
                                 style = "margin: 10px 10px"
                        )
                      )
                    ),
                    
                    style = "background-color: #00778b;
                    color: #FFFFFF;",
                    
                    HTML("<br><br>")
                  ),
                  
          ##############################################
          #Apply secondary footer    
                  tags$footer(
                    tags$br(),
                    
                    tags$div(
                      id = "secondary-footer-docs"
                    ),
                    
                    tags$br(), 
                    
                    style = "background-color: #253746;
                    color: #FFFFFF;"
                  )
                )
              )
      ),


      ##############################################
      #Data analysis page           
      tabItem(tabName = "SLDA",
        
        #Plot full (post-filter) and called allele plots
        fluidRow(
          uiOutput("plotFullUI"),
          box(title = h4("Called Alleles", style = "font-family: Times New Roman"), solidHeader = TRUE, status = "primary", plotOutput("plotCalled"))
        ),
        
        #Select locus
        fluidRow(
          column(2, 
                 align = "left", 
                 offset = 0, 
                 actionButton("toPlot", label = h4("View Full Profile", style = "font-family: Times New Roman"))),
          column(6, 
                 align = "center", 
                 offset = 1,
                 actionButton("backone", label = h4("Previous Locus", style = "font-family: Times New Roman")),
                 actionButton("resetlocus", label = h6("First Locus", style = "font-family: Times New Roman")),
                 actionButton("upone", label = h4("Next Locus", style = "font-family: Times New Roman"))),
          column(2, 
                 align = "right", 
                 offset = 1, 
                 actionButton("save2df", h4("Save Locus", style = "font-family: Times New Roman")))
        ),

        #Auto-save   
        fluidRow(
          column(2, align = "right", offset = 10, checkboxInput("autoSave", "Autosave Data"))
        ),
        
        #Display post-filter table
        fluidRow(
          rHandsontableOutput("hot")
        ),
        
        
        #Clear final data frame 
        fluidRow(
          column(8, 
                 align = "center",
                 offset = 2,
                 div(style = "display: inline-block; padding: 10px 10px 0 0", 
                    actionButton("clearFinal", 
                      label = h4("Clear Final Data Frame", 
                                 style = "font-family: Times New Roman"))),
                 div(style = "display: inline-block; padding: 10px 0 0 10px", 
                     downloadLink("saveDTalt", 
                              label = h4("Save Final Data Frame", 
                                         style = "font-family: Times New Roman")))
          )
        ),
        
        tags$div(
          id = "footers",
          ##############################################
          #Apply primary footer
      
          tags$footer(
      
            tags$a(
              href = "https://www.unthsc.edu/graduate-school-of-biomedical-sciences/laboratory-faculty-and-staff/",
              target = "_blank",
              tags$img(
                src = "UNTCHI_RDU.png",
                width = "75",
                height = "75",
                align = "left",
                style = "margin: 10px 10px"
              )
            ), 
            
            tags$a(
              href = "https://www.untchi.org/",
              target = "_blank",
              tags$img(
                src = "CHI_2018.png",
                width = "75",
                height = "75",
                align = "right",
                style = "margin: 10px 10px"
              )
            ),
            
            tags$div(
              id = "social-media-icons",
              
              tags$br(),
              
              tags$a(
                href ="https://twitter.com/UNTCHI_RDU",
                target = "_blank",
                tags$img(src = "twitter.png",
                         title = "UNTCHI RDU on Twitter",
                         width = "30",
                         height = "30",
                         class = "displayed",
                         style = "margin: 10px 10px"
                )
              ),
              
              tags$a(
                href ="https://www.instagram.com/untchi_rdu/",
                target = "_blank",
                tags$img(src = "instagram.png",
                         title = "UNTCHI RDU on Instagram",
                         width = "30",
                         height = "30",
                         class = "displayed",
                         style = "margin: 10px 10px"
                )
              ),
              
              tags$a(
                href = "https://www.facebook.com/UNTCenterForHumanIdentification",
                target = "_blank",
                tags$img(src = "facebook.png",
                         title = "UNTCHI RDU on Facebook",
                         width = "30",
                         height = "30",
                         class = "displayed",
                         style = "margin: 10px 10px"
                )
              )
            ),
       
            style = "background-color: #00778b;
            color: #FFFFFF;",
            
            HTML("<br><br>")
          ),
          
          ##############################################
          #Apply secondary footer    
          tags$footer(
          tags$br(),
            
          tags$div(
            id = "secondary-footer-docs"
          ),
          
          tags$br(), 
            
          style = "background-color: #253746;
          color: #FFFFFF;"
          )
        )
      ),
      
      ##############################################
      #All loci bar plot page        
      tabItem(tabName = "barplot_all",
              
                fluidRow(
                  box(width = 6,
                    div(style = "display: inline-block; padding: 0 10px 0 0", numericInput("vfullplot", "Vertical Plot Size (px)", 5000, 500, max = NA, step = 500)),
                    div(style = "display: inline-block; padding: 0 10px 0 10px", numericInput("ncolumns", "Number of Columns", 4, 1, 6)),
                    div(style = "display: inline-block; padding: 0 0 0 10px", numericInput("rdt", "Temporary Read Depth Filter", 5, 5, max = NA, step = 5))
                  )
                ),
                
                fluidRow(
                  uiOutput("plot.ui")
                )
      ),

      ##############################################
      #Final allele calls page                
      tabItem(tabName = "final_df",
              
                fluidRow(
                  box(width = 6,
                      div(style = "display: inline-block; padding: 0 10px 0 0", numericInput("vfinalplot", "Vertical Plot Size (px)", 5000, 500, max = NA, step = 500)),
                      div(style = "display: inline-block; padding: 0 10px 0 10px", numericInput("ncolumnsFinal", "Number of Columns", 4, 1, 6)),
                      div(style = "display: inline-block; padding: 0 0 0 10px", numericInput("rdtFinal", "Temporary Read Depth Filter", 10, 0, max = NA, step = 5))),
                  
                  column(width = 3,
                        align = "center",
                      downloadLink("saveDT", 
                                   label = h4("Save Final Data Frame", 
                                              style = "font-family: Times New Roman"), 
                                   style = 'padding: 10px'),
                      downloadLink("saveOGDF", 
                                   label = h4("Save Pre-Interp Data Frame", 
                                              style = "font-family: Times New Roman"), 
                                   style = 'padding: 10px')),
                  
                  column(width = 3,
                        align = "center",
                      downloadLink("saveAlleleSummary", 
                                   label = h4("Save LB-Allele Summary", 
                                              style = "font-family: Times New Roman"), 
                                   style = 'padding: 10px'),
                      downloadLink("saveLocusSummary", 
                                   label = h4("Save Locus Summary", 
                                              style = "font-family: Times New Roman"), 
                                   style = 'padding: 10px'))
                ),
              
                fluidRow(
                  uiOutput("finalPlot.ui")
                ),
              
                fluidRow(
                  style = "padding: 50px 0px 0px 0px;",
                  box(width = 12,
                    rHandsontableOutput("finalTable")
                  )
                )
      ),
      
      
      ##############################################
      #STRidER Upload page                
      # tabItem(tabName = "strider",
      #         fluidRow(
      #           box(width = 4,
      #               div(style = "display: inline-block; padding: 0 0 0 10px", textInput("popdesc", "Population Description*")),
      #               tags$p("* Must contain a description of population(s) reported (e.g., the title of the study), number of samples, geographic origin, and the number of STR loci."),
      #               tags$hr(class = "thiccHr"),
      #               
      #               div(style = "display: inline-block; padding: 0 0 0 10px", textInput("metadata", "MetaPopulation Info [Optional]**")),
      #               tags$p("** Comments or description of the detailed geographic background and the appropriate metapopulation affiliation of the genotypes."),
      #               tags$hr(class = "thiccHr"),
      #               div(style = "display: inline-block; padding: 0 0 0 10px", textInput("name", "Author's Name")),
      #               
      #               div(style = "display: inline-block; padding: 0 0 0 10px", textInput("email", "Author's Email"))
      #           ),
      #           box(width = 4,
      #               div(
      #                 style = "display: inline-block; padding: 0 0 0 0px", actionButton("strider_ext", "Upload Data Set***")
      #               ),
      #               tags$p("*** Data set needs to be a long-form data table. See manual for more details.")
      #           )
      #         ),
      #         
      #         column(2,
      #                align = "center",
      #                offset = 5,
      #                div(
      #                  style = "display: inline-block; padding: 0 0 0 0px", downloadLink("saveSTRidER", label = h4("Save STRidER Data File", style = "font-family: Times New Roman"), style = 'padding: 10px')
      #                )
      #         ),
      #         
      #         fluidRow(
      #           box(width = 12,
      #               tableOutput("strider_table")
      #           )
      #         )
      #         
      # ),
      
      
      ##############################################
      #STRait Razor Analysis page          
      tabItem(tabName = "manual",
           uiOutput("sraManual")   
      ),
      
      ##############################################
      #Contact page          
      tabItem(tabName = "email",
        fluidRow(
          box(width = 6,
            h1("Support"),
            tags$hr(),
            tags$p("If you have any issues, please feel free to contact me."),
            tags$br(),
                    
            tags$p(
              tags$strong("Jonathan King"),
              tags$br(), "Laboratory Manager",
              tags$br(), "817-735-2940",
              tags$br(), tags$a(href = "mailto: jonathan.king@unthsc.edu", "jonathan.king@unthsc.edu")
            ),
                    
            tags$p(
              tags$strong("Center for Human Identification"),
              tags$br(), "3500 Camp Bowie Blvd.",
              tags$br(), "Fort Worth, Tx 76107"
            )
          )
        ),
        HTML("<br><br><br><br><br><br><br><br><br><br><br><br><br><br><br><br>"),
        tags$div(
          id = "footers",
          ##############################################
          #Apply primary footer
          
          tags$footer(
            
            tags$a(
              href = "https://www.unthsc.edu/graduate-school-of-biomedical-sciences/laboratory-faculty-and-staff/",
              target = "_blank",
              tags$img(
                src = "UNTCHI_RDU.png",
                width = "75",
                height = "75",
                align = "left",
                style = "margin: 10px 10px"
              )
            ), 
            
            tags$a(
              href = "https://www.untchi.org/",
              target = "_blank",
              tags$img(
                src = "CHI_2018.png",
                width = "75",
                height = "75",
                align = "right",
                style = "margin: 10px 10px"
              )
            ),
            
            tags$div(
              id = "social-media-icons",
              
              tags$br(),
              
              tags$a(
                href ="https://twitter.com/UNTCHI_RDU",
                target = "_blank",
                tags$img(src = "twitter.png",
                         title = "UNTCHI RDU on Twitter",
                         width = "30",
                         height = "30",
                         class = "displayed",
                         style = "margin: 10px 10px"
                )
              ),
              
              tags$a(
                href ="https://www.instagram.com/untchi_rdu/",
                target = "_blank",
                tags$img(src = "instagram.png",
                         title = "UNTCHI RDU on Instagram",
                         width = "30",
                         height = "30",
                         class = "displayed",
                         style = "margin: 10px 10px"
                )
              ),
              
              tags$a(
                href = "https://www.facebook.com/UNTCenterForHumanIdentification",
                target = "_blank",
                tags$img(src = "facebook.png",
                         title = "UNTCHI RDU on Facebook",
                         width = "30",
                         height = "30",
                         class = "displayed",
                         style = "margin: 10px 10px"
                )
              )
            ),
            
            style = "background-color: #00778b;
            color: #FFFFFF;",
            
            HTML("<br><br>")
          ),
          
          ##############################################
          #Apply secondary footer    
          tags$footer(
            tags$br(),
            
            tags$div(
              id = "secondary-footer-docs"
            ),
            
            tags$br(), 
            
            style = "background-color: #253746;
            color: #FFFFFF;"
          )
        )
      ),
      
      
      ##############################################
      #Settings page
      tabItem(tabName = "settings",
              fluidRow(
                box(width = 12,
                    
                    div(style = "display: inline-block; padding: 0 10px 0 0", numericInput("cores", "Number of cores for v3", 6, 2, max = NA, step = 2)),
                    div(style = "display: inline-block; padding: 0 10px 0 10px", numericInput("gRDT", "Global Read Depth Thresh", 10, 1, NA)),
                    div(style = "display: inline-block; padding: 0 0 0 10px", checkboxInput("recur", "Recursive Directories", value = TRUE)),
                    div(style = "display: inline-block; padding: 0 0 0 10px", checkboxInput("expectedLoci", "Expected Loci", value = TRUE)),
                    div(style = "display: inline-block; padding: 0 0 0 0", checkboxInput("expertCall", "AutoPass Loci", value = FALSE))
                )
              )
      )
      
    )        
  )
)


##############################################
#Define server rules
server <- function(input, output, session) {
  
  options(shiny.maxRequestSize = 1000*1024^2)

  ##############################################
  #Reactive values start-up
  
  #Initialize final data frame for called alleles
  final <- reactiveValues(df = NULL)
  
  #Initialize default locus and locus nav buttons
  values <- reactiveValues(lcount = 0, fti = -1, fileType = 0)
  
  #Initialize db
  db <- reactiveValues()
  
  ##############################################
  #Navigate tabs
  observeEvent(input$toPlot, {updateTabItems(session, "tabs", "barplot_all")})


  ##############################################
  #Read in database files
  sr_db <- eventReactive(input$loaddb, {

    db$dbPath <- read_csv("db/db/DatabasePath.csv")
    db$dipdb<- read_csv("db/db/DiplotypeDatabase.csv")
    db$hapdb <- read_csv("db/db/HaplotypeDatabase.csv")
    db$sraConfig <- read_csv("db/db/STRaitRazorAnalysisConfig.csv")
    db$sr8rzrConfig <- read_csv("db/db/STRaitRazorConfig.csv")
    
  })
  
  
  ##############################################
  #Compile list of kits   
  output$kits <- renderUI({
    
    shiny::req(sr_db())
    
    STRaitRazorAnalysisConfig <- as.data.frame(db$sraConfig)
    kits1 <-sort(unique(STRaitRazorAnalysisConfig[,1]))
    selectInput("kits", "Select Kit", kits1, selected = 1)
    
  })  
  
  
  ##############################################
  #Isolate kit from config
  rt <- eventReactive({
    input$expertCall
    values$df
    1
    }, {
      
    shiny::req(DF1())

    if(input$expertCall){
      DF <- DF1()
      DF %>%
        select(Status, locStatus, Locus) %>%
        distinct(Locus, .keep_all = TRUE)
    }else{
      STRaitRazorAnalysisConfig <- as.data.frame(db$sraConfig) %>%
        filter(Kitid == input$kits)
    }
    
  })
  
  
  ##############################################
  #Compile list of loci for drop-down menu             
  output$loci <- renderUI({
    
    if(input$cLocus){
      selectInput("loci", label = h4("Choose Locus", style = "font-family: Times New Roman"), choices = rt()$Locus)
    }else{}
    
  })
  
  
  ##############################################
  #Move to previous locus
  observeEvent(input$backone, {
    
                maxLcount <- length(rt()$Locus)
                if(values$lcount == 1){
                   values$lcount = maxLcount
                }else{
                   values$lcount <- values$lcount - 1
                }
                
  })
  
  #Move to next locus and, if selected, save data
  observeEvent(input$upone, {

    if(input$autoSave == TRUE){
      
      #Assign interactive table to active data frame
      DF <- hot_to_r(input$hot)
      
      #Select appropriate columns
      newLines <- as.data.frame(select(filter(DF, Status == TRUE), -Status, -HapLen, -ID)) %>% 
        mutate(SLH = paste0(SampleName, "_", Locus, "_", Haplotype))
      
      #isolate unique entries and combine with current data frame
      newLines <- if(is.null(final$df)){
        final$df <- newLines
      }else{
        newLines <- anti_join(newLines, final$df, by = "SLH")
        final$df <- rbind(final$df, newLines)
      }
      
      #Move to next locus
      maxLcount <- length(rt()$Locus)
      if(values$lcount == maxLcount){
        values$lcount = 1
      }else{
        values$lcount <- values$lcount + 1
      }
      
    }else{
      
      maxLcount <- length(rt()$Locus)
      if(values$lcount == maxLcount){
        values$lcount = 1
      }else{
        values$lcount <- values$lcount + 1
      }
      
    }
                
  })
  
  #Return to first locus
  observeEvent({
    input$resetlocus
    values$fti
    1
    }, values$lcount <- 1)
  
  
  ##############################################
  #Clear final data frame
  observeEvent({
    input$expertCall
    input$clearFinal
    values$fti
    1
    }, {
    
    final$df <- NULL
    values$lcount <- 1
    
  })
  
  
  ##############################################
  #Clear final data frame
  observeEvent({
    values$fti
    input$expectedLoci
    1
  }, {
    
    if(is.null(db$sraConfig)){return()} 
    if(input$expectedLoci){
      
      db$loci <- as.data.frame(db$sraConfig) %>% 
        filter(Kitid == input$kits) %>% 
        select(Locus)
      
    }else{
      
      loci <- as.data.frame(values$df) %>% 
        select(Locus)
      db$loci <- unique(loci)  
      
    }
    
  })
  
  
  ##############################################
  #Display profile progress
  output$progressBox <- renderValueBox({
    
    if (values$fileType <= 0){return(
      
      valueBox(
        paste0(round(0 * 100, digits = 1), " %"), "Progress", icon = NULL,
        color = "black"
      )
      
    )}
    
    if(is.null(final$df)){return(
      
      valueBox(
        paste0(round(0 * 100, digits = 1), " %"), "Progress", icon = NULL,
        color = "black"
      )
      
    )}
    
    maxLcount <- nrow(unique(db$loci))
    
    valueBox(
      paste0(round((length(unique(final$df$Locus))/maxLcount)*100, digits = 1), " %"), "Progress", icon = NULL,
      color = "black"
    )
  })
  
  
  ##############################################
  #Display passing loci count
  output$primaryBoxCount <- renderValueBox({
    if (values$fileType <= 0){return(
      valueBox(
        round(0, digits = 0), "Pass QC", icon = icon("thumbs-up", class = "small_icon_test")
      )
    )}
    
    if(input$expertCall){
      DF <- final$df %>% 
        filter(locStatus == "primary")
    }else{
      DF <- DF1() %>% 
        filter(locStatus == "primary")
    }
    
    valueBox(
      length(unique(DF$Locus)), "Pass QC", icon = icon("thumbs-up", class = "small_icon_test")
    )
  })
  
  
  ##############################################
  #Display warning loci count
  output$warningBoxCount <- renderValueBox({
    if (values$fileType <= 0){return(
      valueBox(
        round(0, digits = 0), "Flagged", icon = icon("exclamation", class = "small_icon_test")
      )
    )}
    
    if(input$expertCall){
      DF <- DF1()
    }else{
      DF <- DF1() %>% 
        filter(locStatus == "warning")
    }
    valueBox(
      length(unique(DF$Locus)), "Flagged", icon = icon("exclamation", class = "small_icon_test")
    )
  })
  
  
  ##############################################
  #Pass selected locus to UI
  output$TargetLocus <- renderPrint({
    
    if(input$cLocus){
      HTML(cat("<font>", input$loci), "</font>")
    }else{
      HTML(cat("<font>", rt()[values$lcount, 3]), "</font>")
    }
  })
  
  
  #############################################
  #Define data frame based on the file type selected
  
  observeEvent({
    values$fti
    }, {
      
    shiny::req(values$fileType)
      
    id <- showNotification("Processing data...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)  
      
    #Choose fastq file, run v3 in shell, and process DF
    fastq1 <- {
      if(values$fileType == 1000){
      sr_db <- as.data.frame(db$dbPath)
      
      GlobalRDT <- input$gRDT
      coreUsage <- input$cores
      
      fastqUpload <- input$fastq
      
      SampleName <- if(input$sampleID == ""){
        fastqUpload$name
      } else{
        input$sampleID
      }
      
      allsequences <- STRaitRazorv3(fastqUpload$datapath, SampleName, numcores = coreUsage, Kit = input$kits, KitPath = file.path("db", "db", sr_db[[9, 2]]))
      
      SampleName <- unlist(allsequences[1])
      STRaitRazorIO <- data.frame(allsequences[2])
      
      FilterResults <- STRaitRazorFilter(allsequences = STRaitRazorIO, SampleName = SampleName, Kit = input$kits, GlobalRDT = GlobalRDT, mkallseq_version = "v3", DBPath = file.path("db", "db", sr_db[[10, 2]]))
      
      STRaitRazorIO <- data.frame(FilterResults[1]) %>% 
        mutate(SampleName = SampleName) %>% 
        select(SampleName, everything())
      
      list(STRaitRazorIO, data.frame(FilterResults[2]), data.frame(FilterResults[3]))

      }else{list(data.frame(1:13), data.frame(1:5), data.frame(1:4))}
      
    }
    
    
    #Choose allsequences file and process DF
    allseq1 <- {

      if(values$fileType == 100){
        sr_db <- as.data.frame(db$dbPath)
        GlobalRDT <- input$gRDT

        allseqUpload <- input$allseq
        SampleName <- input$sampleID
        
        STRaitRazorIO <- read_delim(allseqUpload$datapath, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
        
        FilterResults <- STRaitRazorFilter(allsequences = STRaitRazorIO, SampleName = SampleName, Kit = input$kits, GlobalRDT = GlobalRDT, mkallseq_version = input$versionS, DBPath = file.path("db", "db", sr_db[[10, 2]]))
        
        STRaitRazorIO <- data.frame(FilterResults[1]) %>%
          mutate(SampleName = SampleName) %>%
          select(SampleName, everything())

        list(STRaitRazorIO, data.frame(FilterResults[2]), data.frame(FilterResults[3]))
        
        
      }else{list(data.frame(1:13), data.frame(1:5), data.frame(1:4))}
    }
    
    
    #Prepare the appropriate data set
    values$df <- if(values$fileType == 1000){
          fastq1[[1]]} else{
          if(values$fileType == 100){
            allseq1[[1]]} else{"Error"}}
    
    values$allelesummary <- if(values$fileType == 1000){
      fastq1[[2]]} else{
        if(values$fileType == 100){
          allseq1[[2]]} else{"Error"}}
    
    values$locussummary <- if(values$fileType == 1000){
      fastq1[[3]]} else{
        if(values$fileType == 100){
          allseq1[[3]]} else{"Error"}}

  })
  
  
  #############################################
  #Update reactive value to define which file type has been selected
  observeEvent(input$fastq, {
    values$fileType <- 1000
  })
  
  observeEvent(input$allseq, {
    values$fileType <- 100
  })
  
  observeEvent({
    input$allseq
    input$fastq
    1
  },{
    values$fti = 0
    values$fti <- values$fti + 1
  })
  

  ##############################################
  #Clean data and check auto-called alleles 
  DF1 <- eventReactive({
    input$expertCall
    values$df
    1
    }, {

    #Ready config
    sr_db <- as.data.frame(db$sraConfig)
    
    Config <- sr_db %>% 
      filter(Kitid == input$kits) %>% 
      select(Locus, AutoRAPT)

    #Ready the data
    DF <- values$df

    DF <- DF %>% 
      mutate(Status = ifelse(AutoCall == "a", TRUE, FALSE)) %>% 
      select(-AutoCall, -AR) %>% 
      select(Status, everything()) %>% 
      mutate(HapLen = nchar(DF$Haplotype)) %>% 
      filter(Locus != "DYS389I" | Locus == "DYS389I" & Allele <= 25)
    
    #Add threshold for combined RAP to data frame
    DF <- left_join(DF, Config, by = "Locus")
    
    #create identifier and combine called allele RAPs
    DF <- DF %>% 
      mutate(ID = paste0(SampleName, "_", Locus)) 
    
    tempRAP <- summarize(dplyr::group_by(filter(DF, Status == TRUE), ID), RAP_Sum = sum(RAP))

    #Add combined RAP to data frame, classify locus as passing QC, and remove temporary vectors
    DF <- left_join(DF, tempRAP, by = "ID") %>%
      ungroup() %>% 
      mutate(locStatus = if_else(RAP_Sum < AutoRAPT, "warning", "primary")) %>% 
      select(SampleName, Locus, FinalAllele, Nomenclature, HaplotypeSum, HB, RAP, SB, Haplotype, LocusRank, locStatus, Allele, HapLen, Status, ID) 
    
    #Conditionally parse data frame by locus status
    if(input$expertCall){

      #Fast-pass the loci passing QC to final data frame
      final$df <- DF %>% 
        filter(locStatus == "primary", Status == "TRUE") %>% 
        select(-HapLen, -Status, -ID) %>% 
        mutate(SLH = paste0(SampleName, "_", Locus, "_", Haplotype))
      
      #Send the rest of the loci to the UI
      DF %>% 
        filter(locStatus == "warning")
    }else{
      DF
    }
    
  })
  
  
  ##############################################
  #Filter data frame for selected locus
  DF2 <- eventReactive({
    input$loci
    values$lcount
    values$df
    1
    }, {
    
    shiny::req(rt())

    if(input$cLocus){
      DF1() %>% 
        filter(Locus == input$loci)
    }else{
      DF1() %>% 
        filter(Locus == rt()[values$lcount, 3])
    } 
    
    
  })
  
  
  ##############################################
  #Update haplotype-depth initial value dynamically with locus-specific value
  observeEvent({
    input$loci
    values$lcount
    1
  }, {
    
    shiny::req(sr_db())
    
    sr_db <- as.data.frame(db$sraConfig)
    
    LocusInfo <- if(input$cLocus){
      sr_db %>% 
        filter(Kitid == input$kits, Locus == input$loci)
    }else{
      sr_db %>% 
        filter(Kitid == input$kits, Locus == rt()[values$lcount, 3])
    } 
    
    updateNumericInput(session, "HD_Threshold", value = LocusInfo[[1, 18]])
  })
  
  
  ##############################################
  #Update strand-balance initial value dynamically with locus-specific value
  observeEvent({
    input$loci
    values$lcount
    1
  }, {
    
    shiny::req(sr_db())
    
    sr_db <- as.data.frame(db$sraConfig)
    
    LocusInfo <- if(input$cLocus){
      sr_db %>% 
        filter(Kitid == input$kits, Locus == input$loci)
    }else{
      sr_db %>% 
        filter(Kitid == input$kits, Locus == rt()[values$lcount, 3])
    } 
    
    updateNumericInput(session, "SB_Threshold", value = LocusInfo[[1, 19]])
  })  
  
  
  ##############################################
  #Update relative allele proportion initial value dynamically with locus-specific value
  observeEvent({
    input$loci
    values$lcount
    1
  }, {
    
    shiny::req(sr_db())
    
    sr_db <- as.data.frame(db$sraConfig)
    
    LocusInfo <- if(input$cLocus){
      sr_db %>% 
        filter(Kitid == input$kits, Locus == input$loci)
    }else{
      sr_db %>% 
        filter(Kitid == input$kits, Locus == rt()[values$lcount, 3])
    } 
    
    updateNumericInput(session, "RAP_Threshold", value = LocusInfo[[1, 20]])
  })  
    
  ##############################################
  #Render box UI and ggplot
  
  #Create box and insert ggplot
  output$plotFullUI <- renderUI({

    if ((values$fileType <= 0)) {return(invisible())}
    if ((nrow(as_data_frame(DF2())) < 1)) {return(
      box(title = h4(uiOutput('TargetLocus'), style = "font-family: Times New Roman"), solidHeader = TRUE, status = "danger", plotOutput("plotFull"))
    )}
    
    status <- DF2()[1, 11]

    box(title = h4(uiOutput('TargetLocus'), style = "font-family: Times New Roman"), solidHeader = TRUE, status = status, plotOutput("plotFull"))
    
  })
  
  #Create ggplot
  output$plotFull <- renderPlot({

    #Apply thresholds
    DF <- DF2() %>%
      filter(SB >= input$SB_Threshold, HaplotypeSum >= input$HD_Threshold, RAP >= input$RAP_Threshold) 
    
    #Axis labels of bar plots
    HapLen <- select(DF, Allele, HapLen) %>% 
      unique() %>% 
      as.data.frame()
    
    #Name columns
    names(HapLen) <- c("Allele", "HapLen")

    #Plot all alleles
    ggplot(DF, aes(HapLen, HaplotypeSum, fill = as.factor(LocusRank))) + geom_bar(stat = "identity", show.legend = FALSE) + theme_classic() + scale_fill_manual(values = Palette272) + scale_y_continuous(expand = c(-0.000005, 0)) + scale_x_continuous(breaks = HapLen$HapLen,labels = HapLen$Allele) + labs(y = "Read Depth", x = "Alleles")
    
  })
  
  
  ##############################################
  #Output table
  output$hot <- renderRHandsontable({
    shiny::req(DF2())

    DF <- DF2() %>% 
      relocate(Status, .before = SampleName)
    
    if (!is.null(DF))
      
      rhandsontable(filter(DF,SB >= input$SB_Threshold, HaplotypeSum >= input$HD_Threshold, RAP >= input$RAP_Threshold), useTypes = TRUE, stretchH = "all")
    
  })
  
  
  ##############################################
  #ggplot Final
  output$plotCalled <- renderPlot({
    
    shiny::req(input$hot)
    
    #Assign interactive table to active data frame and apply thresholds
    DF <- hot_to_r(input$hot) %>% 
      filter(SB >= input$SB_Threshold, HaplotypeSum >= input$HD_Threshold, RAP >= input$RAP_Threshold) 
    
    #Grab appropriate sample names and filter called alleles
    HapLen <- select(DF, Allele, HapLen, Status) %>% 
      unique() %>% 
      as.data.frame()
    
    #Name columns
    names(HapLen) <- c("Allele", "HapLen", "Status")
    
    #Filter DF for only called alleles
    HapLen <- filter(HapLen, Status == TRUE)
    
    #Plot final allele calls
    ggplot(filter(DF, Status == TRUE), aes(HapLen, HaplotypeSum, fill = as.factor(LocusRank))) + geom_bar(stat = "identity", show.legend = FALSE) + theme_classic() + scale_fill_manual(values = Palette272) + scale_y_continuous(expand = c(-0.000005, 0)) + scale_x_continuous(breaks = HapLen$HapLen, labels = HapLen$Allele) + labs(y = "Read Depth", x = "Alleles")
    
  })
  
  
  ##############################################
  #Save locus to  final data frame
  observeEvent(input$save2df, {
    
    #Assign interactive table to active data frame
    DF <- hot_to_r(input$hot)
    
    #Select appropriate columns
    newLines <- as.data.frame(select(filter(DF, Status == TRUE), -Status, -HapLen, -ID)) %>% 
      mutate(SLH = paste0(SampleName, "_", Locus, "_", Haplotype))
    
    #isolate unique entries and combine with current data frame
    newLines <- if(is.null(final$df)){
      final$df <- newLines
    }else{
      newLines <- anti_join(newLines, final$df, by = "SLH")
      final$df <- rbind(final$df, newLines)
    }
    
  })
  

  ##############################################
  #Plot all STR loci 
  output$allLociplot <- renderPlot({

    #Load the data frame  
    DF <- values$df
    
    #Add length for plotting
    DF <- DF %>% 
      mutate(Status = ifelse(AutoCall == "a", TRUE, FALSE)) %>%
      mutate(HapLen = nchar(DF$Haplotype)) %>% 
      filter(Locus != "DYS389I" | Locus == "DYS389I" & Allele <= 25)
    
    #Extract name
    SampleName <- DF[1, 1]
    
    #Load config file

    STRaitRazorAnalysisConfig <- as.data.frame(db$sraConfig) %>% 
      filter(Kitid == input$kits)
    
    Type <- select(STRaitRazorAnalysisConfig, Locus, Marker_Type)
      
    #Enrich for STRs and filter out noise
    DF <- left_join(DF, Type) %>% 
      filter(Marker_Type %in% c("STR","INDEL")) %>% 
      filter(HaplotypeSum >= input$rdt)
      
    ggplot(DF, aes(HapLen, HaplotypeSum, fill = as.factor(LocusRank))) + geom_bar(stat = "identity", show.legend = FALSE) + facet_wrap(~Locus, ncol = input$ncolumns, scales = "free_x") + theme_classic() + scale_fill_manual(values = Palette272) + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(face = "bold", size = 12, lineheight = 5.0, color = "#FFFFFF"), strip.background = element_rect(fill = "#00778b"), axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_text(vjust = 0.5, hjust = 0.5)) + labs(title = paste0(SampleName, " Sequence-based STR Barplot"), y = "Read Depth", x = "Alleles") + geom_text(aes(x = HapLen, y = (0-(max(HaplotypeSum)/10)), label = round(Allele, 1), angle = 90, hjust = 1, vjust = 0), size = 5) + geom_hline(yintercept = 0) + scale_y_continuous(expand = c(0.15,0))
      
  })
  

  ##############################################
  #Push plot to UI
  output$plot.ui <- renderUI({
    plotOutput("allLociplot", height = input$vfullplot)
  })
  
  
  ##############################################
  #Plot all STR loci 
  output$finalLociplot <- renderPlot({
    
    #Import data frame
    DF <- final$df %>% 
      mutate(HapLen = nchar(final$df$Haplotype))
    
    SampleName <- DF[1, 1]
    
    #Import analysis config

    STRaitRazorAnalysisConfig <- as.data.frame(db$sraConfig) %>% 
      filter(Kitid == input$kits)
    
    Type <- select(STRaitRazorAnalysisConfig, Locus, Marker_Type)
    
    #Enrich for STRs and filter out noise
    DF <- left_join(DF, Type) %>% 
      filter(Marker_Type %in% c("STR","INDEL")) %>% 
      filter(Locus != "DYS389I" | Locus == "DYS389I" & Allele <= 25) %>% 
      filter(HaplotypeSum >= input$rdtFinal) 
    
    ggplot(DF, aes(HapLen, HaplotypeSum, fill = as.factor(LocusRank))) + geom_bar(stat = "identity", show.legend = FALSE) + facet_wrap(~Locus, ncol = input$ncolumnsFinal, scales = "free_x") + theme_classic() + scale_fill_manual(values = Palette272) + theme(panel.spacing = unit(2, "lines"), strip.text = element_text(face = "bold", size = 12, lineheight = 5.0, color = "#FFFFFF"), strip.background = element_rect(fill = "#00778b"), axis.text.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), axis.title.x = element_text(vjust = 0.5, hjust = 0.5)) + labs(title = paste0(SampleName, " Sequence-based STR Barplot"), y = "Read Depth", x = "Alleles") + geom_text(aes(x = HapLen, y = (0-(max(HaplotypeSum)/10)), label = round(Allele, 1), angle = 90, hjust = 1, vjust = 0), size = 5) + geom_hline(yintercept = 0) + scale_y_continuous(expand = c(0.15,0))
    
  })
  
  
  ##############################################
  #Render final data table  
  output$finalTable <- renderRHandsontable({
    shiny::req(final$df)
    
    rhandsontable(select(final$df, 1:11),
      useTypes = TRUE, 
      stretchH = "all") 
    
  })
  
  
  ##############################################
  #Save final data frame
  
  #Final profile tab
  output$saveDT <- downloadHandler(
    filename = function(){
      paste0(final$df[1, 1], "_FinalDataTable.tsv")
    },
    content = function(file){
      write_tsv(as.data.frame(final$df), file)
    }
  )
  
  #Data Analysis tab
  output$saveDTalt <- downloadHandler(
    filename = function(){
      paste0(final$df[1, 1], "_FinalDataTable.tsv")
    },
    content = function(file){
      write_tsv(as.data.frame(final$df),file)
    }
  )
  
  
  ##############################################
  #Save original data frame after filtering
  output$saveOGDF <- downloadHandler(
    filename = function(){
      paste0(values$df[[1, 1]], "_DataTable.tsv")
    },
    content = function(file){
      write_tsv(as.data.frame(values$df), file)
    }
  )
  
  #Save length-based allele summary
  output$saveAlleleSummary <- downloadHandler(
    filename = function(){
      paste0(values$df[[1, 1]], "_AlleleSummary.tsv")
    },
    content = function(file){
      write_tsv(as.data.frame(values$allelesummary), file)
    }
  )
  
  #Save locus summary
  output$saveLocusSummary <- downloadHandler(
    filename = function(){
      paste0(values$df[[1, 1]], "_LocusSummary.tsv")
    },
    content = function(file){
      write_tsv(as.data.frame(values$locussummary), file)
    }
  )
  
  ##############################################
  #Save locus summary
  # output$STRidER <- downloadHandler(
  #   filename = function(){
  #     paste0(values$df[[1, 1]], "STRidER_Data_File.tsv")
  #   },
  #   content = function(file){
  #     write_tsv(as.data.frame(values$STRidER), file)
  #   }
  # )
  
  ##############################################
  #Push plot to UI
  output$finalPlot.ui <- renderUI({
    shiny::req(input$vfinalplot)
    
    plotOutput("finalLociplot", height = input$vfinalplot)
  })
  
  
  ##############################################
  #Open STRait Razor manual
  output$sraManual <- renderUI({tags$iframe(style = "height: 800px; width: 100%", src = "STRaitRazorOnlineManual.pdf")})
 
  
}

##############################################
# Run the application in R script
shinyApp(ui = ui, server = server)

