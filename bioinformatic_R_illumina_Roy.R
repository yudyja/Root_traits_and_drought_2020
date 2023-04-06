install.packages("dada2") # did not work but the next ,,yes!

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("dada2")


packageVersion("dada2")
packageVersion("ShortRead")
packageVersion("Biostrings")
install.packages("icesTAF")

library(dada2)
library(ShortRead)
library(Biostrings)
library(icesTAF)


# FILE PARSING
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
path<-"D:/YUDI LOZANO BERNAL/BIBS project/GRASSLAND_Sp_EXP/Ilumina/Data grasslands Ilumina/data_prueba"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))
if(length(fnFs) != length(fnRs)) stop("Forward and reverse files do not match.")

####????????????Para q los organiza. 

# IDENTIFY PRIMERS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FWD <- "GTGARTCATCGAATCTTTG" 
REV <- "TCCTCCGCTTATTGATATGC"

# identify whether primers are in the good orientation
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
######???? q es todo esto??

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
REV.orients
#####???como sabes q han quedado orientados???

# pre-filter reads for ambiguous bases (Ns), because Ns could be in primers, 
# precluding identification and removal
# Also, some reads might be in the opposite directions:
# it is needed to reorient these reads in the good forward orientation

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0,
              primer.fwd="GTGA", multithread = T)
# now i have files in the "filtN folder
### q hace este filter and trim..porq solo para el forward..porq no para el reverse

# Count the number of times the primers appear in the forward and reverse read, 
# while considering all possible primer orientations. 
# Just for one sample

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# REMOVE PRIMERS
#~~~~~~~~~~~~~~~~
#you have to use cutadapt. see here https://www.biostars.org/p/85788/ 
#first. download anaconda https://www.anaconda.com/distribution/
#download Microsoft Visual C++ Compiler for Python 2.7  
#https://www.microsoft.com/en-us/download/details.aspx?id=44266
#Me salio error en la instalacion. download visual studio 2017. revisar unos parametros que no recuerdo. 
#open console anaconda and type: pip install cutadapt

# Use cutadapt
cutadapt <- "C:/Users/Yudi/Anaconda3/scripts/cutadapt"

system2(cutadapt, args = "--version") # Run shell commands from R

# Create output filenames for the cutadapt-ed files

# Define the parameters to give the cutadapt command
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt


for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-m", 10, # to avoid to have 0 length sequences
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# Check for primer removal
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# DAD2 pipeline
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Forward and reverse fastq filenames (assuming the name has the following format):
cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][2]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

# FILTER AND TRIM
#~~~~~~~~~~~~~~~~~~~
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(1, 2), 
              truncQ = 2, minLen = 50, rm.phix = TRUE, 
              compress = TRUE, multithread = TRUE)  # on windows, set multithread = FALSE

# LEARN ERROR RATE
#~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nread=1e6, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nread=1e6, multithread=TRUE)

# SAMPLE INFERENCE AND MERGER OF PAIRED-END READS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  sam.tmp <- paste(sam,"_",sep="")
  derepF <- derepFastq(filtFs[[grep(sam.tmp,filtFs)]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[grep(sam.tmp,filtRs)]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

mergers <- vector("list", length(sample.names))
for(i in 1:length(sample.names)) {
  cat("Processing:", i, "\n")
  #sam.tmp <- paste(sam,"_",sep="")
  derepF <- derepFastq(filtFs[[i]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[i]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
names(mergers) <- sample.names
names(mergers) <- NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
rownames(seqtab) <- sample.names

# repeat for run2 and merge the tables
seqtab.merged <- mergeSequenceTables(seqtab.pool1,seqtab.pool2)

# remove chimera
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seqtab.merged.nochim <- removeBimeraDenovo(seqtab.merged, method="consensus", multithread=TRUE, verbose=TRUE)

# assign to fungi
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(100) # Initialize random number generator for reproducibility
tax <- assignTaxonomy(seqtab.merged.nochim, 
                      "UNITE database/sh_general_release_dynamic_01.12.2017.fasta", 
                      minBoot=100, multithread=TRUE)

# subset seqtab and taxtab to fungi
seqtab.merged.nochim.fung <- seqtab.merged.nochim[,!is.na(tax[,1])]
taxtab.merged.nochim.fung <- tax[which(rownames(tax)%in%colnames(seqtab.merged.nochim.fung)),]

# set ESVs names
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SEQ.NAMES <- paste("esv", 1:ncol(seqtab.merged.nochim.fung), sep = "")
outfile <- "~/Yudi/seq.merged.nochim.fung.fasta"
uniquesToFasta(seqtab.merged.nochim.fung, outfile, ids = SEQ.NAMES)
# save sequences and names
esvNames_esvSeqRef_Seq_nochim.fung <- data.frame(esv.names =SEQ.NAMES,
                                                 esv.seqRef = colnames(seqtab.merged.nochim.fung))
colnames(seqtab.merged.nochim.fung) <- SEQ.NAMES
rownames(taxtab.merged.nochim.fung) <- SEQ.NAMES