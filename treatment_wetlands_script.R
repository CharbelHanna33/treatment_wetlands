#############################################################################
######################## 				     ########################
########################  TREATMENT WETLANDS SCRIPT  ########################
######################## 				     ########################
#############################################################################

##### Here is the R script of data analyses for the treatment wetlands project. 

##### The first part of the script deals with the bioinformatics,
##### going through the DADA2 pipeline. Prealably demultiplexed 
##### amplicon data were filtered, trimmed, and merged, chimeras were
##### removed and taxonomy was assigned. Small differences in methodology
##### set apart the 16S (bacterial) analyses from the ITS (fungal) analyses,
##### due to the variable length of ITS (fungal) amplicons.

##### The second part of the script deals with the biostatistics,
##### using phyloseq objects for analyses and plots. The data collected
##### encompasses two different experiments that were completed on the same
##### experimental setup. Therefore, the first step, after creating the 
##### phyloseq objects, will be to subset these objects to retain only the 
##### data used specifically in the biochar project. Then, these data will be
##### rarefied, and multiple statistical analyses will be performed.

##### Numerous RDS objects are generated along the way. Each section
##### of the script has its relevant RDS objects listed in the beginning.
##### They are written in a relative path, relative to your own working
##### directory, set by the function setwd(). The line where these various
##### RDS objects are generated remain in the script, as commentary:
##### saveRDS(object,"./object.RDS"). It is not necessary to run these lines,
##### unless you have changed the parameters to test alternative parameters.

##### I have provided one main folder containing two subfolders: 
##### one for 16S, one for ITS. I recommend you set the working directory 
##### to the main folder, wherever you have placed it on your computer. 
##### I have written the paths accordingly, by specifying
##### the relevant subfolders throughout the script.

setwd(choose.dir()) # find and select the main folder provided
			  # or write the path to the folder in setwd().

##### All required packages are listed in the beginning of the script.

#############################################################################

##### Loading libraries and setting theme.

library(dada2) 
library(R.utils)
library(microseq)
library(phyloseq) 
library(Biostrings) 
library(ggplot2)
library(dplyr)
library(tidyr)
library(ShortRead)
library(ggpubr)
library(microbiome) 
library(microbiomeutilities)  
library(RColorBrewer) 
library(vegan)
library(multcomp)
library(microViz)
library(DT) 
library(data.table)
library(metagMisc)
library(vsn)
library(ade4)
library(ape)
library(forcats)
library(micro4all)
library(lme4)
library(lmerTest)
library(car)
library(GUniFrac)

theme_set(theme_bw())

##### Installing packages
## Some packages above are tricky to install. Here is how to install them.

### Installing DADA2

 if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
 BiocManager::install("dada2", version = "3.18")

### Installing phyloseq

 biocpkgs_needed = c("ggplot2", "phyloseq")
 pkgs_needed = c("PMA","dplyr","ade4","ggrepel","vegan", "tidyverse", "ggtree")
 letsinstall = setdiff(pkgs_needed, installed.packages()) 
 letsinstallb = setdiff(biocpkgs_needed, installed.packages()) 
 if (length(letsinstallb) > 0) {
 if (!requireNamespace("BiocManager", quietly = TRUE))
	install.packages("BiocManager")
 BiocManager::install(letsinstallb, version = "3.18")
 	}

### Installing ShortRead

 if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

 BiocManager::install("ShortRead", version = "3.18")

### Installing microbiome

 if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

 BiocManager::install("microbiome")

### Installing microViz

 if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
 BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)

 remotes::install_github("david-barnett/microViz")

 install.packages("ggtext") # for rotated labels on ord_plot() 
 install.packages("ggraph") # for taxatree_plots()
 install.packages("DT") # for tax_fix_interactive()
 install.packages("corncob") # for beta binomial models in tax_model()

### Installing ANCOMBC

 if (!require("BiocManager", quietly = TRUE))
 install.packages("BiocManager")

 BiocManager::install("ANCOMBC")

### Installing DESeq2

 if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

 BiocManager::install("DESeq2")

### Installing metagMisc

 library(remotes)
 remotes::install_github("vmikk/metagMisc")

### Installing DECIPHER

 if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

 BiocManager::install("DECIPHER")

#############################################################################

#############################################################################
################################### DADA2 ###################################
#############################################################################

#################################### 16S ####################################

##### Here begins the bioinformatics part for the 16S (bacterial) data.
##### Tutorial for 16S DADA2 : https://benjjneb.github.io/dada2/tutorial.html

### Here are the RDS objects that will be generated along the way.
out16s <- readRDS("./mesocosmes_16S/out16s.rds")
errF16s <- readRDS("./mesocosmes_16S/errF16s.rds")
errR16s <- readRDS("./mesocosmes_16S/errR16s.rds")
dadaFs16s <- readRDS("./mesocosmes_16S/dadaFs16s.rds")
dadaRs16s <- readRDS("./mesocosmes_16S/dadaRs16s.rds")
mergers16s <- readRDS("./mesocosmes_16S/mergers16s.rds")
seqtab16s <- readRDS("./mesocosmes_16S/seqtab16s.rds")
nochim16s <- readRDS("./mesocosmes_16S/nochim16s.rds")
nochimseqtab16s <- readRDS("./mesocosmes_16S/nochimseqtab16s.rds")
taxa16s <- readRDS("./mesocosmes_16S/taxa16s.rds")

### We set a path variable to where the fastq files are located.
path <- "C:/Users/BaseSpace/mesocosmes/mesocosmes_16S" # CHANGE THIS
# to the directory containing the unzipped fastq files in your computer.
list.files(path) # Verify the path by checking if the files are listed.

### This function allows us to look at the reads directly. This will
### confirm to us that the fastq files still have primers on, which
### will need to be removed.

### Here are the primers used for these amplifications (V4-V5; 515F-926R).
# Forward primer 515F:  GTGYCAGCMGCCGCGGTAA
# Reverse primer 926R:  CCGYCAATTYMTTTRAGTTT

# Forward read: 
reads.R1 <- readFastq("./mesocosmes_16S/SCB1-0-15-G-16S_S19_L001_R1_001.fastq.gz")
sread(reads.R1)
# [1] GTGTCAGCAGCCGCGGTAATACAGAGGGGGCAAG...CAGGATTAGATACCCTGGTAGTCCACGCCGGAA
# Reverse read: 
reads.R2 <- readFastq("./mesocosmes_16S/SCB1-0-15-G-16S_S19_L001_R2_001.fastq.gz")
sread(reads.R2)
# [1] CCGCTAACGCTCTTGAGTGTTAATCACGCGAC...TCTCAACCCATGCAGCACACAAGGCCGCCCCT
### The primers appear in both the forward and reverse reads.

### Now, we perform some string manipulations to get matched 
### lists of the forward and reverse fastq files.
# Forward fastq filenames have format: SAMPLENAME_R1_001.fastq 
# Reverse fastq filenames have format: SAMPLENAME_R2_001.fastq
fnFs16s <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs16s <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names.16s <- sapply(strsplit(basename(fnFs16s), "_"), `[`, 1)

### The filtered files that will be soon generated will be placed
### in a /filtered subdirectory.
filtFs16s <- file.path(path, "filtered", paste0(sample.names.16s, 
	"_F_filt.fastq.gz"))
filtRs16s <- file.path(path, "filtered", paste0(sample.names.16s, 
	"_R_filt.fastq.gz"))
names(filtFs16s) <- sample.names.16s
names(filtRs16s) <- sample.names.16s

### We check the quality profile of the forward and reverse reads.
### This will inform our trimming.

# Forward reads:
plotQualityProfile(fnFs16s[29:30]) # Or any other sample in []
# Reverse reads:
plotQualityProfile(fnRs16s[29:30]) # Or any other sample in []
# Weirdly enough, this function sometimes bugs on my end. If the function
# does not work, try restarting R. That often does the trick.

### We now trim the reads. 

### The argument truncLen is concerned with trimming
### at the distal end of the forward and reverse reads. We want to truncate
### at where the quality of the reads drop. On the other hand, it is
### imperative to ensure that our forward and reverse reads are long enough
### to be merged. 

### In our case, we are working here with reads from the 515F/926R region, 
### which has a length of  926-515 = 511 nucleotides. The Illumina MiSeq 
### technology, used for amplification, reads to up to 250 nucleotides,
### hence why we read in forward and reverse. It goes like this:
# 515 ------------------->
#		   <------------------- 926
### We must truncate the distal end (at the head of the arrow, going both ways)
### at the position where the quality drops, but we must retain enough
### nucleotides to allow the forward and reverse reads to merge later.
# 515 ----------------
# 			 ---------------- 926
### We must keep an overlap of at least 20 nucleotides.
### Based on the quality profiles and the above considerations, 
### we will keep the first 250 nucleotides in the forward reads,
### and 190 nucleotides in the reverse reads. This will allow an overlap
### of 29 nucleotides, which will be enough for merging.

### The argument trimLeft is very important here: it removes the primers.
### The forward primer has 19 nucleotides; the reverse primer has 20.

### The other argument are standard parameters. We can leave them as they are.

# This function takes some time, but you do not need to run it:
# its output is already saved as an RDS file.
# The same is true for all other big functions throughout the DADA2 pipeline.
out16s <- filterAndTrim(fnFs16s, filtFs16s, fnRs16s, filtRs16s, 
		truncLen = c(250,190), trimLeft = c(19,20), maxN = 0, 
		maxEE = c(2,2), truncQ = 2, rm.phix = TRUE, compress = TRUE, 
		multithread = FALSE) # On Mac, set to number of cores.
# saveRDS(out16s, "./mesocosmes_16S/out16s.rds")
# We can see how many reads have been removed after trimming and filtering.
head(out16s)

### The files are now dereplicated: identical sequences are combined together.
### It seems that this step is no longer required, as DADA2 now performs it
### automically; but it does not hurt to do it.

derepFs16s <- derepFastq(filtFs16s, verbose = TRUE)
derepRs16s <- derepFastq(filtRs16s, verbose = TRUE)

### The algorithm now must learn the error rates associated with our dataset.

errF16s <- learnErrors(derepFs16s, multithread = FALSE)
# saveRDS(errF16s, "./mesocosmes_16S/errF16s.rds")
errR16s <- learnErrors(derepRs16s, multithread = FALSE)
# saveRDS(errR16s, "./mesocosmes_16S/errR16s.rds")

plotErrors(errF16s, nominalQ = TRUE)
plotErrors(errR16s, nominalQ = TRUE)
# The estimated error rates (black line) fit the observed rates (points),
# which is a good sign that the model works for our data.

### Now, the algorithm will infer the number of amplicon sequence variants
### (ASVs) from the number of unique sequences in each read, based on the 
### estimated error rates. 

# This function allows multithreading on Windows.
dadaFs16s <- dada(derepFs16s, err=errF16s, multithread = 4)
# saveRDS(dadaFs16s, "./mesocosmes_16S/dadaFs16s.rds")
dadaRs16s <- dada(derepRs16s, err=errR16s, multithread = 4)
# saveRDS(dadaRs16s, "./mesocosmes_16S/dadaRs16s.rds")
# We can inspect the dada-class object for the first sample.
# dadaFs16s[[1]]
# 554 sequence variants were inferred from 12571 input unique sequences.

### We now merge the forward and reverse sequences.

mergers16s <- mergePairs(dadaFs16s, derepFs16s, dadaRs16s, derepRs16s,
	 verbose = TRUE)
# saveRDS(mergers16s, "./mesocosmes_16S/mergers16s.rds")
# We can inspect the merger data.frame from the first sample.
head(mergers16s[[1]])

### We now construct a sequence table.

seqtab16s <- makeSequenceTable(mergers16s)
# saveRDS(seqtab16s, "./mesocosmes_16S/seqtab16s.rds")
dim(seqtab16s) # 82 15016 --> 15016 ASVs in 82 samples.
# Inspect distribution of sequence lengths.
table(nchar(getSequences(seqtab16s))) 
# Most sequences have between 368-376 nucleotides.

### We now remove chimeras.

nochim16s <- removeBimeraDenovo(mergers16s, method = "consensus", 
	multithread = FALSE, verbose = TRUE)
# saveRDS(nochim16s, "./mesocosmes_16S/nochim16s.rds")

nochimseqtab16s <- makeSequenceTable(nochim16s)
# saveRDS(nochimseqtab16s, "./mesocosmes_16S/nochimseqtab16s.rds")
dim(nochimseqtab16s) # 82 7299 --> 7299 ASVs in 82 samples
# Inspect proportion of reads kept after removal of chimeras
sum(nochimseqtab16s)/sum(seqtab16s)
# 87% of reads were preserved, so 13% of reads were removed as chimeras.

### We now do a sanity check, by inspecting the number of reads preserved
### at each step of the process.
getN <- function(x) sum(getUniques(x))
track16s <- cbind(out16s, sapply(dadaFs16s, getN), sapply(dadaRs16s, getN), 
	sapply(mergers16s, getN), rowSums(nochimseqtab16s))
# If processing a single sample, remove the sapply calls: 
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track16s) <- c("input", "filtered", "denoisedF", 
	"denoisedR", "merged", "nonchim")
rownames(track16s) <- sample.names.16s
head(track16s)
#                            input filtered denoisedF denoisedR merged nonchim
# BC-0-15-G-16S              28811    22121     20482     21182  10977   10484
# BC-15-30-G-16S             19148    14114     12859     13500   7176    6756
# ctrl-PCR-neg-ChagnonP7-16S    44        1         1         1      1       1
# ctrl-PCR-pos-ChagnonP7-16S 24523    19598     19286     19552  17892   10015
# GC-0-15-G-16S              20592    15720     14477     15056   8160    7684
# GC-15-30-G-16S             14250    11006      9999     10386   5715    5583

# As we can see, we managed to retain a large number of reads, which is
# reassuring. We did not lose almost all reads at the merging step, which
# indicates that the merging was succesful; the length at which reads were
# truncated during the trimming process (truncLen argument) was adequate.

### We now assign taxonomy, by contrasting the sequences with a reference
### database. Here, we use Silva. This database must normally be downloaded 
### by the user and placed in the directory with the fastq files. You will 
### find it provided in the folder, with the fastq and RDS files, so you
### do not need to download it.

# This function allows multithreading on Windows. It is a very large function,
# that will likely not work at all without multithreading. If you have
# access to a service providing supercomputers, you can use it here, and
# increase the number of cores in the multithread argument up to 20. The
# function should still work on a regular computer, if all 4 cores are used.
taxa16s <- assignTaxonomy(nochimseqtab16s, 
	"./mesocosmes_16S/silva_nr99_v138.1_train_set.fa.gz",
	multithread = 4)
# saveRDS(taxa16s, "./mesocosmes_16S/taxa16s.rds")

### We can now inspect the taxonomic assignments.

taxa16s.print <- taxa16s # Removing sequence rownames for display only
rownames(taxa16s.print) <- NULL
taxa16s.print[1:10,]
dim(taxa16s.print) # 7299    6
# We have 7299 ASVs, most of which seem to have been appropriately assigned.

##### This concludes the DADA2 pipeline for the 16S (bacterial) data.

#################################### ITS ####################################

##### Here begins the bioinformatics part for the ITS (fungal) data.
##### Tutorial for ITS DADA2 : https://benjjneb.github.io/dada2/ITS_workflow.html

### The pipeline is mostly the same, although there are differences in the
### first steps. Because the length of the ITS region is variable, we cannot
### truncate the reads at a fixed length. Also, the primers are often found
### on both ends of both the forward and reverse reads. We need to install a
### specific program, called Cutadapt, to remove the primers on all ends. 

### Here are the RDS objects that will be generated along the way.

outITS <- readRDS("./mesocosmes_ITS/outITS.rds")
errFITS <- readRDS("./mesocosmes_ITS/errFITS.rds")
errRITS <- readRDS("./mesocosmes_ITS/errRITS.rds")
dadaFsITS <- readRDS("./mesocosmes_ITS/dadaFsITS.rds")
dadaRsITS <- readRDS("./mesocosmes_ITS/dadaRsITS.rds")
mergersITS <- readRDS("./mesocosmes_ITS/mergersITS.rds")
seqtabITS <- readRDS("./mesocosmes_ITS/seqtabITS.rds")
nochimITS <- readRDS("./mesocosmes_ITS/nochimITS.rds")
nochimseqtabITS <- readRDS("./mesocosmes_ITS/nochimseqtabITS.rds")
taxaITS <- readRDS("./mesocosmes_ITS/taxaITS.rds")

### We set a path variable to where the fastq files are located.
path <- "C:/Users/BaseSpace/mesocosmes/mesocosmes_ITS" # CHANGE THIS
# to the directory containing the unzipped fastq files in your computer.
list.files(path) # Verify the path by checking if the files are listed.

### This function allows us to look at the reads directly. This will
### confirm to us that the fastq files still have primers on, which
### will need to be removed.

### Here are the primers used for these amplifications (ITS1).
# Forward primer ITS1F:  CTTGGTCATTTAGAGGAAGTAA
# Reverse primer ITS2 :  GCTGCGTTCTTCATCGATGC

# Forward read: 
reads.R1 <- readFastq("./mesocosmes_ITS/SCB1-0-15-G-ITS_S101_L001_R1_001.fastq.gz")
sread(reads.R1)
# [1] CTTGGTCATTTAGAGGAAGTAAAAGTCGTAAC...TCTTGGTTCTGGCATCGATGAAGAACGCAGC
# Reverse read: 
reads.R2 <- readFastq("./mesocosmes_ITS/SCB1-0-15-G-ITS_S101_L001_R2_001.fastq.gz")
sread(reads.R2)
# [1] GCTGCGTTCTTCATCGATGCCAGAACCAAGAG...TTACGACTTTTACTTCCTCTAAATGACCAAG
### The primers appear in both the forward and reverse reads.
### Moreover, we can see here how the primers also appear in the distal end
### of each read. Notice how the last nucleotides of the forward read are
### complementary to the reverse primer; the complementary sequence of the 
### reverse primer is attached to the distal end of the forward read. 
### The same is true for the reverse read, which has the complementary 
### sequence of the forward primer at its end; 
### though not all reverse reads appear to be affected by this.

### Now, we perform some string manipulations to get matched 
### lists of the forward and reverse fastq files.
# Forward fastq filenames have format: SAMPLENAME_R1_001.fastq 
# Reverse fastq filenames have format: SAMPLENAME_R2_001.fastq
fnFsITS <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRsITS <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

### We now must remove the primers on all ends of all reads. The presence of
### the complementary sequence of primers on distal ends of reads complicates 
### the trimming process; this is where Cutadapt comes in handy.

### First, we identify the primers and their complementary sequence.

FWD <- "CTTGGTCATTTAGAGGAAGTAA"  
REV <- "GCTGCGTTCTTCATCGATGC"

allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
        RevComp = Biostrings::reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

### We must remove the sequences containing ambiguous (Ns) bases. The
### filtered files will be stored in a /filtN subdirectory.
fnFsITS.filtN <- file.path(path, "filtN", basename(fnFsITS))
fnRsITS.filtN <- file.path(path, "filtN", basename(fnRsITS))
# Note: as I have provided the filtered files already, you can
# skip the following step.
filterAndTrim(fnFsITS, fnFsITS.filtN, fnRsITS, fnRsITS.filtN, maxN = 0, multithread = FALSE)

### We now count the number of appearances of the forward and reverse primers
### in our forward and reverse reads.

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFsITS.filtN[[1]]),
	FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRsITS.filtN[[1]]), 
	REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFsITS.filtN[[1]]), 
	REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRsITS.filtN[[1]]))

#                  Forward Complement Reverse RevComp
# FWD.ForwardReads   21053          0       0       0
# FWD.ReverseReads       0          0       0    9541
# REV.ForwardReads       0          0       0   17037
# REV.ReverseReads   19205          0       0       0

### We are now ready to remove the primers. The user will need to install
### Cutadapt from: https://github.com/marcelm/cutadapt/releases/tag/v4.0
### I tried putting the program in the folder with the fastq files;
### perhaps that will prove enough, and you will not need to install the
### program from github. At any rate, once installed, we load it
### from the directory where we put it.

# Note: I have also provided the files after having been processed through
# Cutadapt; you may therefore skip the following steps. I will indicate
# where to pick up the script, a few lines below.

cutadapt <- "./mesocosmes_ITS/cutadapt.exe" # From the directory in which I put it.
system2(cutadapt, args = "--version") # Run shell commands from R 
# This command confirms that the program works and has been loaded correctly.
# The version of the program should appear as output of the command, e.g. 4.0.

# We now use Cutadapt to remove the primers. The files will be stored in
# a /cutadapt subdirectory.

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFsITS.cut <- file.path(path.cut, basename(fnFsITS))
fnRsITS.cut <- file.path(path.cut, basename(fnRsITS))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFsITS)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFsITS.cut[i], "-p", fnRsITS.cut[i], # output files
                             fnFsITS.filtN[i], fnRsITS.filtN[i])) # input files 
}

# We do a sanity check: number of primers after cutadapt.

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFsITS.cut[[1]]), 
	FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRsITS.cut[[1]]), 
	REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFsITS.cut[[1]]), 
	REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRsITS.cut[[1]]))

#                  Forward Complement Reverse RevComp
# FWD.ForwardReads       0          0       0       0
# FWD.ReverseReads       0          0       0       0
# REV.ForwardReads       0          0       0       0
# REV.ReverseReads       0          0       0       0

# As we can see, all primers have been removed, on all ends of all reads.
# We can now proceed with the filtering and subsequent steps of the pipeline,
# using the newly created cut files.

# If you skipped the Cutadapt steps, you can pick up the script from here.

path.cut <- file.path(path, "cutadapt")

# Forward and reverse fastq filenames have the format:
cutFsITS <- sort(list.files(path.cut, pattern="_R1_001.fastq", full.names = TRUE))
cutRsITS <- sort(list.files(path.cut, pattern="_R2_001.fastq", full.names = TRUE))
sample.names.its <- sapply(strsplit(basename(cutFsITS), "_"), `[`, 1)

### We check the quality profile of the forward and reverse reads.
### However, this will not really inform the trimming, as we will not
### assign a specific length to be trimmed. In other words, the argument
### TruncLen will not be used. However, to avoid spurious sequences that are
### very short, we will use the argument minLen, set to 50.

# Forward reads:
plotQualityProfile(fnFsITS[34:35]) # Or any other sample in []
# Reverse reads:
plotQualityProfile(fnRsITS[34:35]) # Or any other sample in []
# Reminder: if this function bugs, try restarting R.

filtFsITS <- file.path(path.cut, "filtered", basename(cutFsITS))
filtRsITS <- file.path(path.cut, "filtered", basename(cutRsITS))

outITS <- filterAndTrim(cutFsITS, filtFsITS, cutRsITS, filtRsITS, 
		maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50, 
		rm.phix = TRUE, compress = TRUE, 
		multithread = FALSE) # On Mac, set to number of cores.
# saveRDS(outITS, "./mesocosmes_ITS/outITS.rds")
head(outITS)

### The differences between the 16S and ITS pipelines have now all been
### addressed. The remainder of ITS pipeline is identical to that of 16S.

# Learn the Error Rates

errFITS <- learnErrors(filtFsITS, multithread = FALSE)
# saveRDS(errFITS, "./mesocosmes_ITS/errFITS.rds")
errRITS <- learnErrors(filtRsITS, multithread = FALSE)
# saveRDS(errRITS, "./mesocosmes_ITS/errRITS.rds")
plotErrors(errFITS, nominalQ = TRUE)
plotErrors(errRITS, nominalQ = TRUE)

### Sample inference

dadaFsITS <- dada(filtFsITS, err=errFITS, multithread=4)
# saveRDS(dadaFsITS, "./mesocosmes_ITS/dadaFsITS.rds")
dadaRsITS <- dada(filtRsITS, err=errRITS, multithread=4)
# saveRDS(dadaRsITS, "./mesocosmes_ITS/dadaRsITS.rds")
# dadaFsITS[[1]]
# 167 sequence variants were inferred from 3391 input unique sequences.

### Merging

mergersITS <- mergePairs(dadaFsITS, filtFsITS, dadaRsITS, filtRsITS,
	 verbose = TRUE)
# saveRDS(mergersITS, "./mesocosmes_ITS/mergersITS.rds")
# Inspect the merger data.frame from the first sample
# head(mergersITS[[1]])

### Constructing sequence table

seqtabITS <- makeSequenceTable(mergersITS)
# saveRDS(seqtabITS, "./mesocosmes_ITS/seqtabITS.rds")
dim(seqtabITS) # 81 1628 --> 1628 ASVs in 81 samples.
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtabITS)))
# Sequences vary widely in lengths, between 52 and 537.

### Removing chimeras

nochimITS <- removeBimeraDenovo(mergersITS, method = "consensus", 
	multithread = FALSE, verbose = TRUE)
# saveRDS(nochimITS, "./mesocosmes_ITS/nochimITS.rds")

nochimseqtabITS <- makeSequenceTable(nochimITS)
# saveRDS(nochimseqtabITS, "./mesocosmes_ITS/nochimseqtabITS.rds")
dim(nochimseqtabITS) # 81 1261 --> 1261 ASVs in 81 samples.
# Inspect proportion of reads kept after removal of chimeras
sum(nochimseqtabITS)/sum(seqtabITS)
# 96% of reads were preserved, so 4% of reads were removed as chimeras.

### Inspection

# Inspect number of reads preserved after each step of the overall processing
getN <- function(x) sum(getUniques(x))
trackITS <- cbind(outITS, sapply(dadaFsITS, getN), sapply(dadaRsITS, getN), 
	sapply(mergersITS, getN), rowSums(nochimseqtabITS))
# If processing a single sample, remove the sapply calls: 
# e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(trackITS) <- c("input", "filtered", "denoisedF", 
	"denoisedR", "merged", "nonchim")
rownames(trackITS) <- sample.names.its
head(trackITS)
#                            input filtered denoisedF denoisedR merged nonchim
# BC-0-15-G-ITS              27005    17459     17227     17222  16811   15819
# BC-15-30-G-ITS             21004    10914     10797     10746  10605   10352
# ctrl-PCR-neg-ChagnonP7-ITS  2266     2038      2028      2026   2026    2026
# GC-0-15-G-ITS              28559    17283     17020     16959  16507   16469
# GC-15-30-G-ITS             13284    10000      9779      9769   9537    9537
# P-BC-0-15-G-ITS            23452    14878     14691     14666  14209   13880

# As we can see, we managed to retain a large number of reads, which is
# reassuring. We did not lose almost all reads at the merging step, which
# indicates that the merging was succesful.

### We now assign taxonomy, by contrasting the sequences with a reference
### database. Here, we use UNITE. This database must normally be downloaded 
### by the user and placed in the directory with the fastq files. You will 
### find it provided in the folder, with the fastq and RDS files, so you
### do not need to download it.

# This function allows multithreading on Windows. However, it is so large
# that I could not run it on my computer; you will therefore most likely
# need to use a supercomputer and set the multithread argument to 20,
# if you wish to run the function on your end. Of course, this is not
# necessary, as I have provided the RDS output.

taxaITS <- assignTaxonomy(nochimseqtabITS, 
	"./mesocosmes_ITS/sh_general_release_dynamic_25.07.2023.fasta",
	multithread = 20)
# saveRDS(taxaITS, "./taxaITS.rds")

taxaITS.print <- taxaITS # Removing sequence rownames for display only
rownames(taxaITS.print) <- NULL
taxaITS.print[1:10,]
dim(taxaITS.print)
# We have 1261 ASVs, most of which seem to have been appropriately assigned.

##### This concludes the DADA2 pipeline for the ITS (fungal) data.

#############################################################################

#############################################################################
################################# PHYLOSEQ ##################################
#############################################################################

#############################################################################
################################### 16S #####################################
#############################################################################

##### Here begin the biostatistic analyses for the 16S (bacterial) data.

##### Tutorials :
# https://benjjneb.github.io/dada2/tutorial.html
# https://web.stanford.edu/class/bios221/Pune/Labs/Lab_phyloseq/Phyloseq_Lab.html#introduction_to_(introduction_to_phyloseq)
# https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/introduction.html
# https://astrobiomike.github.io/amplicon/dada2_workflow_ex#analysis-in-r

### First, we need to import our data, obtained through the DADA2 pipeline,
### into Phyloseq: a package that allows us to manipulate and analyse
### amplicon data.

# Note: the final output of the following steps, a phyloseq object
# summarising all the 16S data, has been saved as an RDS object.
# You need not go through the following steps; though they run quickly.

# To construct the phyloseq object, we need the following DADA2 outputs:
nochimseqtab16s <- readRDS("./mesocosmes_16S/nochimseqtab16s.rds")
taxa16s <- readRDS("./mesocosmes_16S/taxa16s.rds")

# We also need a samples dataframe, containing additional information
# about our samples, such as environmental variables or, in our case,
# information about the treatments, substrates, plant species, and depths.
samples <- read.csv("./sample_data_mesocosmes.csv")

# We check the nature of the variables in the samples dataframe.
summary(samples)
# The mesocosm, species, depth, substrate, pesticides, concentration, biochar,
# and the compound concentration_biochar variables must be converted.
samples$Mesocosm <- as.factor(samples$Mesocosm)
samples$Species <- as.factor(samples$Species)
samples$Depth <- as.factor(samples$Depth)
samples$Substrate <- as.factor(samples$Substrate)
samples$Pesticides <- as.factor(samples$Pesticides)
samples$Concentration <- as.numeric(samples$Concentration)
samples$Pesticide_Concentration <- as.factor(samples$Pesticide_Concentration)
samples$Biochar <- as.factor(samples$Biochar)
samples$Concentration_Biochar <- as.factor(samples$Concentration_Biochar)
# We confirm that it worked.
summary(samples)

# We manipulate the various dataframes to ensure they have identical rows.
desired_order <- paste0(samples[,1],"-16S")
rownames(samples) = desired_order
samples <- samples[,-1]
nochimseqtab16s_reordered <- nochimseqtab16s[-c(3,4),][desired_order,]
# [-c(3,4),] removes the negative and positive controls.

# We create the phyloseq object.
ps.meso.16s <- phyloseq(
			otu_table(nochimseqtab16s_reordered, taxa_are_rows = FALSE), 
			sample_data(samples), 
			tax_table(taxa16s))

# This is what a phyloseq object looks like:
print(ps.meso.16s)
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7299 taxa and 80 samples ]
# sample_data() Sample Data:       [ 80 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 7299 taxa by 6 taxonomic ranks ]

### The beauty of a phyloseq object is that it links together various
### dataframes. In our case, we have three linked dataframes.
# otu_table:   site x species dataframe, with sites as rows 
# 		   and sequences as columns;
# sample_data: a dataframe of environmental variables for each site;
# tax_table:   the taxonomic classification of all the sequences found
# 		   in the dataset. 
# More dataframes could be linked, such as a phylogenetic tree, for instance.

# Rename taxa to short strings
dna.meso.16s <- Biostrings::DNAStringSet(taxa_names(ps.meso.16s))
names(dna.meso.16s) <- taxa_names(ps.meso.16s)

# Add the above dataframe to the phyloseq object.
ps.meso.16s <- merge_phyloseq(ps.meso.16s, dna.meso.16s)
taxa_names(ps.meso.16s) <- paste0("ASV", seq(ntaxa(ps.meso.16s)))

# Remove uncharacterised phyla.
ps.meso.16s = subset_taxa(ps.meso.16s, !is.na(Phylum) & 
	!Phylum %in% c("", "uncharacterized"))

# Fix taxonomy to remove further NAs and equivalents.
ps.meso.16s <- ps.meso.16s %>%
 tax_fix(
  min_length = 4,
  unknowns = c(""),
  sep = " ", anon_unique = TRUE,
  suffix_rank = "classified"
 )

# Remove all eukaryota, archaea, chloroplast, and mitochondria taxa.
ps.meso.16s <- subset_taxa(ps.meso.16s, Kingdom != "Eukaryota")
ps.meso.16s <- subset_taxa(ps.meso.16s, Kingdom != "Archaea")
ps.meso.16s <- subset_taxa(ps.meso.16s, Order != "Chloroplast")
ps.meso.16s <- subset_taxa(ps.meso.16s, Family != "Mitochondria")

# Subset the phyloseq object to retain only the data for the biochar project.
ps.meso.16s <- subset_samples(ps.meso.16s, Species == "Scirp.cyp")

# Remove the phyla that have zero counts throughout all samples.
ps.meso.16s <- prune_taxa(taxa_sums(ps.meso.16s) > 0, ps.meso.16s)

### Our phyloseq object is now ready for downstream analyses.
print(ps.meso.16s)
# otu_table()   OTU Table:         [ 4530 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 4530 taxa by 6 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 4530 reference sequences ]

### For future analyses and graphs, it will be very helpful to have
### the concentration levels ordered None-->High, rather than alphabetically.

sample_data(ps.meso.16s)$Pesticide_Concentration <- 
	factor(sample_data(ps.meso.16s)$Pesticide_Concentration, 
	levels = c("None", "Low", "Medium", "High"))

# saveRDS(ps.meso.16s, "./mesocosmes_16S/ps_meso_16s.rds")

### If you did not go through these steps, you can access the phyloseq
### object here:

##########             PROCEED WITH THE FOLLOWING LINES             ##########
##########            	    TO RAREFY THE DATASETS.                 ##########

ps.meso.16s <- readRDS("./mesocosmes_16S/ps_meso_16s.rds")

### We now check the rarefaction curves, to see if we need to rarefy the 
### samples due to large differences in sequencing depth between samples.

### Check number of reads
summary(sample_sums(ps.meso.16s))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     666    4284    6848    7393   10540   16288  

# As we can see, there is a large difference between the
# smallest and largest sequencing depths.

### We will now plot the sequencing depths.
# png("./Figures/16s_reads_per_sample.png",
#	width = 1000, height = 500)
barplot(sample_sums(ps.meso.16s), 
	las = 2, cex.names = 0.5, col = "red",
	main = "Number of reads")
abline(h = 0)
# dev.off()

# There is a large difference in the number of reads.
# We will plot the rarefaction curves.

# png("./Figures/16s_rarefaction_curves.png",
#	width = 1000, height = 500)
rarecurve(t(abundances(ps.meso.16s)),
	step = 100, ylab = "ASVs", label = FALSE)
abline(v = min(sample_sums(ps.meso.16s)), col = "red") # smallest depth
# dev.off()

rar.ps.meso.16s <- rarefy_even_depth(ps.meso.16s,
			sample.size = min(sample_sums(ps.meso.16s)),
			rngseed = 33, replace = FALSE, trimOTUs = TRUE)

# We now check the number of reads per sample and the rarefaction curves.

barplot(sample_sums(rar.ps.meso.16s), 
	las = 2, cex.names = 0.5, col = "red",
	main = "Number of reads")
abline(h = 0)

rarecurve(t(abundances(rar.ps.meso.16s)),
	step = 100, ylab = "ASVs", label = FALSE)
abline(v = min(sample_sums(rar.ps.meso.16s)), 
	col = "red") # smallest depth

# The rarefaction worked.
print(rar.ps.meso.16s)
# otu_table()   OTU Table:         [ 2686 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 2686 taxa by 6 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 2686 reference sequences ]

# We can save it as a RDS object, so that we do not need to run
# the whole rarefaction each time.

# saveRDS(rar.ps.meso.16s, "./mesocosmes_16S/rar_ps_meso_16s.rds")

#############################################################################
#################################### ITS ####################################
#############################################################################

##### Here begin the biostatistic analyses for the ITS (fungal) data.

##### Tutorials :
# https://benjjneb.github.io/dada2/tutorial.html
# https://web.stanford.edu/class/bios221/Pune/Labs/Lab_phyloseq/Phyloseq_Lab.html#introduction_to_(introduction_to_phyloseq)
# https://mibwurrepo.github.io/Microbial-bioinformatics-introductory-course-Material-2018/introduction.html
# https://astrobiomike.github.io/amplicon/dada2_workflow_ex#analysis-in-r

### First, we need to import our data, obtained through the DADA2 pipeline,
### into Phyloseq: a package that allows us to manipulate and analyse
### amplicon data.

# Note: the final output of the following steps, a phyloseq object
# summarising all the ITS data, has been saved as an RDS object.
# You need not go through the following steps; though they run quickly.

# To construct the phyloseq object, we need the following DADA2 outputs:
nochimseqtabITS <- readRDS("./mesocosmes_ITS/nochimseqtabITS.rds")
taxaITS <- readRDS("./mesocosmes_ITS/taxaITS.rds")

# We also need a samples dataframe, containing additional information
# about our samples, such as environmental variables or, in our case,
# information about the treatments, substrates, plant species, and depths.
samples <- read.csv("./sample_data_mesocosmes.csv")

# We check the nature of the variables in the samples dataframe.
summary(samples)
# The mesocosm, species, depth, substrate, pesticides, concentration, biochar,
# and the compound concentration_biochar variables must be converted.
samples$Mesocosm <- as.factor(samples$Mesocosm)
samples$Species <- as.factor(samples$Species)
samples$Depth <- as.factor(samples$Depth)
samples$Substrate <- as.factor(samples$Substrate)
samples$Pesticides <- as.factor(samples$Pesticides)
samples$Concentration <- as.numeric(samples$Concentration)
samples$Pesticide_Concentration <- as.factor(samples$Pesticide_Concentration)
samples$Biochar <- as.factor(samples$Biochar)
samples$Concentration_Biochar <- as.factor(samples$Concentration_Biochar)
# We confirm that it worked.
summary(samples)

# We manipulate the various dataframes to ensure they have identical rows.
desired_order <- paste0(samples[,1],"-ITS")
rownames(samples) = desired_order
samples <- samples[,-1]
rownames(nochimseqtabITS) <- sub("_S.*", "", rownames(nochimseqtabITS))
nochimseqtabITS_reordered <- nochimseqtabITS[-3,][desired_order,]
# [-3,] removes the negative control.

# We create the phyloseq object.
ps.meso.its <- phyloseq(
			otu_table(nochimseqtabITS_reordered, taxa_are_rows = FALSE), 
			sample_data(samples), 
			tax_table(taxaITS))

# Rename taxa to short strings
dna.meso.its <- Biostrings::DNAStringSet(taxa_names(ps.meso.its))
names(dna.meso.its) <- taxa_names(ps.meso.its)

# Add the above dataframe to the phyloseq object.
ps.meso.its <- merge_phyloseq(ps.meso.its, dna.meso.its)
taxa_names(ps.meso.its) <- paste0("ASV", seq(ntaxa(ps.meso.its)))

# Remove uncharacterised phyla
ps.meso.its = subset_taxa(ps.meso.its, !is.na(Phylum) & 
	!Phylum %in% c("", "uncharacterized"))

# Fix taxonomy to remove further NAs and equivalents.
ps.meso.its <- ps.meso.its %>%
 tax_fix(
  min_length = 4,
  unknowns = c(""),
  sep = " ", anon_unique = TRUE,
  suffix_rank = "classified"
 )

### We have one additional cleaning step to do. If you run the following:
tax_table(ps.meso.its)[1,]
### you will notice that all taxonomic rank have a prefix.
#      Kingdom    Phylum          Class                Order            
# ASV1 "k__Fungi" "p__Ascomycota" "c__Dothideomycetes" "o__Pleosporales"
#      Family                               Genus             Species       
# ASV1 "f__Pleosporales_fam_Incertae_sedis" "g__Wettsteinina" "s__lacustris"

### We will now remove these prefixes.
tax_table(ps.meso.its)[, colnames(tax_table(ps.meso.its))] <- 
	gsub(tax_table(ps.meso.its)[, colnames(tax_table(ps.meso.its))], 
		pattern = "[a-z]__", replacement = "")

# Subset the phyloseq object to retain only the data for the biochar project.
ps.meso.its <- subset_samples(ps.meso.its, Species == "Scirp.cyp")

# Remove the phyla that have zero counts throughout all samples.
ps.meso.its <- prune_taxa(taxa_sums(ps.meso.its) > 0, ps.meso.its)

### Our phyloseq object is now ready for downstream analyses.
print(ps.meso.its)
# otu_table()   OTU Table:         [ 511 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 511 taxa by 7 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 511 reference sequences ]

### For future analyses and graphs, it will be very helpful to have
### the concentration levels ordered None-->High, rather than alphabetically.

sample_data(ps.meso.its)$Pesticide_Concentration <- 
	factor(sample_data(ps.meso.its)$Pesticide_Concentration, 
	levels = c("None", "Low", "Medium", "High"))

# saveRDS(ps.meso.its, "./mesocosmes_ITS/ps_meso_its.rds")

### If you did not go through these steps, you can access the phyloseq
### object here:

##########               RUN THE FOLLOWING LINE TO LOAD             ##########
##########                 THE MAIN PHYLOSEQ OBJECT.                ##########

ps.meso.its <- readRDS("./mesocosmes_ITS/ps_meso_its.rds")

##########             PROCEED WITH THE FOLLOWING LINES             ##########
##########            IF YOU WISH TO RAREFY THE DATASETS.           ##########

### We now check the rarefaction curves, to see if we need to rarefy the 
### samples due to large differences in sequencing depth between samples.

### Check number of reads
summary(sample_sums(ps.meso.its))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    2203    7448    9394   11379   14255   28401 

# As we can see, there is a large difference between the
# smallest and largest sequencing depths.

### We will now plot the sequencing depths.
# png("./Figures/ITS_reads_per_sample.png",
#	width = 1000, height = 500)
barplot(sample_sums(ps.meso.its), 
	las = 2, cex.names = 0.5, col = "red",
	main = "Number of reads")
abline(h = 0)
# dev.off()

# There is a large difference in the number of reads.
# We will plot the rarefaction curves.

# png("./Figures/ITS_rarefaction_curves.png",
#	width = 1000, height = 500)
rarecurve(t(abundances(ps.meso.its)),
	step = 100, ylab = "ASVs", label = FALSE)
abline(v = min(sample_sums(ps.meso.its)), col = "red") # smallest depth
# dev.off()

rar.ps.meso.its <- rarefy_even_depth(ps.meso.its,
			sample.size = min(sample_sums(ps.meso.its)),
			rngseed = 33, replace = FALSE, trimOTUs = TRUE)

# We now check the number of reads per sample and the rarefaction curves.

barplot(sample_sums(rar.ps.meso.its), 
	las = 2, cex.names = 0.5, col = "red",
	main = "Number of reads")
abline(h = 0)

rarecurve(t(abundances(rar.ps.meso.its)),
	step = 100, ylab = "ASVs", label = FALSE)
abline(v = min(sample_sums(rar.ps.meso.its)), 
	col = "red") # smallest depth

# The rarefaction worked.
print(rar.ps.meso.its)
# otu_table()   OTU Table:         [ 445 taxa and 48 samples ]
# sample_data() Sample Data:       [ 48 samples by 9 sample variables ]
# tax_table()   Taxonomy Table:    [ 445 taxa by 7 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 445 reference sequences ]

# We can save it as a RDS object, so that we do not need to run
# the whole rarefaction each time.

# saveRDS(rar.ps.meso.its, "./mesocosmes_ITS/rar_ps_meso_its.rds")

##############################################################################
####################             STATISTICAL              ####################
####################              ANALYSES                ####################
##############################################################################

# Load the unrarified and the rarefied 16S and ITS phyloseq objects.

ps.meso.16s <- readRDS("./mesocosmes_16S/ps_meso_16s.rds")
rar.ps.meso.16s <- readRDS("./mesocosmes_16S/rar_ps_meso_16s.rds")

ps.meso.its <- readRDS("./mesocosmes_ITS/ps_meso_its.rds")
rar.ps.meso.its <- readRDS("./mesocosmes_ITS/rar_ps_meso_its.rds")

##########    RUN THESE LINES TO GENERATE THE DIVERSITY INDICES     ##########
##########                  FOR THE RAREFIED DATA                   ##########

### Create a dataframe with many diversity indices for the rarefied data.

rar.div.meso.16s <- alpha(rar.ps.meso.16s, index = "all") # a lot of indices
# We add nmber of reads/samples.
rar.div.meso.16s$ReadsPerSample <- sample_sums(rar.ps.meso.16s)
# We add Hill numbers for Shannon and Simpson
rar.hill.16s <- renyi(otu_table(rar.ps.meso.16s), 
	hill = TRUE, 
	scales = c(1,2))
colnames(rar.hill.16s) <- c("hill_shannon", "hill_simpson")
identical(rownames(rar.hill.16s),rownames(rar.div.meso.16s))
rar.div.meso.16s <- cbind(rar.div.meso.16s, rar.hill.16s)
# get the metadata out as separate object
rar.meta.16s <- meta(rar.ps.meso.16s)
identical(rownames(rar.meta.16s),rownames(rar.div.meso.16s))
rar.div.meso.16s <- cbind(rar.div.meso.16s, rar.meta.16s)
colnames(rar.div.meso.16s)
summary(rar.div.meso.16s)

rar.div.meso.its <- alpha(rar.ps.meso.its, index = "all") # a lot of indices
# We add nmber of reads/samples.
rar.div.meso.its$ReadsPerSample <- sample_sums(rar.ps.meso.its)
# We add Hill numbers for Shannon and Simpson
rar.hill.its <- renyi(otu_table(rar.ps.meso.its), 
	hill = TRUE, 
	scales = c(1,2))
colnames(rar.hill.its) <- c("hill_shannon", "hill_simpson")
identical(rownames(rar.hill.its),rownames(rar.div.meso.its))
rar.div.meso.its <- cbind(rar.div.meso.its, rar.hill.its)
# get the metadata out as separate object
rar.meta.its <- meta(rar.ps.meso.its)
identical(rownames(rar.meta.its),rownames(rar.div.meso.its))
rar.div.meso.its <- cbind(rar.div.meso.its, rar.meta.its)
colnames(rar.div.meso.its)
summary(rar.div.meso.its)

########## Alpha-Diversity ########## 

### We will have separate analyses before and after pesticide treatment.
### before pesticides, we will evaluate the impact of biochar, with a
### linear mixed model. After pesticides, we will evaluate the impact
### of biochar and pesticide, in linear models separately for each depth.
### Data were log-transformed whenever the residuals of the linear model 
### were not normally distributed, or the data were heteroscedastic.

### Impact of biochar before treatment

### 16S

## Gravel

lm.meso.16s.g.b <- lmer(hill_shannon ~ 
	Biochar + Depth + (1 | Mesocosm), 
	data = subset(rar.div.meso.16s, 
		Substrate == "Gravel" & Pesticides == "Pesticides-"))
hist(resid(lm.meso.16s.g.b), col = "red", breaks = 10)
shapiro.test(resid(lm.meso.16s.g.b))
qqPlot(resid(lm.meso.16s.g.b))
anova(lm.meso.16s.g.b, type = "3")
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# Biochar 339.17  339.17     1     4  0.1591 0.7104
# Depth   501.35  501.35     1     5  0.2352 0.6482

## Roots

lm.meso.16s.r.b <- lmer(hill_shannon ~ 
	Biochar + Depth + (1 | Mesocosm), 
	data = subset(rar.div.meso.16s, 
		Substrate == "Roots" & Pesticides == "Pesticides-"))
hist(resid(lm.meso.16s.r.b), col = "red", breaks = 10)
shapiro.test(resid(lm.meso.16s.r.b))
qqPlot(resid(lm.meso.16s.r.b))
anova(lm.meso.16s.r.b, type = "3")
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# Biochar 4799.6  4799.6     1     9  3.5464 0.09234 .
# Depth    999.5   999.5     1     9  0.7385 0.41243 

### ITS

## Gravel

lm.meso.its.g.b <- lmer(hill_shannon ~ 
	Biochar + Depth + (1 | Mesocosm), 
	data = subset(rar.div.meso.its, 
		Substrate == "Gravel" & Pesticides == "Pesticides-"))
hist(resid(lm.meso.its.g.b), col = "red", breaks = 10)
shapiro.test(resid(lm.meso.its.g.b))
qqPlot(resid(lm.meso.its.g.b))
anova(lm.meso.its.g.b, type = "3")
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# Biochar 119.487 119.487     1     9  2.7145 0.1338
# Depth    87.599  87.599     1     9  1.9900 0.1920

## Roots

lm.meso.its.r.b <- lmer(hill_shannon ~ 
	Biochar + Depth + (1 | Mesocosm), 
	data = subset(rar.div.meso.its, 
		Substrate == "Roots" & Pesticides == "Pesticides-"))
hist(resid(lm.meso.its.r.b), col = "red", breaks = 10)
shapiro.test(resid(lm.meso.its.r.b))
qqPlot(resid(lm.meso.its.r.b))
anova(lm.meso.its.r.b, type = "3")
# Type III Analysis of Variance Table with Satterthwaite's method
#         Sum Sq Mean Sq NumDF DenDF F value Pr(>F)
# Biochar 17.007  17.007     1     9  0.2976 0.5986
# Depth   72.346  72.346     1     9  1.2661 0.2896 

### Impact of biochar and pesticide concentration after treatment

### 16S

## Gravel - 0-15 cm

lm.meso.16s.g.15 <- aov(hill_shannon ~
	Biochar * Concentration,
	data = subset(rar.div.meso.16s, 
		Substrate == "Gravel" &
		Depth == "0-15 cm" &
		Pesticides == "Pesticides+")
	)

hist(resid(lm.meso.16s.g.15), col = "red", breaks = 10)
shapiro.test(resid(lm.meso.16s.g.15))
qqPlot(resid(lm.meso.16s.g.15))
anova(lm.meso.16s.g.15)
#                       Df Sum Sq Mean Sq F value Pr(>F)
# Biochar                1   4521    4521   1.303  0.372
# Concentration          1    994     994   0.286  0.646
# Biochar:Concentration  1   6285    6285   1.812  0.311
# Residuals              2   6938    3469 

## Gravel - 15-30 cm

lm.meso.16s.g.30 <- aov(hill_shannon ~
	Biochar * Concentration,
	data = subset(rar.div.meso.16s, 
		Substrate == "Gravel" &
		Depth == "15-30 cm" &
		Pesticides == "Pesticides+")
	)

hist(resid(lm.meso.16s.g.30), col = "red", breaks = 10)
shapiro.test(resid(lm.meso.16s.g.30))
qqPlot(resid(lm.meso.16s.g.30))
anova(lm.meso.16s.g.30)
#                       Df Sum Sq Mean Sq F value Pr(>F)
# Biochar                1   3923    3923   0.454  0.570
# Concentration          1   4357    4357   0.504  0.551
# Biochar:Concentration  1   3412    3412   0.395  0.594
# Residuals              2  17284    8642  

## Roots - 0-15 cm

lm.meso.16s.r.15 <- aov(hill_shannon ~
	Biochar * Concentration,
	data = subset(rar.div.meso.16s, 
		Substrate == "Roots" &
		Depth == "0-15 cm" &
		Pesticides == "Pesticides+")
	)

hist(resid(lm.meso.16s.r.15), col = "red", breaks = 10)
shapiro.test(resid(lm.meso.16s.r.15))
qqPlot(resid(lm.meso.16s.r.15))
anova(lm.meso.16s.r.15)
#                       Df Sum Sq Mean Sq F value Pr(>F)
# Biochar                1  560.3   560.3   0.390  0.596
# Concentration          1  748.2   748.2   0.521  0.545
# Biochar:Concentration  1  152.0   152.0   0.106  0.776
# Residuals              2 2872.0  1436.0  

## Roots - 15-30 cm

lm.meso.16s.r.30 <- aov(hill_shannon ~
	Biochar * Concentration,
	data = subset(rar.div.meso.16s, 
		Substrate == "Roots" &
		Depth == "15-30 cm" &
		Pesticides == "Pesticides+")
	)

hist(resid(lm.meso.16s.r.30), col = "red", breaks = 10)
shapiro.test(resid(lm.meso.16s.r.30))
qqPlot(resid(lm.meso.16s.r.30))
anova(lm.meso.16s.r.30)
#                       Df Sum Sq Mean Sq F value Pr(>F)
# Biochar                1     82    81.6   0.043  0.855
# Concentration          1   1866  1865.5   0.987  0.425
# Biochar:Concentration  1   2235  2235.2   1.183  0.390
# Residuals              2   3778  1889.2   

### ITS

## Gravel - 0-15 cm

lm.meso.its.g.15 <- aov(hill_shannon ~
	Biochar * Concentration,
	data = subset(rar.div.meso.its, 
		Substrate == "Gravel" &
		Depth == "0-15 cm" &
		Pesticides == "Pesticides+")
	)

hist(resid(lm.meso.its.g.15), col = "red", breaks = 10)
shapiro.test(resid(lm.meso.its.g.15))
qqPlot(resid(lm.meso.its.g.15))
anova(lm.meso.its.g.15)
#                       Df Sum Sq Mean Sq F value Pr(>F)
# Biochar                1  93.43   93.43   1.224  0.384
# Concentration          1  30.89   30.89   0.405  0.590
# Biochar:Concentration  1  40.61   40.61   0.532  0.542
# Residuals              2 152.71   76.35

## Gravel - 15-30 cm

lm.meso.its.g.30 <- aov(hill_shannon ~
	Biochar * Concentration,
	data = subset(rar.div.meso.its, 
		Substrate == "Gravel" &
		Depth == "15-30 cm" &
		Pesticides == "Pesticides+")
	)

hist(resid(lm.meso.its.g.30), col = "red", breaks = 10)
shapiro.test(resid(lm.meso.its.g.30))
qqPlot(resid(lm.meso.its.g.30))
anova(lm.meso.its.g.30)
#                       Df Sum Sq Mean Sq F value Pr(>F)
# Biochar                1 124.51  124.51   3.199  0.216
# Concentration          1  11.66   11.66   0.300  0.639
# Biochar:Concentration  1  57.82   57.82   1.486  0.347
# Residuals              2  77.84   38.92 

## Roots - 0-15 cm

lm.meso.its.r.15 <- aov(hill_shannon ~
	Biochar * Concentration,
	data = subset(rar.div.meso.its, 
		Substrate == "Roots" &
		Depth == "0-15 cm" &
		Pesticides == "Pesticides+")
	)

hist(resid(lm.meso.its.r.15), col = "red", breaks = 10)
shapiro.test(resid(lm.meso.its.r.15))
qqPlot(resid(lm.meso.its.r.15))
anova(lm.meso.its.r.15)
#                       Df Sum Sq Mean Sq F value Pr(>F)
# Biochar                1  2.881   2.881   0.771  0.473
# Concentration          1  0.219   0.219   0.059  0.831
# Biochar:Concentration  1  0.105   0.105   0.028  0.882
# Residuals              2  7.476   3.738 

## Roots - 15-30 cm

lm.meso.its.r.30 <- aov(hill_shannon ~
	Biochar * Concentration,
	data = subset(rar.div.meso.its, 
		Substrate == "Roots" &
		Depth == "15-30 cm" &
		Pesticides == "Pesticides+")
	)

hist(resid(lm.meso.its.r.30), col = "red", breaks = 10)
shapiro.test(resid(lm.meso.its.r.30))
qqPlot(resid(lm.meso.its.r.30))
anova(lm.meso.its.r.30)
#                       Df Sum Sq Mean Sq F value Pr(>F)
# Biochar                1   9.38    9.38   0.855  0.453
# Concentration          1  31.70   31.70   2.889  0.231
# Biochar:Concentration  1   1.34    1.34   0.122  0.760
# Residuals              2  21.94   10.97  

### Figure

plot.div.16s <-
   ggboxplot(
	subset(rar.div.meso.16s),
		x = "Pesticide_Concentration",
		y = "hill_shannon",
		alpha = "Biochar",
		fill = "Pesticide_Concentration",
		title = "Bacteria",
		legend = "right") +
	labs(y = "Hill-Shannon Index", x = "Concentration") +
	rotate_x_text() + 
	scale_alpha_manual(values = c(0.35, 1)) +
	scale_fill_manual(
		values = c("#C688D0","#2FDC2A","#394EF5","#F90808")) +
	theme(axis.text = element_text(size = 10)) +
	facet_wrap(vars(Substrate), scales = "free_x", ncol = 2) +
	guides(alpha = "none") +
	geom_vline(xintercept = c(0,1.5), linetype = "dashed") +
	theme(
		legend.position = "none",
		panel.background = element_rect(linewidth = 0.5, colour = "black"))

plot.div.its <-
   ggboxplot(
	subset(rar.div.meso.its),
		x = "Pesticide_Concentration",
		y = "hill_shannon",
		alpha = "Biochar",
		fill = "Pesticide_Concentration",
		title = "Fungi",
		legend = "right") +
	labs(y = "Hill-Shannon Index", x = "Concentration") +
	rotate_x_text() + 
	scale_alpha_manual(values = c(0.35, 1)) +
	scale_fill_manual(
		values = c("#C688D0","#2FDC2A","#394EF5","#F90808")) +
	theme(axis.text = element_text(size = 10)) +
	facet_wrap(vars(Substrate), scales = "free_x", ncol = 2) +
	guides(alpha = "none") +
	geom_vline(xintercept = c(0,1.5), linetype = "dashed") +
	theme(
		legend.position = "none",
		panel.background = element_rect(linewidth = 0.5, colour = "black"))

ggarrange(plot.div.16s, plot.div.its, ncol = 2, nrow = 1)
# ggsave("./Figures/alpha_div.png",
	width = 8, height = 4, unit = "in")
# ggsave("./Figures/Publication/FIG1.pdf",
	width = 190, height = 95, unit = "mm")

########## Taxonomic Composition ########## 

### We will compute a permanova with the Hellinger distance.
### The permutation blocks must be built in a way to take into account the
### random factor of mesocosms. I will use Permute, built-in within Vegan.
### The mesocosm corresponds to the pesticides concentration and to the 
### presence of biochar. From each mesocosm, we have samples taken at two
### time points: before pesticides, and after pesticides. I cannot allow
### permutations between these samples, because I do not want to test for
### the impact of time. In fact, due to the time distance between the two
### sampling periods (before and after pesticides), drastic differences in
### community composition are to be expected, but I am not interested in these
### differences. Therefore, I will set the pesticides factor as Blocks, since
### Blocks are never permuted. I will set the mesocosms at the Plot level; 
### this will forbid permutations of depth levels between different mesocosms,
### ensuring that the mesocosms are permuted as a whole. Individual mesocosms 
### (proxy for pesticide concentration and presence of biochar) will be 
### permuted within each sampling period (before and after pesticides). Also, 
### within each mesocosms, I will allow the permutation of depths.

### 16S

## Gravel

perm.meso.16s.g <- with(
  meta(subset_samples(rar.ps.meso.16s, Substrate == "Gravel")), 
  how(
    nperm = 999,
    blocks = Pesticides,
    plots = Plots(strata = Mesocosm, type = "free"),
    within = Within(type = "free")
  	)
   )
set.seed(100)
adonis2(otu_table(subset_samples(rar.ps.meso.16s, Substrate == "Gravel")) ~ 
		Biochar * Concentration + Depth,
	data = meta(subset_samples(rar.ps.meso.16s, Substrate == "Gravel")),
	permutation = perm.meso.16s.g,
	method = "hellinger",
	by = "terms")
#                       Df SumOfSqs      R2      F Pr(>F)    
# Biochar                1   1.0042 0.06273 1.5593  0.003 ** 
# Concentration          1   1.1846 0.07399 1.8393  0.707    
# Depth                  1   0.9061 0.05660 1.4069  0.001 ***
# Biochar:Concentration  1   0.6780 0.04235 1.0528  0.690    
# Residual              19  12.2363 0.76433                  
# Total                 23  16.0092 1.00000     

## Roots

perm.meso.16s.r <- with(
  meta(subset_samples(rar.ps.meso.16s, Substrate == "Roots")), 
  how(
    nperm = 999,
    blocks = Pesticides,
    plots = Plots(strata = Mesocosm, type = "free"),
    within = Within(type = "free")
  	)
   )
set.seed(100)
adonis2(otu_table(subset_samples(rar.ps.meso.16s, Substrate == "Roots")) ~ 
		Biochar * Concentration + Depth,
	data = meta(subset_samples(rar.ps.meso.16s, Substrate == "Roots")),
	permutation = perm.meso.16s.r,
	method = "hellinger",
	by = "terms")
#                       Df SumOfSqs      R2      F Pr(>F)   
# Biochar                1   0.9760 0.05840 1.5119  0.004 **
# Concentration          1   2.1092 0.12620 3.2673  0.496   
# Depth                  1   0.6407 0.03834 0.9926  0.044 * 
# Biochar:Concentration  1   0.7214 0.04317 1.1175  0.401   
# Residual              19  12.2650 0.73389                 
# Total                 23  16.7123 1.00000  

### ITS

## Gravel

perm.meso.its.g <- with(
  meta(subset_samples(rar.ps.meso.its, Substrate == "Gravel")), 
  how(
    nperm = 999,
    blocks = Pesticides,
    plots = Plots(strata = Mesocosm, type = "free"),
    within = Within(type = "free")
  	)
   )
set.seed(100)
adonis2(otu_table(subset_samples(rar.ps.meso.its, Substrate == "Gravel")) ~ 
		Biochar * Concentration + Depth,
	data = meta(subset_samples(rar.ps.meso.its, Substrate == "Gravel")),
	permutation = perm.meso.its.g,
	method = "hellinger",
	by = "terms")
#                       Df SumOfSqs      R2      F Pr(>F)    
# Biochar                1   0.5026 0.05022 1.3211  0.398    
# Concentration          1   0.7731 0.07725 2.0323  0.525    
# Depth                  1   0.9785 0.09778 2.5723  0.001 ***
# Biochar:Concentration  1   0.5257 0.05253 1.3818  0.266    
# Residual              19   7.2280 0.72223                  
# Total                 23  10.0079 1.00000     

## Roots

perm.meso.its.r <- with(
  meta(subset_samples(rar.ps.meso.its, Substrate == "Roots")), 
  how(
    nperm = 999,
    blocks = Pesticides,
    plots = Plots(strata = Mesocosm, type = "free"),
    within = Within(type = "free")
  	)
   )
set.seed(100)
adonis2(otu_table(subset_samples(rar.ps.meso.its, Substrate == "Roots")) ~ 
		Biochar * Concentration + Depth,
	data = meta(subset_samples(rar.ps.meso.its, Substrate == "Roots")),
	permutation = perm.meso.its.r,
	method = "hellinger",
	by = "terms")
#                       Df SumOfSqs      R2      F Pr(>F)   
# Biochar                1   0.4205 0.04708 1.2371  0.050 *
# Concentration          1   0.7496 0.08392 2.2054  0.015 *
# Depth                  1   0.9862 0.11041 2.9015  0.012 *
# Biochar:Concentration  1   0.3177 0.03557 0.9347  0.201  
# Residual              19   6.4579 0.72302                
# Total                 23   8.9319 1.00000  

### Figures

## Biochar without pesticides

plot.pca.16s.p0.g <-
  subset_samples(rar.ps.meso.16s,
	Pesticides == "Pesticides-" & Substrate == "Gravel") %>%
  tax_transform("hellinger") %>%
  ord_calc(method = "PCA") %>%
    ord_plot(
	shape = "Biochar",
	colour = "#C688D0",
	size = 2,
    	tax_vec_length = 0.5,
	auto_caption = NA,
    	tax_lab_style = tax_lab_style(
		type = "text",
		max_angle = 90, 
		aspect_ratio = 1, 
		justify = "auto",
		size = 4)) +
	scale_shape_manual(values = c(1,16)) +
	scale_colour_manual(
		name = "Concentration",
		values = c("#C688D0")) +
	theme_bw() + ggtitle("Bacteria - Gravel") +
	stat_ellipse(aes(linetype = Biochar), level = 0.95) +
	scale_linetype_manual(values = c(2,1)) +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		legend.text = element_text(size = 11))

plot.pca.16s.p0.r <-
  subset_samples(rar.ps.meso.16s,
	Pesticides == "Pesticides-" & Substrate == "Roots") %>%
  tax_transform("hellinger") %>%
  ord_calc(method = "PCA") %>%
    ord_plot(
	shape = "Biochar",
	colour = "#C688D0",
	size = 2,
    	tax_vec_length = 0.5,
	auto_caption = NA,
    	tax_lab_style = tax_lab_style(
		type = "text",
		max_angle = 90, 
		aspect_ratio = 1, 
		justify = "auto",
		size = 4)) +
	scale_shape_manual(values = c(1,16)) +
	scale_colour_manual(
		name = "Concentration",
		values = c("#C688D0")) +
	theme_bw() + ggtitle("Bacteria - Roots") +
	stat_ellipse(aes(linetype = Biochar), level = 0.95) +
	scale_linetype_manual(values = c(2,1)) +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		legend.text = element_text(size = 11))

plot.pca.its.p0.g <-
  subset_samples(rar.ps.meso.its,
	Pesticides == "Pesticides-" & Substrate == "Gravel") %>%
  tax_transform("hellinger") %>%
  ord_calc(method = "PCA") %>%
    ord_plot(
	shape = "Biochar",
	colour = "#C688D0",
	size = 2,
    	tax_vec_length = 0.5,
	auto_caption = NA,
    	tax_lab_style = tax_lab_style(
		type = "text",
		max_angle = 90, 
		aspect_ratio = 1, 
		justify = "auto",
		size = 4)) +
	scale_shape_manual(values = c(1,16)) +
	scale_colour_manual(
		name = "Concentration",
		values = c("#C688D0")) +
	theme_bw() + ggtitle("Fungi - Gravel") +
	stat_ellipse(aes(linetype = Biochar), level = 0.95) +
	scale_linetype_manual(values = c(2,1)) +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		legend.text = element_text(size = 11))

plot.pca.its.p0.r <-
  subset_samples(rar.ps.meso.its,
	Pesticides == "Pesticides-" & Substrate == "Roots") %>%
  tax_transform("hellinger") %>%
  ord_calc(method = "PCA") %>%
    ord_plot(
	shape = "Biochar",
	colour = "#C688D0",
	size = 2,
    	tax_vec_length = 0.5,
	auto_caption = NA,
    	tax_lab_style = tax_lab_style(
		type = "text",
		max_angle = 90, 
		aspect_ratio = 1, 
		justify = "auto",
		size = 4)) +
	scale_shape_manual(values = c(1,16)) +
	scale_colour_manual(
		name = "Concentration",
		values = c("#C688D0")) +
	theme_bw() + ggtitle("Fungi - Roots") +
	stat_ellipse(aes(linetype = Biochar), level = 0.95) +
	scale_linetype_manual(values = c(2,1)) +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		legend.text = element_text(size = 11))

ggarrange(plot.pca.16s.p0.g, plot.pca.its.p0.g,
	    plot.pca.16s.p0.r, plot.pca.its.p0.r, 
	    ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
# ggsave("./Figures/pca_p0.png",
	width = 8, height = 8, unit = "in")
# ggsave("./Figures/Publication/FIG2.pdf",
	width = 190, height = 190, unit = "mm")

## Biochar with pesticides

plot.pca.16s.p.g <-
  subset_samples(rar.ps.meso.16s,
	Pesticides == "Pesticides+" & Substrate == "Gravel") %>%
  tax_transform("hellinger") %>%
  ord_calc(method = "PCA") %>%
    ord_plot(
	colour = "Pesticide_Concentration",
	shape = "Biochar",
	size = 2,
    	tax_vec_length = 0.5,
	auto_caption = NA,
    	tax_lab_style = tax_lab_style(
		type = "text",
		max_angle = 90, 
		aspect_ratio = 1, 
		justify = "auto",
		size = 4)) +
	scale_shape_manual(values = c(1,16)) +
	scale_colour_manual(
		name = "Concentration",
		values = c("#2FDC2A","#394EF5","#F90808")) +
	theme_bw() + ggtitle("Bacteria - Gravel") +
	stat_ellipse(aes(linetype = Biochar), level = 0.95) +
	scale_linetype_manual(values = c(2,1)) +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		legend.text = element_text(size = 11))

plot.pca.16s.p.r <- 
  subset_samples(rar.ps.meso.16s,
	Pesticides == "Pesticides+" & Substrate == "Roots") %>%
  tax_transform("hellinger") %>%
  ord_calc(method = "PCA") %>%
    ord_plot(
	colour = "Pesticide_Concentration",
	shape = "Biochar",
	size = 2,
    	tax_vec_length = 0.5,
	auto_caption = NA,
    	tax_lab_style = tax_lab_style(
		type = "text",
		max_angle = 90, 
		aspect_ratio = 1, 
		justify = "auto",
		size = 4)) + 
	scale_shape_manual(values = c(1,16)) +
	scale_colour_manual(
		name = "Concentration",
		values = c("#2FDC2A","#394EF5","#F90808")) +
	theme_bw() + ggtitle("Bacteria - Roots") +
	stat_ellipse(aes(linetype = Biochar), level = 0.95) +
	scale_linetype_manual(values = c(2,1)) +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		legend.text = element_text(size = 11))

plot.pca.its.p.g <-
  subset_samples(rar.ps.meso.its,
	Pesticides == "Pesticides+" & Substrate == "Gravel") %>%
  tax_transform("hellinger") %>%
  ord_calc(method = "PCA") %>%
    ord_plot(
	colour = "Pesticide_Concentration",
	shape = "Biochar",
	size = 2,
    	tax_vec_length = 0.5,
	auto_caption = NA,
    	tax_lab_style = tax_lab_style(
		type = "text",
		max_angle = 90, 
		aspect_ratio = 1, 
		justify = "auto",
		size = 4)) +
	scale_shape_manual(values = c(1,16)) +
	scale_colour_manual(
		name = "Concentration",
		values = c("#2FDC2A","#394EF5","#F90808")) +
	theme_bw() + ggtitle("Fungi - Gravel") +
	stat_ellipse(aes(linetype = Biochar), level = 0.95) +
	scale_linetype_manual(values = c(2,1)) +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		legend.text = element_text(size = 11))

plot.pca.its.p.r <- 
  subset_samples(rar.ps.meso.its,
	Pesticides == "Pesticides+" & Substrate == "Roots") %>%
  tax_transform("hellinger") %>%
  ord_calc(method = "PCA") %>%
    ord_plot(
	colour = "Pesticide_Concentration",
	shape = "Biochar",
	size = 2,
    	tax_vec_length = 0.5,
	auto_caption = NA,
    	tax_lab_style = tax_lab_style(
		type = "text",
		max_angle = 90, 
		aspect_ratio = 1, 
		justify = "auto",
		size = 4)) + 
	scale_shape_manual(values = c(1,16)) +
	scale_colour_manual(
		name = "Concentration",
		values = c("#2FDC2A","#394EF5","#F90808")) +
	theme_bw() + ggtitle("Fungi - Roots") +
	stat_ellipse(aes(linetype = Biochar), level = 0.95) +
	scale_linetype_manual(values = c(2,1)) +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		legend.text = element_text(size = 11))

ggarrange(plot.pca.16s.p.g, plot.pca.its.p.g,
	    plot.pca.16s.p.r, plot.pca.its.p.r, 
	    ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
# ggsave("./Figures/pca_p.png",
	width = 8, height = 8, unit = "in")
# ggsave("./Figures/Publication/FIG3.pdf",
	width = 190, height = 190, unit = "mm")

########## Beta-diversity analyses ##########

### Beta-diversity analyses will indicate whether the presence of biochar
### prevents or favours biotic homogenisation. We use similar permutation
### blocks to those used the PERMANOVAs, but without restricting permutations
### between mesocosms. These analyses were not included in the article.

### 16S

## Gravel

beta.16s.g <- betadisper(
 dist(
  decostand(
   otu_table(
    subset_samples(
      rar.ps.meso.16s, Substrate == "Gravel"
	)),
  method = "hellinger"),
 method = "euclidean"),
 group = sample_data(
  		subset_samples(
		rar.ps.meso.16s, 
		Substrate == "Gravel"))$Biochar,
 bias.adjust = TRUE
 )

set.seed(100)
perm.meso.16s.g <- with(
  meta(subset_samples(rar.ps.meso.16s, Substrate == "Gravel")), 
  how(
    nperm = 999,
    blocks = Pesticides,
    within = Within(type = "free")
  	)
   )
permutest(beta.16s.g, permutations = perm.meso.16s.g)
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.006771 0.0067713 3.2085    999  0.072 .
# Residuals 22 0.046429 0.0021104  

## Roots

beta.16s.r <- betadisper(
 dist(
  decostand(
   otu_table(
    subset_samples(
      rar.ps.meso.16s, Substrate == "Roots"
	)),
  method = "hellinger"),
 method = "euclidean"),
 group = sample_data(
  		subset_samples(
		rar.ps.meso.16s, 
		Substrate == "Roots"))$Biochar,
 bias.adjust = TRUE
 )

set.seed(100)
perm.meso.16s.r <- with(
  meta(subset_samples(rar.ps.meso.16s, Substrate == "Roots")), 
  how(
    nperm = 999,
    blocks = Pesticides,
    within = Within(type = "free")
  	)
   )
permutest(beta.16s.r, permutations = perm.meso.16s.r)
#           Df  Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.00445 0.0044502 0.6465    999  0.354
# Residuals 22 0.15143 0.0068832 

### ITS

## Gravel

beta.its.g <- betadisper(
 dist(
  decostand(
   otu_table(
    subset_samples(
      rar.ps.meso.its, Substrate == "Gravel"
	)),
  method = "hellinger"),
 method = "euclidean"),
 group = sample_data(
  		subset_samples(
		rar.ps.meso.its, 
		Substrate == "Gravel"))$Biochar,
 bias.adjust = TRUE
 )

set.seed(100)
perm.meso.its.g <- with(
  meta(subset_samples(rar.ps.meso.its, Substrate == "Gravel")), 
  how(
    nperm = 999,
    blocks = Pesticides,
    within = Within(type = "free")
  	)
   )
permutest(beta.its.g, permutations = perm.meso.its.g)
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.000177 0.0001767 0.0131    999  0.902
# Residuals 22 0.296652 0.0134842  

## Roots

beta.its.r <- betadisper(
 dist(
  decostand(
   otu_table(
    subset_samples(
      rar.ps.meso.its, Substrate == "Roots"
	)),
  method = "hellinger"),
 method = "euclidean"),
 group = sample_data(
  		subset_samples(
		rar.ps.meso.its, 
		Substrate == "Roots"))$Biochar,
 bias.adjust = TRUE
 )

set.seed(100)
perm.meso.its.r <- with(
  meta(subset_samples(rar.ps.meso.its, Substrate == "Roots")), 
  how(
    nperm = 999,
    blocks = Pesticides,
    within = Within(type = "free")
  	)
   )
permutest(beta.its.r, permutations = perm.meso.its.r)
#           Df  Sum Sq  Mean Sq  F N.Perm Pr(>F)
# Groups     1 0.00000 0.000000  0    999  0.994
# Residuals 22 0.99281 0.045128

### Figures

# png("./Figures/beta_div.png",
#	width = 800, height = 800)
par(mfrow = c(2,2))
plot(beta.16s.g, ellipse = TRUE, hull = FALSE,
	main = "Bacteria - Gravel",
	col = c("chocolate3","grey35"))
plot(beta.its.g, ellipse = TRUE, hull = FALSE,
	main = "Fungi - Gravel",
	col = c("chocolate3","grey35"))
plot(beta.16s.r, ellipse = TRUE, hull = FALSE,
	main = "Bacteria - Roots",
	col = c("chocolate3","grey35"))
plot(beta.its.r, ellipse = TRUE, hull = FALSE, xlim = c(-0.5,0.9),
	main = "Fungi - Roots",
	col = c("chocolate3","grey35"))
# dev.off()

########## Differential Abundance Analyses ##########

### We check the effect on biochar on the abundance of specific ASVs
### in the presence of pesticides. ZicoSeq requires using the
### original, non-rarified data. I set strata = Pesticides to confine
### permutations within each sampling period (before and after pesticides).

### 16S

## Gravel

set.seed(100)
zico.16s.g <- ZicoSeq(
	meta.dat = meta(subset_samples(ps.meso.16s, Substrate == "Gravel")), 
	feature.dat = t(
		otu_table(
			prune_taxa(
				taxa_sums(subset_samples(ps.meso.16s, Substrate == "Gravel")) > 0, 
					subset_samples(ps.meso.16s, Substrate == "Gravel")
				)
			)
		),
	grp.name = "Biochar",
	adj.name = c("Concentration","Depth"),
	feature.dat.type = "count",
	prev.filter = 0.2, mean.abund.filter = 0,  
                    max.abund.filter = 0.002, min.prop = 0, 
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                    # Use the square-root transformation
                    link.func = list(function (x) x^0.5), stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 999,  strata = "Pesticides", 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
# ZicoSeq.plot(zico.16s.g, pvalue.type = 'p.adj.fdr', cutoff = 0.1)
which(zico.16s.g$p.adj.fdr <= 0.1)

# We intersect these specific ASVs with the complete list of ASVs.
taxa.16s.g <- as.data.frame(
   tax_table(
	prune_taxa(
		taxa_sums(
			subset_samples(ps.meso.16s, Substrate == "Gravel")) > 0, 
			subset_samples(ps.meso.16s, Substrate == "Gravel")
			)
		)
	)
taxa.16s.g <- cbind(taxa.16s.g[intersect(rownames(taxa.16s.g), 
	rownames(as.data.frame(which(zico.16s.g$p.adj.fdr <= 0.1)))), ])
taxa.16s.g$ASV <- rownames(taxa.16s.g)
sign.16s.g <- t(sign(as.data.frame(zico.16s.g$coef.list)[4,]))
sign.16s.g_ASV <- rownames(sign.16s.g)
sign.16s.g <- as.data.frame(sign.16s.g, sign.16s.g_ASV)
sign.16s.g$ASV <- rownames(sign.16s.g)
colnames(sign.16s.g) <- c("Biochar","ASV")
taxa.16s.g <- merge(sign.16s.g, taxa.16s.g)
rownames(taxa.16s.g) <- taxa.16s.g$ASV
print(taxa.16s.g)

## Roots

set.seed(100)
zico.16s.r <- ZicoSeq(
	meta.dat = meta(subset_samples(ps.meso.16s, Substrate == "Roots")), 
	feature.dat = t(
		otu_table(
			prune_taxa(
				taxa_sums(subset_samples(ps.meso.16s, Substrate == "Roots")) > 0, 
					subset_samples(ps.meso.16s, Substrate == "Roots")
				)
			)
		),
	grp.name = "Biochar",
	adj.name = c("Concentration","Depth"),
	feature.dat.type = "count",
	prev.filter = 0.2, mean.abund.filter = 0,  
                    max.abund.filter = 0.002, min.prop = 0, 
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                    # Use the square-root transformation
                    link.func = list(function (x) x^0.5), stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 999,  strata = "Pesticides", 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
# ZicoSeq.plot(zico.16s.r, pvalue.type = 'p.adj.fdr', cutoff = 0.1)
which(zico.16s.r$p.adj.fdr <= 0.1)

# We intersect these specific ASVs with the complete list of ASVs.
taxa.16s.r <- as.data.frame(
   tax_table(
	prune_taxa(
		taxa_sums(
			subset_samples(ps.meso.16s, Substrate == "Roots")) > 0, 
			subset_samples(ps.meso.16s, Substrate == "Roots")
			)
		)
	)
taxa.16s.r <- cbind(taxa.16s.r[intersect(rownames(taxa.16s.r), 
	rownames(as.data.frame(which(zico.16s.r$p.adj.fdr <= 0.1)))), ])
taxa.16s.r$ASV <- rownames(taxa.16s.r)
sign.16s.r <- t(sign(as.data.frame(zico.16s.r$coef.list)[4,]))
sign.16s.r_ASV <- rownames(sign.16s.r)
sign.16s.r <- as.data.frame(sign.16s.r, sign.16s.r_ASV)
sign.16s.r$ASV <- rownames(sign.16s.r)
colnames(sign.16s.r) <- c("Biochar","ASV")
taxa.16s.r <- merge(sign.16s.r, taxa.16s.r)
rownames(taxa.16s.r) <- taxa.16s.r$ASV
print(taxa.16s.r)

## Combine gravel and roots dataframes.

taxa.16s.diff <- dplyr::union(taxa.16s.g, taxa.16s.r)
taxa.16s.diff$ASV <- gsub(pattern = "ASV", replacement = "", 
	x = taxa.16s.diff$ASV) %>% as.numeric()
taxa.16s.diff <- taxa.16s.diff[order(taxa.16s.diff$ASV),]
taxa.16s.diff <- taxa.16s.diff[,-c(1,3)]
taxa.16s.diff <- taxa.16s.diff[,c(2,3,4,5,6,1)]
taxa.16s.diff <- taxa.16s.diff[order(taxa.16s.diff$Biochar),]
# write.csv(taxa.16s.diff, "taxa_16s_diff.csv")

### ITS

## Gravel

set.seed(100)
zico.its.g <- ZicoSeq(
	meta.dat = meta(subset_samples(ps.meso.its, Substrate == "Gravel")), 
	feature.dat = t(
		otu_table(
			prune_taxa(
				taxa_sums(subset_samples(ps.meso.its, Substrate == "Gravel")) > 0, 
					subset_samples(ps.meso.its, Substrate == "Gravel")
				)
			)
		),
	grp.name = "Biochar",
	adj.name = c("Concentration","Depth"),
	feature.dat.type = "count",
	prev.filter = 0.2, mean.abund.filter = 0,  
                    max.abund.filter = 0.002, min.prop = 0, 
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                    # Use the square-root transformation
                    link.func = list(function (x) x^0.5), stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 999,  strata = "Pesticides", 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
# ZicoSeq.plot(zico.its.g, pvalue.type = 'p.adj.fdr', cutoff = 0.1)
which(zico.its.g$p.adj.fdr <= 0.1)

# We intersect these specific ASVs with the complete list of ASVs.
taxa.its.g <- as.data.frame(
   tax_table(
	prune_taxa(
		taxa_sums(
			subset_samples(ps.meso.its, Substrate == "Gravel")) > 0, 
			subset_samples(ps.meso.its, Substrate == "Gravel")
			)
		)
	)
taxa.its.g <- cbind(taxa.its.g[intersect(rownames(taxa.its.g), 
	rownames(as.data.frame(which(zico.its.g$p.adj.fdr <= 0.1)))), ])
taxa.its.g$ASV <- rownames(taxa.its.g)
sign.its.g <- t(sign(as.data.frame(zico.its.g$coef.list)[4,]))
sign.its.g_ASV <- rownames(sign.its.g)
sign.its.g <- as.data.frame(sign.its.g, sign.its.g_ASV)
sign.its.g$ASV <- rownames(sign.its.g)
colnames(sign.its.g) <- c("Biochar","ASV")
taxa.its.g <- merge(sign.its.g, taxa.its.g)
rownames(taxa.its.g) <- taxa.its.g$ASV
print(taxa.its.g)

## Roots

set.seed(100)
zico.its.r <- ZicoSeq(
	meta.dat = meta(subset_samples(ps.meso.its, Substrate == "Roots")), 
	feature.dat = t(
		otu_table(
			prune_taxa(
				taxa_sums(subset_samples(ps.meso.its, Substrate == "Roots")) > 0, 
					subset_samples(ps.meso.its, Substrate == "Roots")
				)
			)
		),
	grp.name = "Biochar",
	adj.name = c("Concentration","Depth"),
	feature.dat.type = "count",
	prev.filter = 0.2, mean.abund.filter = 0,  
                    max.abund.filter = 0.002, min.prop = 0, 
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                    # Use the square-root transformation
                    link.func = list(function (x) x^0.5), stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 999,  strata = "Pesticides", 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
# ZicoSeq.plot(zico.its.r, pvalue.type = 'p.adj.fdr', cutoff = 0.1)
which(zico.its.r$p.adj.fdr <= 0.1)

# We intersect these specific ASVs with the complete list of ASVs.
taxa.its.r <- as.data.frame(
   tax_table(
	prune_taxa(
		taxa_sums(
			subset_samples(ps.meso.its, Substrate == "Roots")) > 0, 
			subset_samples(ps.meso.its, Substrate == "Roots")
			)
		)
	)
taxa.its.r <- cbind(taxa.its.r[intersect(rownames(taxa.its.r), 
	rownames(as.data.frame(which(zico.its.r$p.adj.fdr <= 0.1)))), ])
taxa.its.r$ASV <- rownames(taxa.its.r)
sign.its.r <- t(sign(as.data.frame(zico.its.r$coef.list)[4,]))
sign.its.r_ASV <- rownames(sign.its.r)
sign.its.r <- as.data.frame(sign.its.r, sign.its.r_ASV)
sign.its.r$ASV <- rownames(sign.its.r)
colnames(sign.its.r) <- c("Biochar","ASV")
taxa.its.r <- merge(sign.its.r, taxa.its.r)
rownames(taxa.its.r) <- taxa.its.r$ASV
print(taxa.its.r)

## Combine gravel and roots dataframes.

taxa.its.diff <- dplyr::union(taxa.its.g, taxa.its.r)
taxa.its.diff$ASV <- gsub(pattern = "ASV", replacement = "", 
	x = taxa.its.diff$ASV) %>% as.numeric()
taxa.its.diff <- taxa.its.diff[order(taxa.its.diff$ASV),]
taxa.its.diff <- taxa.its.diff[,-c(1,3)]
taxa.its.diff <- taxa.its.diff[,c(2,3,4,5,6,7,1)]
taxa.its.diff <- taxa.its.diff[order(taxa.its.diff$Biochar),]
# write.csv(taxa.its.diff, "taxa_its_diff.csv")

### Figures

## 16S - Gravel

df.zico.16s.g <- data.frame(
	pvals = zico.16s.g$p.adj.fdr, 
	prevalence = apply(zico.16s.g$feature.dat, 1, function(x) mean(x > 0)), 
	abundance = apply(t(t(zico.16s.g$feature.dat)/colSums(zico.16s.g$feature.dat)), 
		1, function(x) mean(x)), 
	R2 = (zico.16s.g$R2 + runif(length(zico.16s.g$R2), 0, 1e-10)) * 
		t(sign(as.data.frame(zico.16s.g$coef.list)[4,])),
	taxa = rownames(zico.16s.g$R2)
	)
colnames(df.zico.16s.g) = c("p.adj.fdr","Prevalence","Abundance","R2","taxa")
df.zico.16s.g[df.zico.16s.g$p.adj.fdr > 0.1, 'taxa'] <- ''

plot.zico.16s.g <- 
  ggplot(df.zico.16s.g,
   aes(
	x = R2, 
	y = -log10(p.adj.fdr),
	size = Abundance,
	color = Prevalence)) +
	geom_point() +
	geom_vline(
		aes(xintercept = 0), color = 'gray', linetype = 'dashed') +
	geom_hline(
		aes(yintercept = -log10(0.1)), color = 'gray', linetype = 'dashed') +
	scale_colour_gradient2(low = "white", high = "#006D2C") +
	scale_y_continuous(limits = c(0, max(-log10(df.zico.16s.g$p.adj.fdr)) * 1.2)) +
	ggrepel::geom_text_repel(
		aes(label = taxa), max.overlaps = Inf, color = 'black', size = 3) +
	labs(x = bquote(R^2)) +
	ggtitle("Bacteria - Gravel") +
	theme_bw() +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank())

## 16S - Roots

df.zico.16s.r <- data.frame(
	pvals = zico.16s.r$p.adj.fdr, 
	prevalence = apply(zico.16s.r$feature.dat, 1, function(x) mean(x > 0)), 
	abundance = apply(t(t(zico.16s.r$feature.dat)/colSums(zico.16s.r$feature.dat)), 
		1, function(x) mean(x)), 
	R2 = (zico.16s.r$R2 + runif(length(zico.16s.r$R2), 0, 1e-10)) * 
		t(sign(as.data.frame(zico.16s.r$coef.list)[4,])),
	taxa = rownames(zico.16s.r$R2)
	)
colnames(df.zico.16s.r) = c("p.adj.fdr","Prevalence","Abundance","R2","taxa")
df.zico.16s.r[df.zico.16s.r$p.adj.fdr > 0.1, 'taxa'] <- ''

plot.zico.16s.r <- 
  ggplot(df.zico.16s.r,
   aes(
	x = R2, 
	y = -log10(p.adj.fdr),
	size = Abundance,
	color = Prevalence)) +
	geom_point() +
	geom_vline(
		aes(xintercept = 0), color = 'gray', linetype = 'dashed') +
	geom_hline(
		aes(yintercept = -log10(0.1)), color = 'gray', linetype = 'dashed') +
	scale_colour_gradient2(low = "white", high = "#006D2C") +
	scale_y_continuous(limits = c(0, max(-log10(df.zico.16s.r$p.adj.fdr)) * 1.2)) +
	ggrepel::geom_text_repel(
		aes(label = taxa), max.overlaps = Inf, color = 'black', size = 3) +
	labs(x = bquote(R^2)) +
	ggtitle("Bacteria - Roots") +
	theme_bw() +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank())

## ITS - Gravel

df.zico.its.g <- data.frame(
	pvals = zico.its.g$p.adj.fdr, 
	prevalence = apply(zico.its.g$feature.dat, 1, function(x) mean(x > 0)), 
	abundance = apply(t(t(zico.its.g$feature.dat)/colSums(zico.its.g$feature.dat)), 
		1, function(x) mean(x)), 
	R2 = (zico.its.g$R2 + runif(length(zico.its.g$R2), 0, 1e-10)) * 
		t(sign(as.data.frame(zico.its.g$coef.list)[4,])),
	taxa = rownames(zico.its.g$R2)
	)
colnames(df.zico.its.g) = c("p.adj.fdr","Prevalence","Abundance","R2","taxa")
df.zico.its.g[df.zico.its.g$p.adj.fdr > 0.1, 'taxa'] <- ''

plot.zico.its.g <- 
  ggplot(df.zico.its.g,
   aes(
	x = R2, 
	y = -log10(p.adj.fdr),
	size = Abundance,
	color = Prevalence)) +
	geom_point() +
	geom_vline(
		aes(xintercept = 0), color = 'gray', linetype = 'dashed') +
	geom_hline(
		aes(yintercept = -log10(0.1)), color = 'gray', linetype = 'dashed') +
	scale_colour_gradient2(low = "white", high = "#006D2C") +
	scale_y_continuous(limits = c(0, max(-log10(df.zico.its.g$p.adj.fdr)) * 1.2)) +
	ggrepel::geom_text_repel(
		aes(label = taxa), max.overlaps = Inf, color = 'black', size = 3) +
	labs(x = bquote(R^2)) +
	ggtitle("Fungi - Gravel") +
	theme_bw() +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank())

## ITS - Roots

df.zico.its.r <- data.frame(
	pvals = zico.its.r$p.adj.fdr, 
	prevalence = apply(zico.its.r$feature.dat, 1, function(x) mean(x > 0)), 
	abundance = apply(t(t(zico.its.r$feature.dat)/colSums(zico.its.r$feature.dat)), 
		1, function(x) mean(x)), 
	R2 = (zico.its.r$R2 + runif(length(zico.its.r$R2), 0, 1e-10)) * 
		t(sign(as.data.frame(zico.its.r$coef.list)[4,])),
	taxa = rownames(zico.its.r$R2)
	)
colnames(df.zico.its.r) = c("p.adj.fdr","Prevalence","Abundance","R2","taxa")
df.zico.its.r[df.zico.its.r$p.adj.fdr > 0.1, 'taxa'] <- ''

plot.zico.its.r <- 
  ggplot(df.zico.its.r,
   aes(
	x = R2, 
	y = -log10(p.adj.fdr),
	size = Abundance,
	color = Prevalence)) +
	geom_point() +
	geom_vline(
		aes(xintercept = 0), color = 'gray', linetype = 'dashed') +
	geom_hline(
		aes(yintercept = -log10(0.1)), color = 'gray', linetype = 'dashed') +
	scale_colour_gradient2(low = "white", high = "#006D2C") +
	scale_y_continuous(limits = c(0, max(-log10(df.zico.its.r$p.adj.fdr)) * 1.2)) +
	ggrepel::geom_text_repel(
		aes(label = taxa), max.overlaps = Inf, color = 'black', size = 3) +
	labs(x = bquote(R^2)) +
	ggtitle("Fungi - Roots") +
	theme_bw() +
	theme(
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank())

ggarrange(plot.zico.16s.g, plot.zico.its.g,
	    plot.zico.16s.r, plot.zico.its.r, 
	    ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")
# ggsave("./Figures/zicoseq.png",
	width = 8, height = 8, unit = "in")
# ggsave("./Figures/Publication/FIG4.pdf",
	width = 190, height = 190, unit = "mm")

### I now want to identify microbial ASVs that are differentially abundant
### in response to pesticides. I will compare the relative abundances
### before and after pesticides application. Of course, these results are
### confounded by the normal microbial succession that occurs throughout the
### summer. Hence, not all ASVs singled out by the analysis will have
### responded specifically to pesticides. Nonetheless, this analysis will
### allow me to identify specific ASVs, to look them up in the literature,
### and to single out a few ASVs that are known to either decrease or increase
### in the presence of pesticides. I will set the biochar as strata, to
### restrict the comparisons between mesocosms with or without biochar.

### 16S

## Gravel

set.seed(100)
zico.16s.g.p <- ZicoSeq(
	meta.dat = meta(subset_samples(ps.meso.16s, Substrate == "Gravel")), 
	feature.dat = t(
		otu_table(
			prune_taxa(
				taxa_sums(subset_samples(ps.meso.16s, Substrate == "Gravel")) > 0, 
					subset_samples(ps.meso.16s, Substrate == "Gravel")
				)
			)
		),
	grp.name = "Pesticides",
	adj.name = c("Depth"),
	feature.dat.type = "count",
	prev.filter = 0.2, mean.abund.filter = 0,  
                    max.abund.filter = 0.002, min.prop = 0, 
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                    # Use the square-root transformation
                    link.func = list(function (x) x^0.5), stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 999,  strata = "Biochar", 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
# ZicoSeq.plot(zico.16s.g.p, pvalue.type = 'p.adj.fdr', cutoff = 0.1)
which(zico.16s.g.p$p.adj.fdr <= 0.1)

# We intersect these specific ASVs with the complete list of ASVs.
taxa.16s.g.p <- as.data.frame(
   tax_table(
	prune_taxa(
		taxa_sums(
			subset_samples(ps.meso.16s, Substrate == "Gravel")) > 0, 
			subset_samples(ps.meso.16s, Substrate == "Gravel")
			)
		)
	)
taxa.16s.g.p <- cbind(taxa.16s.g.p[intersect(rownames(taxa.16s.g.p), 
	rownames(as.data.frame(which(zico.16s.g.p$p.adj.fdr <= 0.1)))), ])
taxa.16s.g.p$ASV <- rownames(taxa.16s.g.p)
sign.16s.g.p <- t(sign(as.data.frame(zico.16s.g.p$coef.list)[3,]))
sign.16s.g.p_ASV <- rownames(sign.16s.g.p)
sign.16s.g.p <- as.data.frame(sign.16s.g.p, sign.16s.g.p_ASV)
sign.16s.g.p$ASV <- rownames(sign.16s.g.p)
colnames(sign.16s.g.p) <- c("Pesticides","ASV")
taxa.16s.g.p <- merge(sign.16s.g.p, taxa.16s.g.p)
rownames(taxa.16s.g.p) <- taxa.16s.g.p$ASV
print(taxa.16s.g.p)

## Roots

set.seed(100)
zico.16s.r.p <- ZicoSeq(
	meta.dat = meta(subset_samples(ps.meso.16s, Substrate == "Roots")), 
	feature.dat = t(
		otu_table(
			prune_taxa(
				taxa_sums(subset_samples(ps.meso.16s, Substrate == "Roots")) > 0, 
					subset_samples(ps.meso.16s, Substrate == "Roots")
				)
			)
		),
	grp.name = "Pesticides",
	adj.name = c("Depth"),
	feature.dat.type = "count",
	prev.filter = 0.2, mean.abund.filter = 0,  
                    max.abund.filter = 0.002, min.prop = 0, 
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                    # Use the square-root transformation
                    link.func = list(function (x) x^0.5), stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 999,  strata = "Biochar", 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
# ZicoSeq.plot(zico.16s.r.p, pvalue.type = 'p.adj.fdr', cutoff = 0.1)
which(zico.16s.r.p$p.adj.fdr <= 0.1)

# We intersect these specific ASVs with the complete list of ASVs.
taxa.16s.r.p <- as.data.frame(
   tax_table(
	prune_taxa(
		taxa_sums(
			subset_samples(ps.meso.16s, Substrate == "Roots")) > 0, 
			subset_samples(ps.meso.16s, Substrate == "Roots")
			)
		)
	)
taxa.16s.r.p <- cbind(taxa.16s.r.p[intersect(rownames(taxa.16s.r.p), 
	rownames(as.data.frame(which(zico.16s.r.p$p.adj.fdr <= 0.1)))), ])
taxa.16s.r.p$ASV <- rownames(taxa.16s.r.p)
sign.16s.r.p <- t(sign(as.data.frame(zico.16s.r.p$coef.list)[3,]))
sign.16s.r.p_ASV <- rownames(sign.16s.r.p)
sign.16s.r.p <- as.data.frame(sign.16s.r.p, sign.16s.r.p_ASV)
sign.16s.r.p$ASV <- rownames(sign.16s.r.p)
colnames(sign.16s.r.p) <- c("Pesticides","ASV")
taxa.16s.r.p <- merge(sign.16s.r.p, taxa.16s.r.p)
rownames(taxa.16s.r.p) <- taxa.16s.r.p$ASV
print(taxa.16s.r.p)

## Combine gravel and roots dataframes.

taxa.16s.diff.p <- dplyr::union(taxa.16s.g.p, taxa.16s.r.p)
taxa.16s.diff.p$ASV <- gsub(pattern = "ASV", replacement = "", 
	x = taxa.16s.diff.p$ASV) %>% as.numeric()
taxa.16s.diff.p <- taxa.16s.diff.p[order(taxa.16s.diff.p$ASV),]
taxa.16s.diff.p <- taxa.16s.diff.p[,-c(1,3)]
taxa.16s.diff.p <- taxa.16s.diff.p[,c(2,3,4,5,6,1)]
taxa.16s.diff.p <- taxa.16s.diff.p[order(taxa.16s.diff.p$Pesticides),]
# write.csv(taxa.16s.diff.p, "taxa_16s_diff_p.csv")

### ITS

## Gravel

set.seed(100)
zico.its.g.p <- ZicoSeq(
	meta.dat = meta(subset_samples(ps.meso.its, Substrate == "Gravel")), 
	feature.dat = t(
		otu_table(
			prune_taxa(
				taxa_sums(subset_samples(ps.meso.its, Substrate == "Gravel")) > 0, 
					subset_samples(ps.meso.its, Substrate == "Gravel")
				)
			)
		),
	grp.name = "Pesticides",
	adj.name = c("Depth"),
	feature.dat.type = "count",
	prev.filter = 0.2, mean.abund.filter = 0,  
                    max.abund.filter = 0.002, min.prop = 0, 
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                    # Use the square-root transformation
                    link.func = list(function (x) x^0.5), stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 999,  strata = "Biochar", 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
# ZicoSeq.plot(zico.its.g.p, pvalue.type = 'p.adj.fdr', cutoff = 0.1)
which(zico.its.g.p$p.adj.fdr <= 0.1)

# We intersect these specific ASVs with the complete list of ASVs.
taxa.its.g.p <- as.data.frame(
   tax_table(
	prune_taxa(
		taxa_sums(
			subset_samples(ps.meso.its, Substrate == "Gravel")) > 0, 
			subset_samples(ps.meso.its, Substrate == "Gravel")
			)
		)
	)
taxa.its.g.p <- cbind(taxa.its.g.p[intersect(rownames(taxa.its.g.p), 
	rownames(as.data.frame(which(zico.its.g.p$p.adj.fdr <= 0.1)))), ])
taxa.its.g.p$ASV <- rownames(taxa.its.g.p)
sign.its.g.p <- t(sign(as.data.frame(zico.its.g.p$coef.list)[3,]))
sign.its.g.p_ASV <- rownames(sign.its.g.p)
sign.its.g.p <- as.data.frame(sign.its.g.p, sign.its.g.p_ASV)
sign.its.g.p$ASV <- rownames(sign.its.g.p)
colnames(sign.its.g.p) <- c("Pesticides","ASV")
taxa.its.g.p <- merge(sign.its.g.p, taxa.its.g.p)
rownames(taxa.its.g.p) <- taxa.its.g.p$ASV
print(taxa.its.g.p)

## Roots

set.seed(100)
zico.its.r.p <- ZicoSeq(
	meta.dat = meta(subset_samples(ps.meso.its, Substrate == "Roots")), 
	feature.dat = t(
		otu_table(
			prune_taxa(
				taxa_sums(subset_samples(ps.meso.its, Substrate == "Roots")) > 0, 
					subset_samples(ps.meso.its, Substrate == "Roots")
				)
			)
		),
	grp.name = "Pesticides",
	adj.name = c("Depth"),
	feature.dat.type = "count",
	prev.filter = 0.2, mean.abund.filter = 0,  
                    max.abund.filter = 0.002, min.prop = 0, 
                    # Winsorization to replace outliers
                    is.winsor = TRUE, outlier.pct = 0.03, winsor.end = 'top',
                    # Use the square-root transformation
                    link.func = list(function (x) x^0.5), stats.combine.func = max,
                    # Permutation-based multiple testing correction
                    perm.no = 999,  strata = "Biochar", 
                    # Reference-based multiple stage normalization
                    ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                    # Family-wise error rate control
                    is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)
# ZicoSeq.plot(zico.its.r.p, pvalue.type = 'p.adj.fdr', cutoff = 0.1)
which(zico.its.r.p$p.adj.fdr <= 0.1)

# We intersect these specific ASVs with the complete list of ASVs.
taxa.its.r.p <- as.data.frame(
   tax_table(
	prune_taxa(
		taxa_sums(
			subset_samples(ps.meso.its, Substrate == "Roots")) > 0, 
			subset_samples(ps.meso.its, Substrate == "Roots")
			)
		)
	)
taxa.its.r.p <- cbind(taxa.its.r.p[intersect(rownames(taxa.its.r.p), 
	rownames(as.data.frame(which(zico.its.r.p$p.adj.fdr <= 0.1)))), ])
taxa.its.r.p$ASV <- rownames(taxa.its.r.p)
sign.its.r.p <- t(sign(as.data.frame(zico.its.r.p$coef.list)[3,]))
sign.its.r.p_ASV <- rownames(sign.its.r.p)
sign.its.r.p <- as.data.frame(sign.its.r.p, sign.its.r.p_ASV)
sign.its.r.p$ASV <- rownames(sign.its.r.p)
colnames(sign.its.r.p) <- c("Pesticides","ASV")
taxa.its.r.p <- merge(sign.its.r.p, taxa.its.r.p)
rownames(taxa.its.r.p) <- taxa.its.r.p$ASV
print(taxa.its.r.p)

## Combine gravel and roots dataframes.

taxa.its.diff.p <- dplyr::union(taxa.its.g.p, taxa.its.r.p)
taxa.its.diff.p$ASV <- gsub(pattern = "ASV", replacement = "", 
	x = taxa.its.diff.p$ASV) %>% as.numeric()
taxa.its.diff.p <- taxa.its.diff.p[order(taxa.its.diff.p$ASV),]
taxa.its.diff.p <- taxa.its.diff.p[,-c(1,3)]
taxa.its.diff.p <- taxa.its.diff.p[,c(2,3,4,5,6,7,1)]
taxa.its.diff.p <- taxa.its.diff.p[order(taxa.its.diff.p$Pesticides),]
# write.csv(taxa.its.diff.p, "taxa_its_diff_p.csv")

########## Chi-Square Tests ##########

### We now compare the mean relative abundances of the phyla with or without
### biochar, for each substrate and depth, with a Chi-squared test.

### 16S

## Before pesticides
## Gravel

rar.ps.meso.16s.p0.g.phyla <- tax_glom(
	subset_samples(rar.ps.meso.16s, 
		Pesticides == "Pesticides-" & Substrate == "Gravel"), "Phylum")

rar.ps.meso.16s.p0.g.phyla.15 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.16s.p0.g.phyla, Depth == "0-15 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.16s.p0.g.phyla.15 <- prune_taxa(
	taxa_sums(rar.ps.meso.16s.p0.g.phyla.15) > 0, 
	rar.ps.meso.16s.p0.g.phyla.15)
rar.ps.meso.16s.p0.g.phyla.15.otu <- otu_table(rar.ps.meso.16s.p0.g.phyla.15)
colnames(rar.ps.meso.16s.p0.g.phyla.15.otu) <- tax_table(rar.ps.meso.16s.p0.g.phyla.15)[,2]
chisq.test(rar.ps.meso.16s.p0.g.phyla.15.otu)
# X-squared = 168.06, df = 21, p-value < 2.2e-16

rar.ps.meso.16s.p0.g.phyla.30 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.16s.p0.g.phyla, Depth == "15-30 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.16s.p0.g.phyla.30 <- prune_taxa(
	taxa_sums(rar.ps.meso.16s.p0.g.phyla.30) > 0, 
	rar.ps.meso.16s.p0.g.phyla.30)
rar.ps.meso.16s.p0.g.phyla.30.otu <- otu_table(rar.ps.meso.16s.p0.g.phyla.30)
colnames(rar.ps.meso.16s.p0.g.phyla.30.otu) <- tax_table(rar.ps.meso.16s.p0.g.phyla.30)[,2]
chisq.test(rar.ps.meso.16s.p0.g.phyla.30.otu)
# X-squared = 181.36, df = 15, p-value < 2.2e-16

## Before pesticides
## Roots

rar.ps.meso.16s.p0.r.phyla <- tax_glom(
	subset_samples(rar.ps.meso.16s, 
		Pesticides == "Pesticides-" & Substrate == "Roots"), "Phylum")

rar.ps.meso.16s.p0.r.phyla.15 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.16s.p0.r.phyla, Depth == "0-15 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.16s.p0.r.phyla.15 <- prune_taxa(
	taxa_sums(rar.ps.meso.16s.p0.r.phyla.15) > 0, 
	rar.ps.meso.16s.p0.r.phyla.15)
rar.ps.meso.16s.p0.r.phyla.15.otu <- otu_table(rar.ps.meso.16s.p0.r.phyla.15)
colnames(rar.ps.meso.16s.p0.r.phyla.15.otu) <- tax_table(rar.ps.meso.16s.p0.r.phyla.15)[,2]
chisq.test(rar.ps.meso.16s.p0.r.phyla.15.otu)
# X-squared = 296.93, df = 16, p-value < 2.2e-16

rar.ps.meso.16s.p0.r.phyla.30 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.16s.p0.r.phyla, Depth == "15-30 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.16s.p0.r.phyla.30 <- prune_taxa(
	taxa_sums(rar.ps.meso.16s.p0.r.phyla.30) > 0, 
	rar.ps.meso.16s.p0.r.phyla.30)
rar.ps.meso.16s.p0.r.phyla.30.otu <- otu_table(rar.ps.meso.16s.p0.r.phyla.30)
colnames(rar.ps.meso.16s.p0.r.phyla.30.otu) <- tax_table(rar.ps.meso.16s.p0.r.phyla.30)[,2]
chisq.test(rar.ps.meso.16s.p0.r.phyla.30.otu)
# X-squared = 44.203, df = 11, p-value = 6.695e-06

## After pesticides
## Gravel

rar.ps.meso.16s.p.g.phyla <- tax_glom(
	subset_samples(rar.ps.meso.16s, 
		Pesticides == "Pesticides+" & Substrate == "Gravel"), "Phylum")

rar.ps.meso.16s.p.g.phyla.15 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.16s.p.g.phyla, Depth == "0-15 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.16s.p.g.phyla.15 <- prune_taxa(
	taxa_sums(rar.ps.meso.16s.p.g.phyla.15) > 0, 
	rar.ps.meso.16s.p.g.phyla.15)
rar.ps.meso.16s.p.g.phyla.15.otu <- otu_table(rar.ps.meso.16s.p.g.phyla.15)
colnames(rar.ps.meso.16s.p.g.phyla.15.otu) <- tax_table(rar.ps.meso.16s.p.g.phyla.15)[,2]
chisq.test(rar.ps.meso.16s.p.g.phyla.15.otu)
# X-squared = 92.387, df = 23, p-value = 2.816e-10

rar.ps.meso.16s.p.g.phyla.30 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.16s.p.g.phyla, Depth == "15-30 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.16s.p.g.phyla.30 <- prune_taxa(
	taxa_sums(rar.ps.meso.16s.p.g.phyla.30) > 0, 
	rar.ps.meso.16s.p.g.phyla.30)
rar.ps.meso.16s.p.g.phyla.30.otu <- otu_table(rar.ps.meso.16s.p.g.phyla.30)
colnames(rar.ps.meso.16s.p.g.phyla.30.otu) <- tax_table(rar.ps.meso.16s.p.g.phyla.30)[,2]
chisq.test(rar.ps.meso.16s.p.g.phyla.30.otu)
# X-squared = 106.94, df = 24, p-value = 1.923e-12

## After pesticides
## Roots

rar.ps.meso.16s.p.r.phyla <- tax_glom(
	subset_samples(rar.ps.meso.16s, 
		Pesticides == "Pesticides+" & Substrate == "Roots"), "Phylum")

rar.ps.meso.16s.p.r.phyla.15 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.16s.p.r.phyla, Depth == "0-15 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.16s.p.r.phyla.15 <- prune_taxa(
	taxa_sums(rar.ps.meso.16s.p.r.phyla.15) > 0, 
	rar.ps.meso.16s.p.r.phyla.15)
rar.ps.meso.16s.p.r.phyla.15.otu <- otu_table(rar.ps.meso.16s.p.r.phyla.15)
colnames(rar.ps.meso.16s.p.r.phyla.15.otu) <- tax_table(rar.ps.meso.16s.p.r.phyla.15)[,2]
chisq.test(rar.ps.meso.16s.p.r.phyla.15.otu)
# X-squared = 125.49, df = 19, p-value < 2.2e-16

rar.ps.meso.16s.p.r.phyla.30 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.16s.p.r.phyla, Depth == "15-30 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.16s.p.r.phyla.30 <- prune_taxa(
	taxa_sums(rar.ps.meso.16s.p.r.phyla.30) > 0, 
	rar.ps.meso.16s.p.r.phyla.30)
rar.ps.meso.16s.p.r.phyla.30.otu <- otu_table(rar.ps.meso.16s.p.r.phyla.30)
colnames(rar.ps.meso.16s.p.r.phyla.30.otu) <- tax_table(rar.ps.meso.16s.p.r.phyla.30)[,2]
chisq.test(rar.ps.meso.16s.p.r.phyla.30.otu)
# X-squared = 175.84, df = 21, p-value < 2.2e-16

### ITS

## Before pesticides
## Gravel

rar.ps.meso.its.p0.g.phyla <- tax_glom(
	subset_samples(rar.ps.meso.its, 
		Pesticides == "Pesticides-" & Substrate == "Gravel"), "Phylum")

rar.ps.meso.its.p0.g.phyla.15 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.its.p0.g.phyla, Depth == "0-15 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.its.p0.g.phyla.15 <- prune_taxa(
	taxa_sums(rar.ps.meso.its.p0.g.phyla.15) > 0, 
	rar.ps.meso.its.p0.g.phyla.15)
rar.ps.meso.its.p0.g.phyla.15.otu <- otu_table(rar.ps.meso.its.p0.g.phyla.15)
colnames(rar.ps.meso.its.p0.g.phyla.15.otu) <- tax_table(rar.ps.meso.its.p0.g.phyla.15)[,2]
chisq.test(rar.ps.meso.its.p0.g.phyla.15.otu)
# X-squared = 268.79, df = 6, p-value < 2.2e-16

rar.ps.meso.its.p0.g.phyla.30 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.its.p0.g.phyla, Depth == "15-30 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.its.p0.g.phyla.30 <- prune_taxa(
	taxa_sums(rar.ps.meso.its.p0.g.phyla.30) > 0, 
	rar.ps.meso.its.p0.g.phyla.30)
rar.ps.meso.its.p0.g.phyla.30.otu <- otu_table(rar.ps.meso.its.p0.g.phyla.30)
colnames(rar.ps.meso.its.p0.g.phyla.30.otu) <- tax_table(rar.ps.meso.its.p0.g.phyla.30)[,2]
chisq.test(rar.ps.meso.its.p0.g.phyla.30.otu)
# X-squared = 703.52, df = 6, p-value < 2.2e-16

## Before pesticides
## Roots

rar.ps.meso.its.p0.r.phyla <- tax_glom(
	subset_samples(rar.ps.meso.its, 
		Pesticides == "Pesticides-" & Substrate == "Roots"), "Phylum")

rar.ps.meso.its.p0.r.phyla.15 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.its.p0.r.phyla, Depth == "0-15 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.its.p0.r.phyla.15 <- prune_taxa(
	taxa_sums(rar.ps.meso.its.p0.r.phyla.15) > 0, 
	rar.ps.meso.its.p0.r.phyla.15)
rar.ps.meso.its.p0.r.phyla.15.otu <- otu_table(rar.ps.meso.its.p0.r.phyla.15)
colnames(rar.ps.meso.its.p0.r.phyla.15.otu) <- tax_table(rar.ps.meso.its.p0.r.phyla.15)[,2]
chisq.test(rar.ps.meso.its.p0.r.phyla.15.otu)
# X-squared = 364.24, df = 6, p-value < 2.2e-16

rar.ps.meso.its.p0.r.phyla.30 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.its.p0.r.phyla, Depth == "15-30 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.its.p0.r.phyla.30 <- prune_taxa(
	taxa_sums(rar.ps.meso.its.p0.r.phyla.30) > 0, 
	rar.ps.meso.its.p0.r.phyla.30)
rar.ps.meso.its.p0.r.phyla.30.otu <- otu_table(rar.ps.meso.its.p0.r.phyla.30)
colnames(rar.ps.meso.its.p0.r.phyla.30.otu) <- tax_table(rar.ps.meso.its.p0.r.phyla.30)[,2]
chisq.test(rar.ps.meso.its.p0.r.phyla.30.otu)
# X-squared = 189.39, df = 5, p-value < 2.2e-16

## After pesticides
## Gravel

rar.ps.meso.its.p.g.phyla <- tax_glom(
	subset_samples(rar.ps.meso.its, 
		Pesticides == "Pesticides+" & Substrate == "Gravel"), "Phylum")

rar.ps.meso.its.p.g.phyla.15 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.its.p.g.phyla, Depth == "0-15 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.its.p.g.phyla.15 <- prune_taxa(
	taxa_sums(rar.ps.meso.its.p.g.phyla.15) > 0, 
	rar.ps.meso.its.p.g.phyla.15)
rar.ps.meso.its.p.g.phyla.15.otu <- otu_table(rar.ps.meso.its.p.g.phyla.15)
colnames(rar.ps.meso.its.p.g.phyla.15.otu) <- tax_table(rar.ps.meso.its.p.g.phyla.15)[,2]
chisq.test(rar.ps.meso.its.p.g.phyla.15.otu)
# X-squared = 179.02, df = 6, p-value < 2.2e-16

rar.ps.meso.its.p.g.phyla.30 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.its.p.g.phyla, Depth == "15-30 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.its.p.g.phyla.30 <- prune_taxa(
	taxa_sums(rar.ps.meso.its.p.g.phyla.30) > 0, 
	rar.ps.meso.its.p.g.phyla.30)
rar.ps.meso.its.p.g.phyla.30.otu <- otu_table(rar.ps.meso.its.p.g.phyla.30)
colnames(rar.ps.meso.its.p.g.phyla.30.otu) <- tax_table(rar.ps.meso.its.p.g.phyla.30)[,2]
chisq.test(rar.ps.meso.its.p.g.phyla.30.otu)
# X-squared = 204.68, df = 7, p-value < 2.2e-16

## After pesticides
## Roots

rar.ps.meso.its.p.r.phyla <- tax_glom(
	subset_samples(rar.ps.meso.its, 
		Pesticides == "Pesticides+" & Substrate == "Roots"), "Phylum")

rar.ps.meso.its.p.r.phyla.15 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.its.p.r.phyla, Depth == "0-15 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.its.p.r.phyla.15 <- prune_taxa(
	taxa_sums(rar.ps.meso.its.p.r.phyla.15) > 0, 
	rar.ps.meso.its.p.r.phyla.15)
rar.ps.meso.its.p.r.phyla.15.otu <- otu_table(rar.ps.meso.its.p.r.phyla.15)
colnames(rar.ps.meso.its.p.r.phyla.15.otu) <- tax_table(rar.ps.meso.its.p.r.phyla.15)[,2]
chisq.test(rar.ps.meso.its.p.r.phyla.15.otu)
# X-squared = 180.41, df = 7, p-value < 2.2e-16

rar.ps.meso.its.p.r.phyla.30 <- phyloseq::merge_samples(
	subset_samples(rar.ps.meso.its.p.r.phyla, Depth == "15-30 cm"),
	group = "Biochar", fun = "mean")
rar.ps.meso.its.p.r.phyla.30 <- prune_taxa(
	taxa_sums(rar.ps.meso.its.p.r.phyla.30) > 0, 
	rar.ps.meso.its.p.r.phyla.30)
rar.ps.meso.its.p.r.phyla.30.otu <- otu_table(rar.ps.meso.its.p.r.phyla.30)
colnames(rar.ps.meso.its.p.r.phyla.30.otu) <- tax_table(rar.ps.meso.its.p.r.phyla.30)[,2]
chisq.test(rar.ps.meso.its.p.r.phyla.30.otu)
# X-squared = 164.55, df = 7, p-value < 2.2e-16

### Figures

## First, we create colour palettes to use throughout the graphs.
## We will assign a specific colour for each phyla that we expect to
## see in the graphs. This will ensure that each phylum has its own colour.

myPalettePhyla16S <- tax_palette(
  data = rar.ps.meso.16s, rank = "Phylum", n = 25, pal = "brewerPlus",
  add = c(Other = "white")
	)
# tax_palette_plot(myPalettePhyla16S)

myPalettePhylaITS <- tax_palette(
  data = rar.ps.meso.its, rank = "Phylum", n = 25, pal = "brewerPlus",
  add = c(Other = "white")
	)
# tax_palette_plot(myPalettePhylaITS)

### 16S

## Biochar without pesticides

comp_barplot(subset_samples(rar.ps.meso.16s, Pesticides == "Pesticides-"),
	tax_level = "Phylum",
	n_taxa = 10,
	palette = myPalettePhyla16S,
	label = "Concentration_Biochar",
	sample_order = rev(rownames(sample_data(subset_samples(
				rar.ps.meso.16s, Pesticides == "Pesticides-")))),
	facet_by = c("Substrate", "Depth"),
		) + coord_flip()
# ggsave("./Figures/16s_comp_p0_phyla.png",
	width = 9, height = 6, units = "in")

## Biochar with pesticides

comp_barplot(subset_samples(rar.ps.meso.16s, Pesticides == "Pesticides+"),
	tax_level = "Phylum",
	n_taxa = 10,
	palette = myPalettePhyla16S,
	label = "Concentration_Biochar",
	sample_order = rev(rownames(sample_data(subset_samples(
				rar.ps.meso.16s, Pesticides == "Pesticides+")))),
	facet_by = c("Substrate", "Depth"),
		) + coord_flip()
# ggsave("./Figures/16s_comp_p_phyla.png",
	width = 9, height = 6, units = "in")

### ITS

## Biochar without pesticides

comp_barplot(subset_samples(rar.ps.meso.its, Pesticides == "Pesticides-"),
	tax_level = "Phylum",
	n_taxa = 10,
	palette = myPalettePhylaITS,
	label = "Concentration_Biochar",
	sample_order = rev(rownames(sample_data(subset_samples(
				rar.ps.meso.its, Pesticides == "Pesticides-")))),
	facet_by = c("Substrate", "Depth"),
		) + coord_flip()
# ggsave("./Figures/its_comp_p0_phyla.png",
	width = 9, height = 6, units = "in")

## Biochar with pesticides

comp_barplot(subset_samples(rar.ps.meso.its, Pesticides == "Pesticides+"),
	tax_level = "Phylum",
	n_taxa = 10,
	palette = myPalettePhylaITS,
	label = "Concentration_Biochar",
	sample_order = rev(rownames(sample_data(subset_samples(
				rar.ps.meso.its, Pesticides == "Pesticides+")))),
	facet_by = c("Substrate", "Depth"),
		) + coord_flip()
# ggsave("./Figures/its_comp_p_phyla.png",
	width = 9, height = 6, units = "in")

#############################################################################