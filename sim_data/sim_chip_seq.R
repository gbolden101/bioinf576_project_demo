#Tested on R 4.4.0
# Simulated Data from the ChIPsm documentation
#https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPsim/inst/doc/ChIPsimIntro.pdf

if(!require(dplyr)) install.packages("dplyr", dependencies=TRUE)
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("ChIPsim")
if(!require(purrr)) install.packages("purrr", dependencies=TRUE)

library(ChIPsim)
library(tidyverse)

set.seed(1)

#big note still not done because the referance genome needs to be mmu
#creating the reference genome is crucial for simulation step
#creating reference genome, controls the length of reference genome
#passed to the n function
chrLen <- c(2e5, 1e5)
# Paste DNA bases for n number of times collapse into no spaces
DNA_BASES <- c("A", "T", "G", "C")
chromosomes <- sapply(chrLen, function(n) paste(sample(DNA_BASES, n, replace = TRUE), collapse = ""))

#paste the CHR for each n  of chromosomes
names(chromosomes) <- paste("CHR", seq_along(chromosomes))
# makes a DNAString object out of the chromosomes
genome <- Biostrings::DNAStringSet(chromosomes)



randomQuality <- function(read){
  paste(sample(unlist(strsplit(rawToChar(as.raw(64:104)),"")), 
               nchar(read), replace = TRUE), collapse="")
}

dfReads <- function(readPos, readNames, sequence, readLen, ...){
  
  ## create vector to hold read sequences and qualities
  readSeq <- character(sum(sapply(readPos, sapply, length)))
  readQual <- character(sum(sapply(readPos, sapply, length)))
  
  idx <- 1
  ## process read positions for each chromosome and strand
  for(k in length(readPos)){ ## chromosome
    for(i in 1:2){ ## strand
      for(j in 1:length(readPos[[k]][[i]])){
        ## get (true) sequence
        readSeq[idx] <- as.character(readSequence(readPos[[k]][[i]][j], sequence[[k]], 
                                                  strand=ifelse(i==1, 1, -1), readLen=readLen))
        ## get quality
        readQual[idx] <- randomQuality(readSeq[idx])
        ## introduce sequencing errors
        readSeq[idx] <- readError(readSeq[idx], decodeQuality(readQual[idx]))
        idx <- idx + 1
      }
    }
  }
  data.frame(name=unlist(readNames), sequence=readSeq, quality=readQual, 
             stringsAsFactors = FALSE)
}
###########################

# use default function to obtain the default 
#list of functions and replace the one for the 
#final step of the simulation

# In this case we are changing the mean 
#frequent length of 150bp and the default
myFunctions <- defaultFunctions()

myFunctions$readSequence <- dfReads 
nReads <- 100

simulated <- simChIP(nReads, genome, 
                     file = "test", functions = myFunctions, 
                     control = defaultControl(readDensity=list(meanLength = 150)))
######################################
simmed_reads <- simulated$readSequence %>% 
                  tibble() %>% 
                dplyr::select(sequence) %>% 
                filter(sequence != "" & !is.na(sequence)) %>% 
                  as.list()

simmed_qscores <- simulated$readSequence %>% 
                  tibble() %>% 
                  dplyr::select(quality) %>% 
                  filter(quality != "" & !is.na(quality)) %>% 
                  as.list()

simmed_names <- tibble(name = simulated$readSequence$name[1:length(unlist(simmed_qscores))]) %>% 
                as.list() %>% 
                lapply(function(x) gsub("^(.*)", "@\\1", x))
#########################################
library(purrr)

# Create structured reads list using pmap()
reads <- pmap(list(simmed_names$name, simmed_reads$sequence, simmed_qscores$quality), 
              function(name, seq, quality) {
                list(id = name, seq = seq, plus = "+", quality = quality)
              })

# Print FASTQ format
lines <- walk(reads, function(read) {
  cat(read$id, "\n", read$seq, "\n", read$plus, "\n", read$quality, "\n\n", sep="")
})

file_name <- "sim_reads"

fastq_entry <- unlist(lines)

writeLines(fastq_entry, file_name)



