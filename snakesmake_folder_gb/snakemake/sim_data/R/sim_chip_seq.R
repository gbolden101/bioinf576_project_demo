#Tested on R 4.4.0
# Simulated Data from the ChIPsm documentation
#https://www.bioconductor.org/packages/devel/bioc/vignettes/ChIPsim/inst/doc/ChIPsimIntro.pdf
#https://github.com/humburg/ChIPsim/tree/master

if(!require(dplyr)) install.packages("dplyr", dependencies=TRUE)
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("ChIPsim")
if(!require(purrr)) install.packages("purrr", dependencies=TRUE)
if(!require(purrr)) install.packages("actuar", dependencies=TRUE)
if(!require(testthat)) install.packages("testthat", dependencies=TRUE)

library(testthat)
library(ChIPsim)
library(tidyverse)

set.seed(1)

#creating the reference genome is crucial for simulation step
#creating reference genome, controls the length of reference genome
#passed to the n function

fasta_file <- '../inputs/chr1_limit.fa'

fasta = readDNAStringSet(fasta_file)

genome <- fasta

#need to add in background and binding features for macs2 to work
 
#Binding site sim divides the genome into a # of non-overlapping binding sites.

#Simulation of BngSites is driven by a markov that with two states BngSites, and background respectively

#lets start by specifying the transition probabilities
transition <- list(Binding=c(Background=1), Background=c(Binding=0.05, Background=0.95))

# Q: why include background regions?
# A: A little more realistic to experimental 

#which is simply a numeric vector with the no-zero transition probabilities to other states
#
transition <- lapply(transition, "class<-","StateDistribution")

#We will not allow the sequence to start with a binding site
#the 'stateDistribution'
initial <- c(Binding=0, Background=1)
#assign 
class(initial) <- "StateDistribution"

# define fxn to generate parameters for binding and background region--> (start, length), sampling weight

#arguments: start_pos , 
#           length of region ,
#           shape:
#           scale:

#modeling background
backgroundFeature <- function(start, length=1000, shape=1, scale=20){
  #weight of the background feature is generated randomly adding somevariability
  #rgamma sample from a gamma dist with parameters shape=1, and scale=20
  #gdist. is right skewed
  #mean = scale * shape, var = shape * scale^2
  #
  weight <- rgamma(1, shape=shape,   scale=scale)
  params <- list(start = start, length = length, weight = weight)
  class(params) <- c("Background", "SimulatedFeature")
  
  params
}

#modeling bindingFeatures
#Uses the Pareto distribution
# heavy tail: small proportion of values contribute to most,of the distribution's mass
#Mean: exists only if alpha > 1
#arguments:
#pos_start: 
# length: length of binding site
# shape: weight estimation input
# scale: weight estimation input ->
# (shape * scale): mean of the gamma 
# enrichment: used to increase the AvWeight of the binding site --> amplifie binding intensity
# r(): shape of the parameter of pareto dist.

#Using pareto dist (rpareto()), will mean most binding, sites will have small weights, but a few will be extremely high (heavy tail dist)

bindingFeature <- function(start, length=500, shape=1, scale=20, enrichment=5, r=1.5){
  stopifnot(r > 1)
  
  avgWeight <- shape * scale * enrichment
  lowerBound <- ((r - 1) * avgWeight)
  
  weight <- actuar::rpareto1(1, r, lowerBound)
  
  params <- list(start = start, length = length, weight = weight)
  class(params) <- c("Binding", "SimulatedFeature")
  
  params
}

# we now have two functions that creates obect that assigns class to object. S3 is essentially OOP in R
  


# save the new feature functions to a vector which we will implement into the ChIPsim later on

generator <- list(Binding=bindingFeature,Background=backgroundFeature)
  

 
features <- ChIPsim::placeFeatures(generator, 
                                   transition, 
                                   initial, 
                                   start = 0, 
                                   length = Biostrings::width(genome),
                                   globals=list(shape=1, scale=20),
                                   experimentType="TFExperiment",      lastFeat=c(Binding = FALSE, Background = TRUE))

bindIdx <- sapply(features, inherits, "Binding")

 
constRegion <- function(weight, length) rep(weight, length)
featureDensity.Binding <- function(feature, ...) constRegion(feature$weight, feature$length)
featureDensity.Background <- function(feature, ...) constRegion(feature$weight, feature$length)
  

 
featureDensity.Binding <- function(feature, ...){
  featDens <- numeric(feature$length)
  featDens[floor(feature$length/2)] <- feature$weight
  featDens
}
  



dens <- ChIPsim::feat2dens(features, length = Biostrings::width(genome))
  


fragLength <- function(x, minLength, maxLength, meanLength, ...){
  sd <- (maxLength - minLength)/4
  prob <- dnorm(minLength:maxLength, mean = meanLength, sd = sd)
  prob <- prob/sum(prob)
  prob[x - minLength + 1]
}
  


readDens <- ChIPsim::bindDens2readDens(dens, fragLength, bind = 50, minLength = 150, maxLength = 250,
                                       meanLength = 200)
  


featureArgs <- list(generator=generator, transition=transition, init=initial, start = 0, 
                    length = Biostrings::width(genome), globals=list(shape=1, scale=20), experimentType="TFExperiment", 
                    lastFeat=c(Binding = FALSE, Background = TRUE), control=list(Binding=list(length=50)))


readDensArgs <- list(fragment=fragLength, bind = 50, minLength = 150, maxLength = 250,
                     meanLength = 200)
  


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
        #readSeq[idx] <- readError(readSeq[idx], decodeQuality(readQual[idx]))
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
nReads <- 50000

simulated <- simChIP(nReads, genome, file = "test", functions = myFunctions, 
                     control = ChIPsim::defaultControl(features=featureArgs, readDensity=readDensArgs))

#-------------------------------------
#original coding starts here
#unit test for simChIP function
test_that("defaultFunctions returns a valid function list", {
  myFunctions <- defaultFunctions()
  
  # Ensure defaultFunctions() returns a list
  expect_type(myFunctions, "list")
  
  # Ensure readSequence function exists in the list
  expect_true("readSequence" %in% names(myFunctions))
  
  # Replace readSequence function and check if replacement is successful
  myFunctions$readSequence <- dfReads
  expect_identical(myFunctions$readSequence, dfReads)
})
# I only had time to create a unit-test for the simChIP fxn.
##--------------------------------------

#formating reads, read_names and qscores for file output
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
#exporting FASTQ

writeLines(unlist(lines), "../../data/sim_reads.fastq")


#----------------------------------------
# N's insertion :distinguish_simulated_reads function here
#--------------------------------------
distinguish_simulated_reads <- function(reads, reference_genome) {
  results <- lapply(reads, function(read) {
    print(paste("Processing read:", read$id))
    
    seq <- read$seq
    #takes the symbol read qual and makes a number
    quality <- utf8ToInt(read$quality) - 33
    
    #######################################################
    start_pos <- gsub("@read_", "", read$id)
    start_pos <- as.numeric(start_pos)
    #######################################################
    
    #Explicitly prints the sequence
    print(paste("Read sequence:", seq))
    seq <- trimws(toupper(seq))
    
    #Explicitly prints the reference seq found by parsing the file assuming proper formatting (this could be the issue)
    ref_seq <- as.character(subseq(reference_genome[[1]], start = start_pos, width = nchar(seq)))
    print(paste("Ref sequence:", ref_seq)) # Modified to print ref_seq
    ref_seq <- trimws(toupper(ref_seq))
    
    print(paste("subseq start:", start_pos, ", width:", nchar(seq), ", ref_seq length:", nchar(ref_seq)))
    
    #matches should not be where the seq matches the ENTIRE ref seq, but where the bases match. i think this is the core of the problem.
    matches <- seq == ref_seq
    print(paste("Matches:", paste(matches, collapse = ",")))
    mismatch_qualities <- quality[!matches]
    match_qualities <- quality[matches]
    
    print(paste("Mismatch qualities:", paste(mismatch_qualities, collapse = ","))) # Added
    
    # we are getting the mean quality of all the bases that are matches here and all the bases that aren't matches
    # this is only executing the else branch and mismatch_qualities is not the problem. 
    if (length(mismatch_qualities) > 0 && length(match_qualities) > 0) {
      mean_mismatch_quality <- mean(mismatch_qualities)
      mean_match_quality <- mean(match_qualities)
      diff <- mean_match_quality - mean_mismatch_quality
      return(list(read_id = read$id, quality_diff = diff))
    } else {
      print("No mismatches or matches found.")
      return(list(read_id = read$id, quality_diff = NA))
    }
  })
  
  quality_diffs <- sapply(results, function(x) x$quality_diff)
  
  #Determine a threshold.
  #threshold <- quantile(quality_diffs, 0.95, na.rm = TRUE) # Original: Quantile-based threshold
  threshold <- 5  # jank hard-coded threshold
  print(paste("Threshold:", threshold))
  
  is_simulated <- quality_diffs < threshold
  print(paste("Is Simulated:", paste(is_simulated, collapse = ",")))
  
  return(data.frame(read_id = sapply(results, function(x) x$read_id), quality_diff = quality_diffs, is_simulated = is_simulated))
}
#---------------------
results <- distinguish_simulated_reads(reads, genome)
print(results)

#----------------------------------------------------------------------
# End of N's contribution
#---------------------------------------------------------





