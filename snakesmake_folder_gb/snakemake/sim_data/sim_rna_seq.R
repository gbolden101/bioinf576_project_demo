if(!require(polyester)) devtools::install_github("alyssafrazee/polyester-release")
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("Biostrings")


library(polyester)
library(Biostrings)

#Polyester github: https://github.com/alyssafrazee/polyester github install: 
# devtools::install_github("alyssafrazee/polyester-release")


#args
#file name of ANNOTATED fasta file
#name of the out directory of where you want the reads to go
# this will generate reads from 20 transcripts
rnaseq_sim_reads <- function(input_fasta_file, read_outdir) {
  
  
  
  #The Github says that the simulated_reads function takes an FC Matrix
  #However, simulated_experiment only takes a vector
  
  fold_changes <- matrix(c(rep(c(1,4,4,1,1),4),rep(c(4,1,1,4,1),4)), nrow = 20)
  #Inputs need for the simulate experiment function
  
  #How to make example inputs for the function
  # chr1_ann.fa made using the following commands cmdline
  # extract chr1 from the the mm10 gtf file: mm10.ncbiRefSeq.gtf
  # awk '$1 == "chr1"' mm10.ncbiRefSeq.gtf > chr1_mm10_ncbirefseq.gtf 
  # samtools faidx GRCm38.p6.genome.fa chr1 > chr1.fa
  # bedtools getfasta -name -fi chr1.fa -bed chr1_mm10_ncbirefseq.gtf > chr1_ann.fa
  
  
  
  fasta_file <- as.character(input_fasta_file)
  
  fasta = readDNAStringSet(fasta_file)
  
  # This lets us cut down the number of transcripts we want
  # I could ask the user but that requires conditionals just in case
  small_fasta = fasta[1:20]
  
  #the small fasta will be the 20 transcripts of annotated mm10_chr1
  #write out the small_fasta_file
  writeXStringSet(small_fasta, './inputs/small_fasta.fa')
  
  #~ 20X coverage reads per transcript = transcriptlength/readlength * 20
  # here all transcripts will have ~equal FPKM
  # This will give us equal read lengths and paired-end
  readspertx <- round(20 * width(small_fasta) / 100)
  
  #call simulated_experiment
  #file, readpertx, num reps, 2 group with 10 replicates each, outpath
  
  # fasta: small fasta
  # num_reps --> number of replicates per group: in this case here we have 2 groups with 10 replicates each, but you can have as many groups as you like)
  
  #outdirectory: user provided
  return(simulate_experiment(fasta = './inputs/small_fasta.fa',
                             reads_per_transcript = readspertx, 
                             num_reps=c(10,10), 
                             fold_changes = fold_changes[,1],
                             outdir=read_outdir))
  
}

input_fasta_file <- './inputs/chr1_ann.fa'
read_outdir <- 'sim_rnaseq_reads'


rnaseq_sim_reads(input_fasta_file, read_outdir)
