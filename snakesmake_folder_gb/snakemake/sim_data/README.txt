Input of sim_chip_seq : Nothing

Output: sim_reads.txt (FASTQ file)
Should be less than 3Kb

Purpose: for the development of a ChIP-Seq preprocessing pipeline

Intrusctions
1. Create a directory with the Rscript inside
2. Navigate to the directory
3. run the Rscript in the terminal using Rscript sim_chip_seq.R

You could either use the inputs provided on the polyester GitHub for sim_rnaseq reads or you could make your own


#How to make example inputs for the function

#Making the referance for polyester and CHIPsim
# working directory at this point: sim_data/inputs
#download GRCm38.p6.genome.fa and mm10.ncbiRefSeq.gtf from NCBI
  #wget ftp://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/Mus_musculus.GRCm38.110.gtf.gz
  #wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz
  # chr1_ann.fa made using the following commands cmdline
  # extract chr1 from the the mm10 gtf file: mm10.ncbiRefSeq.gtf
  # awk '$1 == "chr1"' mm10.ncbiRefSeq.gtf > chr1_mm10_ncbirefseq.gtf
  # samtools faidx GRCm38.p6.genome.fa chr1 > chr1.fa
  # bedtools getfasta -name -fi chr1.fa -bed chr1_mm10_ncbirefseq.gtf > chr1_ann.fa
# start=$(awk '$1 !~ /^#/ {print $4}' chr1_mm10_ncbirefseq.gtf | sort -n | head -1)
# end=$(awk '$1 !~ /^#/ {print $5}' chr1_mm10_ncbirefseq.gtf | sort -n | tail -1)
# end will give a result of 195mil, just use 6million for the end point. Loading times for ChIPSim will be long
# samtools faidx chr1.fa chr1:$start-$end > chr1_limit.fa
# you can use this file as the referance genome for ChIP script

#Running RSscripts
#Work directory should be R
#Rscript "sim_chip file"
#Rscript "sim_rna file"

Next Goals: 
1. Start constructing the snamke_make pipeline to access out sim_data
2. 3/8/25 goals convert sam to bam.
3. implement real vs sim conditional on bowtie 
4.Ultimate Goal: finish ChIP-PP before spring break is over and get the RNA-Seq one done (YOU GOT THIS))
