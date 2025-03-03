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

#Making the referance for polyester
#download GRCm38.p6.genome.fa and mm10.ncbiRefSeq.gtf from NCBI
  # chr1_ann.fa made using the following commands cmdline
  # extract chr1 from the the mm10 gtf file: mm10.ncbiRefSeq.gtf
  # awk '$1 == "chr1"' mm10.ncbiRefSeq.gtf > chr1_mm10_ncbirefseq.gtf
  # samtools faidx GRCm38.p6.genome.fa chr1 > chr1.fa
  # bedtools getfasta -name -fi chr1.fa -bed chr1_mm10_ncbirefseq.gtf > chr1_ann.fa

#Making the referance genome for ChIPSim
# This referance genome should be just the sequences made from chr1_ann.fa
#download GRCm38.p6.genome.fa and mm10.ncbiRefSeq.gtf from NCBI
# awk '$1 == "chr1"' mm10.ncbiRefSeq.gtf > chr1_mm10_ncbirefseq.gtf
# samtools faidx GRCm38.p6.genome.fa chr1 > chr1.fa
# start=$(awk '$1 !~ /^#/ {print $4}' chr1_mm10_ncbirefseq.gtf | sort -n | head -1)
# end=$(awk '$1 !~ /^#/ {print $5}' chr1_mm10_ncbirefseq.gtf | sort -n | tail -1)
# samtools faidx chr1.fa chr1:$start-$end > chr1_limit.fa
# you can use this file as the referance genome for ChIP script



Next Goals: 
1. Start constructing the snamke_make pipeline to access out sim_data
