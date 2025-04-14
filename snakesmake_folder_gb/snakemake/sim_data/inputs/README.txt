#How to make example inputs for the function
#download GRCm38.p6.genome.fa and mm10.ncbiRefSeq.gtf from NCBI
  # chr1_ann.fa made using the following commands cmdline
  # extract chr1 from the the mm10 gtf file: mm10.ncbiRefSeq.gtf
  # awk '$1 == "chr1"' mm10.ncbiRefSeq.gtf > chr1_mm10_ncbirefseq.gtf
  # samtools faidx GRCm38.p6.genome.fa chr1 > chr1.fa
  # bedtools getfasta -name -fi chr1.fa -bed chr1_mm10_ncbirefseq.gtf > chr1_ann.fa