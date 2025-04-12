#!/bin/bash
#Script to attain the simulated data for example uses
cd ./sim_data/inputs

if ! command -v samtools &> /dev/null; then
    echo "Error: samtools is not installed or not in PATH."
    exit 1
fi

echo "samtools found. Continuing..."


#Setting references

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr1.fa.gz
gunzip chr1.fa.gz

wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/genes/mm10.ncbiRefSeq.gtf.gz
gunzip mm10.ncbiRefSeq.gtf.gz
awk '$1 == "chr1"' mm10.ncbiRefSeq.gtf > chr1_mm10_ncbirefseq.gtf

rm mm10.ncbiRefSeq.gtf


#trimming the gtf file
awk '$1 == "chr1" && $3 == "exon" && $4 >= 1 && $5 <= 10000000' chr1_mm10_ncbirefseq.gtf > chr1_mm10_1M.gtf

rm chr1_mm10_ncbirefSeq.gtf


#need to install samtools
samtools faidx chr1.fa chr1:1-10000000 > chr1_10M.fa

rm chr1.fa
rm chr1.fa.fai















