DDS Chart Summaries

ChIP seq Pre-Processing Pipeline Overview
FASTQ --(FastQC)--> FASTQ --(Trimmomatic/Cutadapt)--> 
Trimmed FASTQ --(Bowtie2)--> SAM --(SAMtools view)--> 
BAM --(SAMtools sort)--> Sorted BAM --(SAMtools rmdup)--> 
Deduplicated BAM --(MACS2)--> Peak Files (narrowPeak/broadPeak) --(HOMER/ChIPseeker)-->
Annotated Peaks (CSV/TSV)


RNA seq Pre-Processing Pipeline Overview
FASTQ --(FastQC)--> FASTQ --(Trimmomatic/Cutadapt)--> 
Trimmed FASTQ --(FastQC)--> Trimmed FASTQ --(STAR/HISAT2)--> 
SAM --(SAMtools view)--> BAM --(SAMtools sort)--> Sorted BAM --(SAMtools index)--> 
Indexed BAM --(HTSeq-count/featureCounts)--> Raw Count Matrix (TXT/CSV) --(DESeq2/edgeR)--> 
Differential Expression Results (CSV/TSV)


Rstudio Pipeline Analysis Pipeline Overview

Differential Expression Results + Annotated Peaks --(ClusterProfiler, DAVIDenrichr())--> GSEA results --(GSEA library package)--> Relavant GSEA results --(Scatter Plot + Volcano fxns)--> GSEA Visuals
			|
			|
			|(Motiff Analysis)
			|
			v
			Use Granges, MMU Gen Ref, GTF, for motif finding based on ChIP experiment
			|
			|
			|
			|
			|
			V
		Final Table: TF, Nucleotide Motiff


DDS Chart Descriptions**
Pre-Amble:
During actual development I will be using small toy examples so that I do not need to depend on long wait times
during trial and error phase. Each of the Preprocessing steps will be deposited into a separate directory.
---> Phase1/Module 1: Preprocessing

bulk-RNA-Seq Pre-processing

FastQC assesses the quality of raw sequencing reads (e.g., per base sequence quality, GC content). 
It doesn't modify the FASTQ file but generates quality reports.

We will be trimming low quality bases and adapters using Trimmomatic or Cutadapt. The result will be paired Fwd an paired Rev file.This will improve mapping efficiency in the next steps.

Alignment to the Reference genome using the Trimmed FASTQ file using STAR to create SAM file. STAR are splice-aware aligners optimized for RNA-Seq data
Converting from SAM to BAM using SAMtools (view, sort, index) to get an indexed BAM file. They are compressed an easier to manage
Sorting and indexing will be essential for downstream visualization

featureCounts or HTSeq will convert the indexed BAM into the counts matrix for DESeq2 analysis in RStudio. The ChIP-seq FASTA will have gone through almost the exact same processing up until we get to SAMtool where we will be using SamTools to remove duplicated peaks. MACS2 will then be used to call the peaks. HOMER is then used annotate the Peaks.

The final step for RNA count Matrix and the ChIP-Seq matrix is to go through a statistical analysis to ensure that we truly have statistically significant peaks and DEGs.

Thresholds are not determined at this time. The best statistical validation for the data sets is the correct number of gene samples and 

---> Phase2/Module 2: RStudio Analysis + GSEA 

Once we have determined that there a statistically significant DEGs and Peaks we can combine them into one data table and use clusterProfilr to determine the enriched pathways.

From here we must determine what genes appear in all three categories of gene oncology Cell Component, Biological Process, and Molecular function and are significantly differentially expressed.



---> Phase3/Module 3: Motiff Calling using RS functions

As of now we do not have much research on the Motif calling. I am aware that there are motif functions that directly retrieve motifs using the reference gene.

Tasks are addressed in the WBS


**ATAC-Seq data Inclusion is uncertain at this point. Something must be built first before adding
