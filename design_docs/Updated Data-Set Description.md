Real Data-Set Description (4/11/25)
Flores JC, Sidoli S, Dawlaty MM. Tet2 regulates Sin3a recruitment at active enhancers in embryonic stem cells. iScience 2023 Jul 21;26(7):107170. PMID: 37456851
https://pmc.ncbi.nlm.nih.gov/articles/PMC10338317/

Short-Description:
These dataset were taken from WT and KO loss of function CUT&TAG and RNA-Seq experiments of Tet2's influence on the epigenomic landscape and the of mmu-ESCs. The Tet family of enzymes is important TF for mammalian development and are known for iterative oxidation of 5mC and other 5mC variations. Each of these histone alterations signals other regulatory factors and promote demthylation thus opening up chromatin for expression.
Tet2 has distinct catalytic and noncatalytic target genes in ESCs based on the histone landscape (hetero/euchromatin). In the experiment the Dawlaty et. al. use the high throughput loss of function experiments listed earlier to find additional signs of regulatory feedback due to Tet2 inactivity.



RNA-Seq Datasets Total: 8Gb
Tet2_wildtype_RNAseq_rep1, GSM6585594, SAMN30862176 (Paired-End)
wget s3://sra-pub-src-13/SRR21586232/Tet2_wildtype_RNAseq_rep1_2.fq.gz.1
wget s3://sra-pub-src-13/SRR21586232/Tet2_wildtype_RNAseq_rep1_1.fq.gz.1


Tet2_knockout_RNAseq_rep1, GSM6585599, SAMN30862171(Paired-End)
wget s3://sra-pub-src-12/SRR21586226/Tet2_knockout_RNAseq_rep1_2.fq.gz.1
wget s3://sra-pub-src-16/SRR21586226/Tet2_knockout_RNAseq_rep1_1.fq.gz.1	


ChIP-Seq Datasets Total: 6.6Gb

Tet2_wildtype_H3K27ac_CUTnTag,SAMN30862166, SRR21586222 (Paired-End)
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR21586222/SRR21586222
~3.3Gb


Tet2_knockout_H3K27ac_CUTnTag, SAMN30862164, SRR21586220, (Paired-End)
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR21586220/SRR21586220
~3.3Gb


use sra-toolskits to unpack the fastqc files
