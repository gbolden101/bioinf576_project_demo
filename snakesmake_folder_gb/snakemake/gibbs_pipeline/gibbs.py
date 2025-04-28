import numpy as np
import seqlogo

# Import functions from previous class for building sequence motif & identifying seqs matching the motif
from data_readers import *
from seq_ops import get_seq
from motif_ops import *

# Some of the code in this notebook uses other programs installed in the Python
# environment.  We need to add the b529 Conda env to our $PATH environment
# variable so that code knows where to find those programs.
import sys, pathlib, os


def GibbsMotifFinder(seqs, k):
    '''
    Function to find a pfm from a list of strings using a Gibbs sampler
    
    Args: 
        seqs (str list): a list of sequences, not necessarily in same lengths
        k (int): the length of motif to find

    Returns:
        pfm (numpy array): dimensions are 4xlength
        
    '''

    N = len(seqs) #number of sequences input
    
    # First we randomly pick k-length sequences (kmer) from each sequence and place in list motifs to start
    motifs = [] 
    for i in range(N):
        if len(seqs[i]) >= k:
            # Generate random k-mer and strand
            index = np.random.randint(len(seqs[i])-k+1)
            strand = np.random.randint(2)
            if strand == 0:
                motifs.append(seqs[i][index:index+k])
            else:
                motifs.append(reverse_complement(seqs[i][index:index+k]))
    
    # Next 20000 iterations or to convergence
    # For my implementation, I define convergence as when the information content of the motif
    # does not change between two measurements that are 100 iterations apart.
    # The logic behind this is that, while there are some fluxuations in the random walk
    # there is a very low chance that nothing will have changed in 100 iterations unless
    # convergence has occured.
    last_ic = None # Store IC for testing convergence later

    for j in range(20000):
        i = np.random.randint(N) # Randomly pick a sequence to be excluded when building motif
                
        # Check for convergence every 100 iterations
        if (j % 100 == 0):
            pfm = build_pfm(motifs,k)  # Build PFM from all motifs
            print(pfm_ic(pfm))         # Output IC
            if last_ic == pfm_ic(pfm): # End if IC from 100 iterations ago is the same as now
                break
            last_ic = pfm_ic(pfm)
        
        # Build a PFM and PWM without the ith sequence that was randomly selected
        pfm = build_pfm(motifs[:i]+motifs[i+1:],k)
        pwm = build_pwm(pfm)

        # Score all kmers in the sequence using our score_kmer() function from last class
        # We do this forward and then backwards and just append the scores and strand to a list.
        # This is similar to the score_nmer() function from last class, but instead of
        # identifying the top score, we just keep track of all of the scores.
        # Using the function score_nmer() from last class would be akin to the hill climbing version
        # of this algorithm.
        scores = []
        strands = []
        for pos in range(len(seqs[i])-k+1):                 # For all of the positions
            score_temp = score_kmer(seqs[i][pos:pos+k],pwm) # Score for the forward strand
            scores.append(score_temp)
            strands.append(0)                               # Strand 0 to represent forward
        
        #Also reverse strand
        for pos in range(len(seqs[i])-k+1):
            score_temp = score_kmer(reverse_complement(seqs[i][pos:pos+k]), pwm) # Score for the negative strand
            scores.append(score_temp)
            strands.append(1)                                                    # Strand 1 to represent negative

        # Now convert scores to 'probabilities' as defined above
        scores = np.exp2(scores) # 2^score to move out of log2 space
        scores_norm = [score/np.sum(scores) for score in scores] # Score / sum scores to generate 'probabilities'
        
        # Choose a new motif location in a probabilistic fashion
        # To do this we draw an index from distribution of motif scores
        # We can just pass our 'probabilities' above to the np.random.choice() function along with the matching
        # indexes for each position
        index_new = np.ndarray.item(np.random.choice((len(seqs[i])-k+1)*2,1,p=scores_norm))
        strand_new = strands[index_new]
        
        index_new = index_new % (len(seqs[i])-k+1) # Because our index list is twice the length (forward and backward)
                                                   # the new index is actually modulo the length of the number of kmers
                                                   # that were tested in the sequence above.

        # Finally, replace the motif for ith sequence with the probabilistically selected kmer
        if strand_new == 0:
            motifs[i] = seqs[i][index_new:index_new+k]
        else:
            motifs[i] = reverse_complement(seqs[i][index_new:index_new+k])
    
    return pfm


fasta_file = "../data/gibbs_data/tet2_seq.fa"
gff_file = "../data/gibbs_data/extracted.gtf"

seqs = []

for name, seq in get_fasta(fasta_file):
    chrom = name.split(':')[0]
    offset = int(name.split(':')[1].split('-')[0])  # 36770235
    
    for gff_entry in get_gff(gff_file):
        if gff_entry.type == 'CDS' and gff_entry.seqid == chrom:
            # Convert GTF genomic positions to relative coordinates
            relative_start = gff_entry.start - offset
            relative_end = gff_entry.end - offset

            # Skip invalid coordinates
            if relative_start < 0 or relative_end > len(seq):
                continue

            promoter_seq = get_seq(seq, relative_start, relative_end, gff_entry.strand, 50)

            # Keep only sequences long enough for motif finding
            if promoter_seq and len(promoter_seq) >= 10:
                seqs.append(promoter_seq)


# Run the Gibbs sampler
seqs = [s for s in seqs if len(s) >= 10]
promoter_pfm = GibbsMotifFinder(seqs, 10)
# Plot the final PFM that is generated
seqlogo.seqlogo(
    seqlogo.CompletePm(pfm=promoter_pfm.T),
    format='pdf',
    filename='../gibbs_output/motif_logo.pdf'
)
