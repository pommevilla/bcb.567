import doctest
import numpy as np

def read_fasta(seq):
    """
    seq - a .txt file in FASTA format
    
    returns a string representing a DNA sequence
    
    >>> read_fasta("A.short.txt")
    'GATCGTAGAGTGAGACCTAGTGTTTG'
    
    >>> read_fasta("B.short.txt")
    'CTCGTAGGTGAGATTCCTAGTGCC'
    """
    
    with open(seq) as fin:
        header = fin.readline()
        fasta_string = ""
        for line in fin:
            fasta_string += line.strip()
            
    return fasta_string
    
def delta(a, b, match, mismatch):
    """
    a, b - two characters 
    match - a positive integer score for a character match
    mismatch - a negative integer score for a character mismatch
    
    >>> delta("C", "C", 10, -20)
    10
    
    >>> delta("C", "T", 10, -30)
    -30
    
    """
    return match if a is b else mismatch
        
    
def populate_matrices(A, B, MATCH_SCORE, mismatch_score, gap_open, gap_extend):
    """
    A, B - two DNA sequences 
    mismatch - a negative integer indicating mismatch scare
	gap_open - a non-negative integer indicating gap_open penalty
	gap_extend - a positive integer indicating gap extension penalty
    
    returns matrices of integers S, D, and I, where:
        -S[i][j] is the maximum score of all local alignments starting at A[i + 1] and B[j + 1]
        -D[i][j] is the maximum score of those local alignments that begin with a deletion gap
        -I[i][j] is the maximum score of those local alignments that begin with an insertion gap
        
    >>> populate_matrices('GATCGTAGAGTGAGACCTAGTGTTTG', 'CTCGTAGGTGAGATTCCTAGTGCC', 10, -20, 40, 2)
    'SO FAR SO GOOD'
        
    """
    m = len(A)
    n = len(B)
    S, D, I = [np.zeros((m + 1, n + 1)) for _ in range(3)]
    print("S is a {} matrix".format(S.shape))
    
    row_first = m
    column_first = n
    
    gap_init_penalty = gap_open + gap_extend
    
    score = 0
    
    #S[m][n] = 0
    D[m][n] = -(gap_open + gap_extend)
    I[m][n] = -(gap_open + gap_extend)
    
    for j in reversed(range(n - 1)):
        #S[m][j] = 0
        D[m][j] = -(gap_open + gap_extend)
        I[m][j] = -(gap_open + gap_extend)

        
    for i in reversed(range(m - 1)):
        #S[i][n] = 0
        D[i][n] = -(gap_open + gap_extend)
        I[i][n] = -(gap_open + gap_extend)
        
        for j in reversed(range(n - 1)):
            D[i][j] = max(D[i + 1][j] - gap_extend,
                            S[i + 1][j] - gap_open - gap_extend)
            I[i][j] = max(I[i][j + 1] - gap_extend, 
                            S[i][j + 1] - gap_open - gap_extend)
            S[i][j] = max(0, S[i + 1][j + 1] + delta(A[i + 1], B[j + 1], MATCH_SCORE, mismatch_score),
                            D[i][j], I[i][j])
            if score < S[i][j]:
                score = S[i][j]
                row_first = i
                column_first = j 
                
    print(S)
        
    
    
    
    
def traceback(S, D, I):
    """
    S, D, I - the three matrices returned by populate_matrices(A, B)
    
    returns integers i, j, m, n corresponding to the intervals A[i:j] and B[m:n] that give optimal local alignment
    """
    print("Not yet implemented.")

    
    
"""
def local_alignment(seq1, seq2, mismatch_score, gap_open, gap_extend):
	seq1, seq2 - .txt files in FASTA format
	mismatch - a negative integer value for mismatch score
	gap_open - a non-negative integer value for gap_open penalty
	gap_extend - a positive integer value for gap extension penalty
    
    
    MATCH_SCORE = 10
    # read in FASTA strings, return DNA sequences
    dna_string1 = read_fasta(seq1)
    dna_string2 = read_fasta(seq2)
    
    # compute dynamic programming matrices
    S, D, I = populatate_matrices(dna_string1, dna_string2, MATCH_SCORE, mismatch_score, gap_open, gap_extend)
    
    # trace back optimal alignment through matrices
    opt_align = traceback(S, D, I)
    
    # print out alignment
    print("Not yet implemented.")
"""

if __name__ == '__main__':
    #doctest.testmod()
    populate_matrices('GATCGTAGAGTGAGACCTAGTGTTTG', 'CTCGTAGGTGAGATTCCTAGTGCC', 10, -20, 40, 2)