import doctest

def read_fasta(fasta_file):
    """
    fasta_file - a .txt file in FASTA format
    
    returns a string representing a DNA sequence
    
    >>> read_fasta("A.short.txt")
    ('Ashort', 'GATCGTAGAGTGAGACCTAGTGTTTG')
    
    >>> read_fasta("B.short.txt")
    ('Bshort', 'CTCGTAGGTGAGATTCCTAGTGCC')
    """
    
    with open(fasta_file) as fin:
        header = fin.readline()[1:].strip()
        fasta_string = ""
        for line in fin:
            fasta_string += line.strip()
            
    return header, fasta_string
    
def populate_matrices(A, B, MATCH_SCORE, mismatch_score, gap_open, gap_extend):
    pass
    
def traceback(S, D, I, A, B, gap_init_penalty):
    pass

class DNAString:

    def __init__(self, header, seq):
        self.header, self.seq = header, seq
       
       
    @property
    def length(self):
        """
        >>> string1 = read_fasta("A.short.txt")
        >>> A = DNAString(*string1)
        >>> A.length
        26
        """
        return len(self.seq)
        
    # @property
    # def header(self):
        # """
        # >>> string1 = read_fasta("A.short.txt")
        # >>> A = DNAString(*string1)
        # >>> A.header
        # Ashort
        # """
        # return self.header
        
    def __iter__(self):
        self.current = 0
        return self
        
    def __next__(self):
        if self.current == self.length - 1:
            raise StopIteration
        else:
            self.current = self.current + 1
            return self.seq[self.current]
    
    def __str__(self):
        """
        >>> string1 = read_fasta("A.short.txt")
        >>> A = DNAString(*string1)
        >>> print(A)
        Ashort
        GATCGTAGAGTGAGACCTAGTGTTTG
        """
        return "{}\n{}".format(self.header, self.seq)

class LocalAlignment:
    def __init__(self, dnastring1, dnastring2, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty):
        self.dnastring1 = dnastring1 
        self.dnastring2 = dnastring2
        self.match_score = match_score 
        self.mismatch_score = mismatch_score 
        self.gap_open_penalty = gap_open_penalty
        self.gap_extend_penalty = gap_extend_penalty
        self.S_matrix, self.D_matrix, self.I_matrix, self.max_alignment_score, self.max_score_row, self.max_score_column = populate_matrices(self.dnastring1, self.dnastring2, self.match_score, mismatch_score, gap_open_penalty, gap_extend_penalty)
        self.optimal_alignment = traceback(self.S_matrix, self.D_matrix, self.S_matrix, self.dnastring1, self.dnastring2, self.gap_init_penalty)
        
    @property
    def length(self):
        """
            HUH
        """
        return len(self.optimal_alignment)
        
    def __str__(self):
        """
        >>> string1 = read_fasta("A.short.txt")
        >>> A = DNAString(*string1)
        >>> string2 = read_fasta("B.short.txt")
        >>> B = DNAString(*string2)
        >>> AB_local_alignment = LocalAlignment(A, B, 10, -20, 40, 2)
        >>> print(AB_local_alignment)
        Something
        """
        score_info = "Math Score \t Mismatch Score \t Gap Open Penalty \t Gap Extension Penalty\n{}\t{}\t{}\t{}\n".format(
            self.match_score, self.mismatch_score, self.gap_open_penalty, self.gap_extend_penalty)
            
        seq_info = "Sequence A: {}\n Length: {}\nSequence B: {}\n Length: {}\n".format(
            dnastring1.header, dnastring1.length, dnastring2.header, dnastring2.length)
            
        align_info = "Alignment Score: {}\n Length: {}\n Start in A: {}\nStart in B: {}\nEnd in A: {}\nEnd in B: {}".format(
            self.max_alignment_score, self.length, self.max_score_row, self.max_score_column, 22, 22)
            
        aligned_string_A = ''
        aligned_string_mid = ''
        aligned_string_B = ''
        
        for i in self.optimal_alignment:
            aligned_string_A += i[0]
            aligned_string_B += i[1]
            if i[0] == '-' or i[1] == '-':
                aligned_string_mid += "-"
            else:
                aligned_string_mid += "|"
    
        aligned_strings = "\t{}\n\t{}\n\t{}".format(aligned_string_A, aligned_string_mid, aligned_string_B)
        
        return "".join(score_info, seq_info, align_info, aligned_strings)
        
if __name__ == "__main__":
    doctest.testmod()
    string1 = read_fasta("A.short.txt")
    print(*string1)
    A = DNAString(*string1)
    print(A)
    A.header = "Something"
    print(A)
    print(A.header)

    
        