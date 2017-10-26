'''
    File name: bcb567_project2_classes.py
    Author: Paul Villanueva
    Date created: 10/7/2017
    Date last modified: 10/26/2017
    Python Version: 2.7
'''

import doctest
import bcb567_project2_utils as bcbutils

class DNAString:
    """
    Represents a DNA sequence with a header.  
    """

    def __init__(self, header, seq):
        self.header = header
        self.seq = seq
       
    @property
    def length(self):
        """
        >>> string1 = bcbutils.read_fasta_file("A.short.txt")
        >>> A = DNAString(*string1)
        >>> A.length
        26
        """
        return len(self.seq)
        
    def __len__(self):
        return len(self.seq)

    def __iter__(self):
        self.current = 0
        return self
        
    def __next__(self):
        if self.current == self.length - 1:
            raise StopIteration
        else:
            self.current = self.current + 1
            return self.seq[self.current]
            
    def __getitem__(self, i):
        """
        >>> string1 = bcbutils.read_fasta_file("A.short.txt")
        >>> A = DNAString(*string1)
        >>> print(A[0])
        G
        """
        return self.seq[i]
    
    def __str__(self):
        """
        >>> string1 = bcbutils.read_fasta_file("A.short.txt")
        >>> A = DNAString(*string1)
        >>> print(A)
        Ashort
        GATCGTAGAGTGAGACCTAGTGTTTG
        """
        return "{}\n{}".format(self.header, self.seq)

class SuperwordArray:

    def __init__(self, dnastring, word_model, wlcut):
        self.dnastring = dnastring
        self.word_code_array = bcbutils.create_word_code_array(dnastring[1], word_model) 
        self.word_model = word_model
        self.wlength = len(word_model)
        self.wlcut = wlcut
        self.sorted_superwords = bcbutils.create_superword_array(self.word_code_array, self.wlength, self.wlcut)
        self.max_block_start, self.max_block_end, self.max_block_size, self.max_block_superword = self.find_largest_block()
        
        
        

        
    def get_superword_code(self, pos):
        '''
        >>> wca = [4, 3, 13, 6, 11, 14, 8, -1, -1, 3, 13, 6, 11, 14, 11, -1, -1]
        >>> get_superword_code(wca, 1, 2, 3)
        (3, 6, 14)
        
        >>> get_superword_code(wca, 5, 3, 2)
        (14, -1)
        
        >>> get_superword_code(wca, 7, 2, 4)
        (-1, 3, 6, 14)
        
        >>> get_
        '''
        

        
        wlevel = 0
        n = len(self.word_code_array)
        superword_code = []
        
        for wlevel in range(self.wlcut):
            if pos + self.wlength * wlevel < n:
                superword_code.append(self.word_code_array[pos + self.wlength * wlevel])
            else:
                superword_code.append(-1)
                
        return tuple(superword_code)
        
    def find_largest_block(self):
            
        max_block_size = 0
        max_block_start = 0
        max_block_end = 0
        max_block_superword = ()
        
        current_block_size = 0
        current_block_start = 0
        current_block_end = 0  
        
        i = 0
        while i < len(self.word_code_array):
            sw_i = self.sorted_superwords[i] - 1
            current_superword = self.get_superword_code(sw_i)
            if -1 not in current_superword:
                next_superword = self.get_superword_code(self.sorted_superwords[i + 1])
                current_block_start = i 
                current_block_end = i 
                while current_superword == next_superword:
                        
                    current_block_end += 1
                    sw_i = self.sorted_superwords[current_block_end]
                    next_superword = self.get_superword_code(current_block_end + 1)
                    
                    
                current_block_size = current_block_end - current_block_start + 1
                
                if current_block_size > max_block_size:
                    max_block_start = current_block_start
                    max_block_end = current_block_end
                    max_block_size = current_block_size
                    max_block_superword = current_superword
                    
                i = current_block_end
            i += 1
        
        return max_block_start, max_block_end, max_block_size, max_block_superword
        
        
    def __str__(self):
        inputinfo = "{:>40}: {}\n{:>40}: {}".format("Word model", self.word_model, "wlcut", self.wlcut)
        
        blockinfo = "{:>40}: {}".format("Number of positions in largest block", self.max_block_size)
        
        positioninfo = "{:>40}: {}".format("Positions in largest block", ' '.join([str(self.sorted_superwords[i]) for i in range(self.max_block_start - 1, self.max_block_end)]))
        superwordinfo = "{:>40}: {}".format("Superword of largest block", self.max_block_superword)
        
        return '\n'.join([inputinfo, blockinfo, positioninfo, superwordinfo])
 