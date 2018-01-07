'''
    File name: bcb567_project2_classes.py
    Author: Paul Villanueva
    Date created: 10/7/2017
    Date last modified: 10/26/2017
    Python Version: 2.7
'''

import doctest
import bcb567_project2_utils as bcbutils

class SuperwordArray:

    def __init__(self, dnastring, word_model, wlcut):
        '''
            dnastring - a tuple (h, s) containing a header h and a dna string s
            word_model - a sequence of 1s and 0s
            wlcut - an integer
        '''
        self.dnastring = dnastring
        self.word_code_array = bcbutils.create_word_code_array(dnastring[1], word_model) 
        self.word_model = word_model
        self.wlength = len(word_model)
        self.wlcut = wlcut
        self.sorted_superwords = bcbutils.create_superword_array(self.word_code_array, self.wlength, self.wlcut)
        self.max_block_start, self.max_block_end, self.max_block_size, self.max_block_superword = self.find_largest_block()
        
    def get_superword_code(self, pos):
        '''
        inputs:
            pos - the position that we want the superword array from
            
        output:
            returns the superword code of word level wlcut starting at position pos in self.word_code_array
        
        >>> wca = [4, 3, 13, 6, 11, 14, 8, -1, -1, 3, 13, 6, 11, 14, 11, -1, -1]
        >>> get_superword_code(wca, 1, 2, 3)
        (3, 6, 14)
        
        >>> get_superword_code(wca, 5, 3, 2)
        (14, -1)
        
        >>> get_superword_code(wca, 7, 2, 4)
        (-1, 3, 6, 14)
        
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
        '''
        A block in superword array is a subsequence of self.sorted_superwords where each index corresponds
        to the same superword.
        
        find_largest_block returns information about the first largest block found in self.sorted_superwords.
        
        output:
            Let M be the largest block of self.sorted_superwords
            max_block_start - the beginning in self.sorted_superwords of M
            max_block_end - the end in self.sorted_superwords of M
            max_block_size - the number of elements in M
            max_block_superword - the superword code associated with M
            
        '''
            
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
        inputinfo = "{:>40}: {:>4}\n{:>40}: {:>4}".format("Word model", self.word_model, "wlcut", self.wlcut)
        
        blockinfo = "{:>40}: {:>4}".format("Number of positions in largest block", self.max_block_size)
        
        superwordarrayinfo = "{:>40}: {}".format("Positions in superword array", ''.join(['{:>4}'.format(str(i + 1)) for i in range(self.max_block_start - 1, self.max_block_end)]))
        
        positioninfo = "{:>40}: {}".format("Positions in largest block", ''.join(['{:>4}'.format(str(self.sorted_superwords[i])) for i in range(self.max_block_start - 1, self.max_block_end)]))
        
        superwordinfo = "{:>40}:   {}".format("Superword of largest block", self.max_block_superword)
        
        return '\n'.join([inputinfo, blockinfo, superwordarrayinfo, positioninfo, superwordinfo])
 