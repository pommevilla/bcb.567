'''
    File name: bcb567_project2_classes.py
    Author: Paul Villanueva
    Date created: 11/29/2017
    Date last modified: 11/29/2017
    Python Version: 2.7
'''


import bcb567_project3_utils as bcbutils
        
class MultiSequenceAlignment:

    def __init__(self, dnastrings, word_model, wlcut):
        '''
            dnastrings - a list of tuples (h, s) of header h and dna string s
            word_model - a sequence of 1s and 0s
            wlcut - an integer
        '''
        
        self.dnastrings = dnastrings
        self.num_sequences = len(dnastrings)
        self.superstring = bcbutils.create_superstring(dnastrings)

        self.superword_array = SuperwordArray(self.superstring, word_model, wlcut)
        self.alignment_block_list = self.find_alignment_blocks()
        self.base_pairs_aligned = self.calc_bp_aligned()

    def find_alignment_blocks(self):
        
        candidate_blocks = self.superword_array.find_sized_blocks(self.num_sequences)
        

        alignment_blocks = []
        check_list = [i + 1 for i in range(self.num_sequences)]

        current_list = []
        for i in candidate_blocks:
            indices = self.superword_array.sorted_superwords[i:i + self.num_sequences]

            for j in indices:
                current_list.append(self.id_string(j))
            if current_list == check_list:
                alignment_blocks.append(i)
            current_list = []
                
        return alignment_blocks
            
    def id_string(self, idx):
        """
            idx - the position in the superstring
            
            outputs the dna string that the index belongs to 
        
        >>> l = bcbutils.read_fasta_file_new("./test_files/fasta_files/multi_2.txt")
        >>> msa = MultiSequenceAlignment(l, '1', 1)
        >>> msa.id_string(15)
        1
        
        """
        idx -= 1
        if self.superstring[idx] == '#':
            return 0
        for i, dnastr in enumerate(self.dnastrings):
            if idx - len(dnastr[1]) < 0:
                return i + 1
            else:
                idx -= len(dnastr[1])
                
        return 0
        
    def calc_bp_aligned(self):
        superword_length = self.superword_array.wlcut * self.superword_array.wlength
        bp_total = superword_length
        for i in range(len(self.alignment_block_list) - 1):
            bp_total += superword_length
            overlap = self.alignment_block_list[i + 1] - self.alignment_block_list[i]
            if overlap < superword_length:
                bp_total -= overlap
        
        return bp_total
            
    def rel_str_position(self, idx):
        """
            idx = the absolute position in the super string
            outputs the relative position in the substring that idx belongs to 
        """
        home_str = self.id_string(idx)
        idx -= 1
        for i in range(home_str - 1):
            idx -= len(self.dnastrings[i][1]) + 1
            
        return idx + 1
    
        
    def __str__(self):
        dna_sequences_info = "\n"
        for dnastr in self.dnastrings:
            temp = ">{}\n{}\n".format(*dnastr)
            dna_sequences_info += temp
            
        word_model_info = "Word model: {}".format(self.superword_array.word_model)
        wlcut_info = "wlcut: {}\n".format(self.superword_array.wlcut)
        
        longest_chain_info = "The length of a longest chain of superword blocks: {}".format(self.base_pairs_aligned)
        block_chain_info = "The number of superword blocks in the chain: {}\n".format(len(self.alignment_block_list))
        
        
        
        block_info = ""
        
        superword_length = self.superword_array.wlcut * self.superword_array.wlength
        
        
        for i in range(len(self.alignment_block_list)):
            temp_block_info = "Block {}\n".format(i  + 1)
            
            block_start = self.alignment_block_list[i]

            for j in range(self.num_sequences):
                seq_name = self.dnastrings[j][0]
                rel_pos = self.rel_str_position(self.superword_array.sorted_superwords[block_start + j])
                sw = self.dnastrings[j][1][rel_pos - 1:rel_pos + superword_length - 1]
                temp_block_info += "{:<10} {:<9}{}\n".format(seq_name, rel_pos, sw)
                
            block_info += temp_block_info + "\n"
                
        
        return '\n'.join([dna_sequences_info, word_model_info, wlcut_info, longest_chain_info, block_chain_info, block_info])
        
        
class SuperwordArray:

    def __init__(self, dnastring, word_model, wlcut):
        '''
            dnastring - a string representing a dna sequence
            word_model - a sequence of 1s and 0s
            wlcut - an integer
        '''
        self.dnastring = dnastring
        self.word_code_array = bcbutils.create_word_code_array(dnastring, word_model) 
        self.word_model = word_model
        self.wlength = len(word_model)
        self.wlcut = wlcut
        self.sorted_superwords = bcbutils.create_superword_array(self.word_code_array, self.wlength, self.wlcut)
        
    def get_superword_code(self, pos):
        '''
        inputs:
            pos - the position that we want the superword array from
            
        output:
            returns the superword code of word level wlcut starting at position pos in self.word_code_array
        
        >>> wca = [4, 3, 13, 6, 11, 14, 8, -1, -1, 3, 13, 6, 11, 14, 11, -1, -1]
        >>> self.get_superword_code(wca, 1, 2, 3)
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
        
        
    def find_sized_blocks(self, n):
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
            
        current_block_size = 0
        current_block_start = 0
        current_block_end = 0  
        
        block_list = []
        
        i = 0
        while i < len(self.sorted_superwords) - n + 1:
            
            current_superword = self.get_superword_code(self.sorted_superwords[i] - 1)

            if -1 not in current_superword:
                next_superword = self.get_superword_code(self.sorted_superwords[i + 1] - 1)
                
                current_block_start = i
                current_block_end = i 
                
                while current_superword == next_superword and current_block_end < len(self.sorted_superwords) - 1:
                        
                    current_block_end += 1
                    try:
                        next_superword = self.get_superword_code(self.sorted_superwords[current_block_end + 1] - 1)
                    except IndexError:
                        break
                    
                    
                current_block_size = current_block_end - current_block_start + 1
                
                
                if current_block_size == n:
                    
                    block_list.append(current_block_start)
                    
                i = current_block_end
            i += 1

        return block_list
        
        
    def __str__(self):
    
    
        inputinfo = "{:>40}:    {}\n{:>40}:    {}".format("Word model", self.word_model, "wlcut", self.wlcut)
        
        blockinfo = "{:>40}:    {}".format("Number of positions in largest block", self.max_block_size)
        
        superwordarrayinfo = "{:>40}: {}".format("Positions in superword array", ''.join(['{:>4}'.format(str(i + 2)) for i in range(self.max_block_start - 1, self.max_block_end)]))
        
        positioninfo = "{:>40}: {}".format("Positions in DNA sequence", ''.join(['{:>4}'.format(str(self.sorted_superwords[i])) for i in range(self.max_block_start , self.max_block_end + 1)]))
        
        if int(self.wlcut) * int(self.wlength) == 1:
            superwordlength = 1
        else:   
            superwordlength = int(self.wlcut) * int(self.wlength)
        
        
        
        wordcodeinfo = "{:>40}:    {}".format("Largest block string", ' '.join(i for i in self.dnastring[1][self.sorted_superwords[self.max_block_start] - 1:self.sorted_superwords[self.max_block_start] + superwordlength - 1]))
        
        superwordinfo = "{:>40}:    {}".format("Superword of largest block", ' '.join(str(i) for i in self.max_block_superword))
        
        return '\n'.join([inputinfo, blockinfo, superwordarrayinfo, positioninfo, wordcodeinfo, superwordinfo])
        

        
if __name__ == "__main__":
    pass
    
    # l = bcbutils.read_fasta_file("./test_files/fasta_files/multi_2.txt")
    # wm = bcbutils.read_word_model("./test_files/word_models/wm_101.txt")
    # wlcut = 2
    # msa = MultiSequenceAlignment(l, wm, wlcut)
    
    # print msa

    

 