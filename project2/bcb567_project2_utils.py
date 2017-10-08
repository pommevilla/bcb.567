'''
	File name: bcb567_project1_utils.py
	Author: Paul Villanueva
	Date created: 9/13/2017
	Date last modified: 9/20/2017
	Python Version: 3.6
'''

import doctest

def sectioner(seq, n = 70):
	"""
	seq - any sequence- or list-like object
	n - desired size of chunks.  Defaulted to 70 per assignment directions.
	
	returns a list of subsequences of seq, each of size n.  
	"""
	return (seq[i:i + n] for i in range(0, len(seq), n))

def read_word_model(word_model_file):
    """
    word_model_file - a .txt file containing a string of 1s and 0s on a single line 
    
    returns a string of 1s and 0s
    >>> read_word_model('word_model_1.txt')
    '1101010100011111111010111010101'
    
    >>> read_word_model('word_model_2.txt')
    '11111101011110001110'
    """
    with open(word_model_file) as fin:
        return fin.readline()


def read_fasta_file(fasta_file):
	"""
	fasta_file - a .txt file in FASTA format
	
	returns a string representing a DNA sequence
	
	>>> read_fasta_file("A.short.txt")
	('Ashort', 'GATCGTAGAGTGAGACCTAGTGTTTG')
	
	>>> read_fasta_file("B.short.txt")
	('Bshort', 'CTCGTAGGTGAGATTCCTAGTGCC')
	
	>>> read_fasta_file("C.short.txt")
	('Cshort', 'CTCGTAGGTGAGATTCCTAGTGCCCTCGTAGGTGAGATTCCTAGTGCCCTCGTAGGTGAGATTCCTAGTGCCCTCGTAGGTGAGATTCCTAGTGCCCTCGTAGGTGAGATTCCTAGTGCC')
	"""
	
	with open(fasta_file) as fin:
		header = fin.readline()[1:].strip()
		fasta_string = ""
		for line in fin:
			fasta_string += line.strip()
			
	return header, fasta_string
	
def calc_word_code(dnastring):
    """
    dnastring - some dna string 
    
    returns a number representing the decimal value of the dnastring
    
    >>> calc_word_code('ACTA')
    28
    
    >>> calc_word_code('ACTNGTA')
    -1
    
    >>> calc_word_code('TCGAGC')
    3465
    
    >>> calc_word_code('ACATTN')
    -1
    
    >>> calc_word_code('AN')
    -1
    
    >>> calc_word_code('TT')
    15
    
    >>> calc_word_code('CA')
    4
    
    >>> calc_word_code('TC')
    13
    """
    value_map = {   'A': 0,
                    'C': 1,
                    'G': 2,
                    'T': 3
                }
    
    values = []
    for char in dnastring:
        if char not in value_map:
            return -1
        else:
            values.append(value_map[char])
            
    word_code = 0
    for n, value in enumerate(values):
        word_code += value * (4 ** (len(values) - n - 1))

    return word_code 
def create_word_code_array(dnastring, wlength):
    '''
    >>> create_word_code_array('CATCGTGANATCGTGTN', 2)
    [4, 3, 13, 6, 11, 14, 8, -1, -1, 3, 13, 6, 11, 14, 11, -1, -1]
    '''
    word_code_array = []
    for i in range(0, len(dnastring)):
        chunk = dnastring[i: i + wlength]
        code = -1 if len(chunk) < wlength else calc_word_code(chunk)
        word_code_array.append(code)
        
    return word_code_array
if __name__ == "__main__":
	doctest.testmod()