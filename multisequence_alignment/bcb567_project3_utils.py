'''
    File name: bcb567_project2_utils.py
    Author: Paul Villanueva
    Date created: 11/29/2017
    Date last modified: 11/29/2017
    Python Version: 2.7
'''

import bcb567_project3_classes as bcbclasses
import doctest

def create_superstring(dnastrings):
    """
    input:
        dnastrings - a list of tuples (h, s), with header h and dna sequence s
        
    output:
        a string of the form s_1#s_2#...#s_len(dnastrings)
    
    >>> l = bcbutils.read_fasta_file_new("./test_files/fasta_files/multi_1.txt")
    >>> create_superstring(l)
    'AAC#AAAACA#CATCGTGANATCGTGTN'
    
    >>> l = bcbutils.read_fasta_file_new("./test_files/fasta_files/multi_2.txt")
    >>> create_superstring(l)
    'ANACNCGNGNNTNTTNTN#NANACNCGNGNNNTNTTNT#NANACNCGNGNNTNTTNTNN'
    """
    return '#'.join([dnastring[1] for dnastring in dnastrings])

def process_dna_string(dnastring, wordmodel):
    """
    dnastring - some dna string
    wordmodel - a string of 1s and 0s indicating checked and unchecked positions
    
    returns a new dnastring made by removing char[j] for each j such that 
        wordmodel[j] == 0
    
    >>> process_dna_string('ACTT', '1101')
    'ACT'
    
    >>> process_dna_string('GGATGAT', '1101100')
    'GGTG'
    
    >>> process_dna_string('CGTATGCAA', '010101011')
    'GAGAA'
    
    >>> process_dna_string('TGCACGATCGATGCA', '111111110010110')
    'TGCACGATAGC'
    """
            
    return ''.join([dnastring[i] for i in range(len(dnastring)) if int(wordmodel[i]) == 1])

def read_word_model(word_model_file):
    """
    word_model_file - a .txt file containing a string of 1s and 0s on a single line 
    
    returns a string of 1s and 0s
    >>> read_word_model('./test_files/word_models/word_model_1.txt')
    '1101010100011111111010111010101'
    
    >>> read_word_model('./test_files/word_models/word_model_2.txt')
    '11111101011110001110'
    """
    with open(word_model_file) as fin:
        return fin.readline()


   

def read_fasta_file(fasta_file):
    """
    fasta_file - a .txt file in FASTA format
    
    returns a list containing tuples of (header, dnastring)
    
    >>> read_fasta_file_new("./test_files/fasta_files/multi_1.txt")
    [(seq1, AAC), (seq2, AAAACA), (seq3, CATCGTGANATCGTGTN)]
    
    >>> read_fasta_file_new("./test_files/fasta_files/multi_2.txt")
    [(Seq1, ANACNCGNGNNTNTTNTN), (Seq2, NANACNCGNGNNNTNTTNT), (Seq3, NANACNCGNGNNTNTTNTNN)]
    """
    
    with open(fasta_file) as fin:
    
        dnastring_list = []
        header, seq = None, ""
        
        for line in fin:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    dnastring_list.append((header, seq))
                header, seq = line[1:], ""
            else:
                seq += line.rstrip()
        dnastring_list.append((header, seq))
                    
        return dnastring_list

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
    value_map = {   
                    'A': 0,
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
   
def create_word_code_array(dnastring, wordmodel):
    '''
    >>> create_word_code_array('CATCGTGANATCGTGTN', '11')
    [4, 3, 13, 6, 11, 14, 8, -1, -1, 3, 13, 6, 11, 14, 11, -1, -1]
    
    >>> create_word_code_array('TCGCATCNT', '101')
    [14, 5, 8, 7, 1, -1, 7, -1, -1]
    
    >>> create_word_code_array('ACNTNTCAT', '0101')
    [7, -1, 15, -1, 12, 7, -1, -1, -1]
    '''
    
    wlength = len(wordmodel)
    word_code_array = []
    
    for i in range(len(dnastring) - wlength + 1):
        raw_chunk = dnastring[i: i + wlength]
        processed_chunk = process_dna_string(raw_chunk, wordmodel)
        code = calc_word_code(processed_chunk)
        word_code_array.append(code)
    
    for i in range(wlength - 1):
        word_code_array.append(-1)
    
    return word_code_array
    
   
def create_superword_array(word_code_array, wlength, wlcut):
    """
    inputs:
        word_code_array - an array of word codes returned by create_word_code_array
        wlength - the number of characters in each word 
        wlcut - the number of words in each superword 
        
    Create a sorted array of superwords from word_code_array.
    
    >>> create_superword_array([4, 3, 13, 6, 11, 14, 8, -1, -1, 3, 13, 6, 11, 14, 11, -1, -1], 2, 3)
    [8, 9, 16, 17, 2, 10, 1, 4, 12, 7, 15, 5, 13, 3, 11, 6, 14]
    """
    n = len(word_code_array)
    
    score_bucket = [[] for _ in range(4 ** wlength + 1)]
    sorted_superwords = [i for i in range(n)]

    for word_level in range(wlcut):
        for i in range(n):
            pos = sorted_superwords[i] - wlength
            if pos >= 0:
                code = word_code_array[pos]
                score_bucket[code].append(pos)

        for pos in range(n - wlength, n):
            code = word_code_array[pos]
            score_bucket[code].append(pos)
   
        sorted_superwords = [i for i in score_bucket[-1]]
    
        for bucket in score_bucket[:-1]:
            for pos in bucket:
                sorted_superwords.append(pos)

        score_bucket = [[] for _ in range(4 ** wlength + 1)]
       
    return [i + 1 for i in sorted_superwords]