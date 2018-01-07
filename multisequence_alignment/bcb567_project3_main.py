'''
    File name: bcb567_project2_main.py
    Author: Paul Villanueva
    Date created: 11/29/2017
    Date last modified: 11/29/2017
    Python Version: 2.7
'''

import bcb567_project3_classes as bcbclasses
import bcb567_project3_utils as bcbutils
from sys import argv



def multisequence_alignment_driver(args):
    '''
    inputs:
        args - user input of the form "seqs_file, word_model, wlcut", where
        seqs_file - a .txt file in fasta format 
        word_model - a .txt file containing a string of 1s and 0s on a single line
        wlcut - an integer
	
	returns a SuperwordArray object create from the inputs
    '''
    dnastrings = bcbutils.read_fasta_file(args[1])
    word_model = bcbutils.read_word_model(args[2])
    wlcut = int(args[3])

    return bcbclasses.MultiSequenceAlignment(dnastrings, word_model, wlcut)

if __name__ == "__main__":
    multisequence_alignment_solution = multisequence_alignment_driver(argv)
    print multisequence_alignment_solution
    
    
