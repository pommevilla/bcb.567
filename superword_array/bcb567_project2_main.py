'''
    File name: bcb567_project2_main.py
    Author: Paul Villanueva
    Date created: 10/7/2017
    Date last modified: 10/26/2017
    Python Version: 2.7
'''

import bcb567_project2_classes as bcbclasses
import bcb567_project2_utils as bcbutils
from sys import argv



def superword_array_driver(args):
    '''
    inputs:
        args - user input of the form "file1, file2, wlcut", where
        file1 - a .txt file in fasta format 
        file2 - a .txt file containing a string of 1s and 0s on a single line
        wlcut - an integer
	
	returns a SuperwordArray object create from the inputs
    '''
    dnastring = bcbutils.read_fasta_file(args[1])
    word_model = bcbutils.read_word_model(args[2])
    wlcut = int(args[3])

    return bcbclasses.SuperwordArray(dnastring, word_model, wlcut)

if __name__ == "__main__":
    superword_array_solution = superword_array_driver(argv)
    print(superword_array_solution)
    
    
