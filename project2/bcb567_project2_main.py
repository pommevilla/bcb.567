'''
    File name: bcb567_project2_main.py
    Author: Paul Villanueva
    Date created: 10/7/2017
    Date last modified: 10/7/2017
    Python Version: 2.7
'''

import bcb567_project2_classes as bcbclasses
import bcb567_project2_utils as bcbutils
from sys import argv

def superword_driver(args):
	"""
	args - user input of the form "seq_file word_model wlcut"
	
	prints out local alignment information as specified in assignment.
	"""
	
	dnastring = bcbclasses.DNAString(*bcbutils.read_fasta_file(args[1]))
    
	wlcut = int(args[3])
	
	
if __name__ == "__main__":
	optimal_alignment = local_alignment_driver(argv)
	print(swarray)