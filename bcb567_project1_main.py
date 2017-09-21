'''
    File name: bcb567_project1_main.py
    Author: Paul Villanueva
    Date created: 9/13/2017
    Date last modified: 9/20/2017
    Python Version: 3.6
'''

import bcb567_project1_classes as bcbclasses
import bcb567_project1_utils as bcbutils
from sys import argv

def local_alignment_driver(args):
	"""
	args - user input of the form "file1, file2, mismatch_score, gap_open_penalty, gap_extend_penalty"
	
	prints out local alignment information as specified in assignment.
	"""
	# Default match score as specified in assignment.
	MATCH_SCORE = 10
	
	A = bcbclasses.DNAString(*bcbutils.read_fasta_file(args[1]))
	B = bcbclasses.DNAString(*bcbutils.read_fasta_file(args[2]))
	mismatch_score = int(args[3])
	gap_open_penalty = int(args[4])
	gap_extend_penalty = int(args[5])
	
	AB_optimal_alignment = bcbclasses.LocalAlignment(A, B, MATCH_SCORE, mismatch_score, gap_open_penalty, gap_extend_penalty)
	return AB_optimal_alignment
	
if __name__ == "__main__":
	optimal_alignment = local_alignment_driver(argv)
	print(optimal_alignment)