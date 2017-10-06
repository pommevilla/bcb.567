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

def delta(a, b, match, mismatch):
	"""
	a, b - two characters 
	match - a positive integer score for a character match
	mismatch - a negative integer score for a character mismatch
	
	>>> delta("C", "C", 10, -20)
	10
	
	>>> delta("C", "T", 10, -30)
	-30
	
	"""
	return match if a is b else mismatch

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
	
def calculate_matrices(A, B, MATCH_SCORE, mismatch_score, gap_open, gap_extend):
	"""
	A, B - two DNA sequences 
	mismatch - a negative integer indicating mismatch scare
	gap_open - a non-negative integer indicating gap_open penalty
	gap_extend - a positive integer indicating gap extension penalty
	
	returns:
			matrices of integers S, D, and I, where:
				-S[i][j] is the maximum score of all local alignments starting at A[i + 1] and B[j + 1]
				-D[i][j] is the maximum score of those local alignments that begin with a deletion gap
				-I[i][j] is the maximum score of those local alignments that begin with an insertion gap
			integers max_score, max_score_row, max_score_column, where:
				-max_score is the highest alignment score in the S matrix
				-max_score_row and max_score_column are max_score's row and column position in the S matrix, respectively
		
	"""
	m = len(A)
	n = len(B)
	S, D, I = [[[0 for _ in range(n + 1)] for _ in range(m + 1)] for _ in range(3)]
	
	max_score_row = m
	max_score_column = n
	
	gap_init_penalty = gap_open + gap_extend
	
	max_score = 0
	
	D[m][n] = -gap_init_penalty
	I[m][n] = -gap_init_penalty
	
	for i in reversed(range(n)):
		D[m][i] = -gap_init_penalty
		I[m][i] = -gap_init_penalty

		
	for i in reversed(range(m)):
		D[i][n] = -gap_init_penalty
		I[i][n] = -gap_init_penalty
		
		for j in reversed(range(n)):
			D[i][j] = max(D[i + 1][j] - gap_extend,
							S[i + 1][j] - gap_init_penalty)
			I[i][j] = max(I[i][j + 1] - gap_extend, 
							S[i][j + 1] - gap_init_penalty)
			S[i][j] = max(0, S[i + 1][j + 1] + delta(A[i], B[j], MATCH_SCORE, mismatch_score),
							D[i][j], I[i][j])
			if max_score < S[i][j]:
				max_score = S[i][j]
				max_score_row = i
				max_score_column = j 
	
	return S, D, I, max_score_row, max_score_column
	
def traceback(S, D, I, A, B, gap_init_penalty, max_score_row, max_score_column):
	"""
	A, B - two DNA sequences
	S, D, I - the three matrices returned by calculate_matrices(A, B)
	gap_init_penalty - a positive integer indicating a penalty for opening a gap.  This quantity is equal to the sum gap_open_penalty + gap_extend_penalty that is subtracted in the traceback algorithm.
	max_score_row, max_score_column - integers indicating the location of the max_score in the S matrix
	
	"""

	m = len(S)
	n = len(S[0])

	opt_align = []
	i = max_score_row
	j = max_score_column
	current_mat = 'S'

	while i <= m and j <= n:
		if current_mat == 'S':
			if i == m - 1 or j == n - 1 or S[i][j] == 0:
				break
			if S[i][j] == D[i][j]:
				current_mat = 'D'
				continue
			if S[i][j] == I[i][j]:
				current_mat = 'I'
				continue
			opt_align.append((A[i], B[j]))
			i += 1
			j += 1
			
		if current_mat == 'D':
			opt_align.append((A[i], ' '))
			if i == m - 1 or D[i][j] == S[i + 1][j] - gap_init_penalty:
				current_mat = 'S'
			i += 1
			continue
			
		if current_mat == 'I':
			opt_align.append((' ', B[j]))
			if j == n - 1 or I[i][j] == S[i][j + 1] - gap_init_penalty:
				current_mat = 'S'
			j += 1
			continue
			

	row_last = i
	col_last = j
	
	aligned_string_A = ''
	aligned_string_mid = ''
	aligned_string_B = ''
	
	for pair in opt_align:
		aligned_string_A += pair[0]
		aligned_string_B += pair[1]
		if pair[0] == ' ' or pair[1] == ' ':
			aligned_string_mid += '-'
		else:
			aligned_string_mid += "|"
			
	return opt_align
		
if __name__ == "__main__":
	doctest.testmod()