import doctest
from sys import argv

def read_fasta(seq):
	"""
	seq - a .txt file in FASTA format
	
	returns a string representing a DNA sequence
	
	>>> read_fasta("A.short.txt")
	'GATCGTAGAGTGAGACCTAGTGTTTG'
	
	>>> read_fasta("B.short.txt")
	'CTCGTAGGTGAGATTCCTAGTGCC'
	"""
	
	with open(seq) as fin:
		header = fin.readline()
		fasta_string = ""
		for line in fin:
			fasta_string += line.strip()
			
	return fasta_string
	
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
	
def print_matrix(A):
	"""
	A - a 2d matrix
	
	Output - a not-ugly matrix
	"""
	# print('\n\t'.join(['\t'.join(['{}'.format(int(item)) for item in row]) for row in A]))
	with open('output.txt', 'w') as fout:
		print('\n'.join(['\t'.join(['{}'.format(int(item)) for item in row]) for row in A]))
		fout.write('\n'.join(['\t'.join(['{}'.format(int(item)) for item in row]) for row in A]))
		
	
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
		
	>>> calculate_matrices('GATCGTAGAGTGAGACCTAGTGTTTG', 'CTCGTAGGTGAGATTCCTAGTGCC', 10, -20, 40, 2)
	'SO FAR SO GOOD'
		
	"""
	m = len(A)
	n = len(B)
	S, D, I = [[[0 for _ in range(n + 1)] for _ in range(m + 1)] for _ in range(3)]
	
	# S = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
	# D = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
	# I = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
	
	# S, D, I = [np.zeros((m + 1, n + 1)) for _ in range(3)]
	
	# print("S is a {} by {} matrix.".format(len(S), len(S[0])))
	
	
	
	max_score_row = m
	max_score_column = n
	
	gap_init_penalty = gap_open + gap_extend
	
	max_score = 0
	
	D[m][n] = -gap_init_penalty
	I[m][n] = -gap_init_penalty
	
	for i in reversed(range(n)):
		D[m][i] = -gap_init_penalty
		I[m][i] = -gap_init_penalty

                
def traceback(S, D, I, A, B, gap_init_penalty, max_score, max_score_row, max_score_column):
	"""
	A, B - two DNA sequences
	S, D, I - the three matrices returned by calculate_matrices(A, B)
	gap_init_penalty - a positive integer indicating a penalty for opening a gap
   
	
	"""
	# m, n = S.shape
	m = len(S)
	n = len(S[0])
	# print("m -> {}, n -> {}".format(m, n))
	
	opt_align = []
	i = max_score_row
	j = max_score_column
	current_mat = 'S'
	# print("A -> {}\nB -> {}".format(A, B))
	
	while i <= m and j <= n:
		if current_mat == 'S':
			# print("Current matrix is S.")
			# print("i -> {}, m -> {}, j -> {}, n -> {}, S[i][j] -> {}".format(i, m, j, n, S[i][j]))
			# print("S[i][j] = {} is D[i][j] = {} - {}".format(S[i][j], D[i][j], S[i][j] is D[i][j]))
			if i == m - 1 or j == n - 1 or S[i][j] == 0:
				break
			if S[i][j] == D[i][j]:
				# print("Current matrix is D.")
				current_mat = 'D'
				continue
			if S[i][j] == I[i][j]:
				# print("Current matrix is I.")   
				current_mat = 'I'
				continue
			opt_align.append((A[i], B[j]))
			# print("Appended {}".format((A[i], B[j])))
			i += 1
			j += 1
			
		if current_mat == 'D':
			opt_align.append((A[i], ' '))
			# print("Appended {}".format((A[i], ' ')))
			# print("S[i][j] - gap_init_penalty = {} is D[i][j] = {} - {}".format(S[i + 1][j] - gap_init_penalty, D[i][j], S[i + 1][j] - gap_init_penalty is D[i][j]))
			if i == m - 1 or D[i][j] == S[i + 1][j] - gap_init_penalty:
				# print("Current matrix is S.")
				current_mat = 'S'
			i += 1
			continue
			
		if current_mat == 'I':
			opt_align.append((' ', B[j]))
			# print("Appended {}".format((' ', B[j])))
			if j == n - 1 or I[i][j] == S[i][j + 1] - gap_init_penalty:
				# print("Current matrix is S.")
				current_mat = 'S'
			j += 1
			continue
			
	# print(opt_align)
	row_last = i
	col_last = j
	
	# print(S[2][col_last])
	# print(S[row_first][col_last]) CHANGE LATER AHHHH
	
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
		
	print("\t{}\n\t{}\n\t{}".format(aligned_string_A, aligned_string_mid, aligned_string_B))
	return "\t{}\n\t{}\n\t{}".format(aligned_string_A, aligned_string_mid, aligned_string_B)
	
	# return (m, n)
	
def results(mismatch_score, gap_open, gap_extend, ):

	
	print("Match Score\tMismatch Score\tGap-Open Penalty\tGap-Extension Penalty")
	print("---------------------------------------------------------------------")
	print("\t{}\t\t{}\t\t{}\t\t\t{}".format(10, -20, 40, 2))
	
	print("\n     Sequence A: Ashort")
	print("         Length: 26")
	print("     Sequence B: Bshort")
	print("         Length: 24")
	
	print("\nAlignment Score: 104")
	print("         Length: 22")
	print("     Start in A: 3")
	print("     Start in B: 2")
	print("       End in A: 22")
	print("       End in B: 22")
	
def local_alignment(seq1, seq2, mismatch_score, gap_open, gap_extend):
	"""
	seq1, seq2 - .txt files in FASTA format
	mismatch - a negative integer value for mismatch score
	gap_open - a non-negative integer value for gap_open penalty
	gap_extend - a positive integer value for gap extension penalty
	"""
	
	
	MATCH_SCORE = 10
	# read in FASTA strings, return DNA sequences
	A = read_fasta(seq1)
	B = read_fasta(seq2)
	
	# compute dynamic programming matrices
	S, D, I, max_score, max_score_row, max_score_column = calculate_matrices(A, B, MATCH_SCORE, mismatch_score, gap_open, gap_extend)
	
	gap_init_penalty = gap_open + gap_extend
	
	# trace back optimal alignment through matrices
	traceback(S, D, I, A, B, gap_init_penalty, max_score, max_score_row, max_score_column)
	

if __name__ == '__main__':
	# doctest.testmod()
	#calculate_matrices('GATCGTAGAGTGAGACCTAGTGTTTG', 'CTCGTAGGTGAGATTCCTAGTGCC', 10, -20, 40, 2)
	# A = read_fasta("A.short.txt")
	# B = read_fasta("B.short.txt")
	# S, D, I, max_score, max_score_row, max_score_column = calculate_matrices(A, B, 10, -20, 40, 2)
	# gap_init_penalty = 42
	# traceback(S, D, I, A, B, gap_init_penalty, max_score, max_score_row, max_score_column)
	#results()
	print(argv)
	argv[3:] = [int(arg) for arg in argv[3:]]
	print(argv)
	local_alignment(*argv[1:])
	
	
	dna_string1 = read_fasta(seq1)
	dna_string2 = read_fasta(seq2)
	
	# compute dynamic programming matrices
	S, D, I = calculate_matrices(dna_string1, dna_string2, MATCH_SCORE, mismatch_score, gap_open, gap_extend)
	
	# trace back optimal alignment through matrices
	opt_align = traceback(S, D, I)
	
	# print out alignment
	print("Not yet implemented.")
