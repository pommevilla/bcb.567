'''
    File name: bcb567_project1_classes.py
    Author: Paul Villanueva
    Date created: 9/13/2017
    Date last modified: 9/20/2017
    Python Version: 3.6
'''

import doctest
import bcb567_project1_utils as bcbutils

class DNAString:
	"""
	Represents a DNA sequence with a header.  
	"""

	def __init__(self, header, seq):
		self.header = header
		self.seq = seq
	   
	@property
	def length(self):
		"""
		>>> string1 = bcbutils.read_fasta_file("A.short.txt")
		>>> A = DNAString(*string1)
		>>> A.length
		26
		"""
		return len(self.seq)
		
	def __len__(self):
		return len(self.seq)

	def __iter__(self):
		self.current = 0
		return self
		
	def __next__(self):
		if self.current == self.length - 1:
			raise StopIteration
		else:
			self.current = self.current + 1
			return self.seq[self.current]
			
	def __getitem__(self, i):
		"""
		>>> string1 = bcbutils.read_fasta_file("A.short.txt")
		>>> A = DNAString(*string1)
		>>> print(A[0])
		G
		"""
		return self.seq[i]
	
	def __str__(self):
		"""
		>>> string1 = bcbutils.read_fasta_file("A.short.txt")
		>>> A = DNAString(*string1)
		>>> print(A)
		Ashort
		GATCGTAGAGTGAGACCTAGTGTTTG
		"""
		return "{}\n{}".format(self.header, self.seq)

class LocalAlignment:
	"""
	Represents an optimal local alignment between two DNA strings.
	"""
	
	def __init__(self, dnastring1, dnastring2, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty):
		self.dnastring1 = dnastring1 
		self.dnastring2 = dnastring2
		
		self.match_score = match_score 
		self.mismatch_score = mismatch_score 
		self.gap_open_penalty = gap_open_penalty
		self.gap_extend_penalty = gap_extend_penalty
		
		
		S_matrix, D_matrix, I_matrix, self.max_score_row, self.max_score_column = bcbutils.calculate_matrices(self.dnastring1, self.dnastring2, self.match_score, mismatch_score, gap_open_penalty, gap_extend_penalty)
		
		self.max_alignment_score = S_matrix[self.max_score_row][self.max_score_column]
		
		self.optimal_alignment = bcbutils.traceback(S_matrix, D_matrix, I_matrix, self.dnastring1, self.dnastring2, self.gap_open_penalty + self.gap_extend_penalty, self.max_score_row, self.max_score_column)
		
		del S_matrix, D_matrix, I_matrix
		
		
	@property
	def length(self):
		"""
		>>> string1 = bcbutils.read_fasta_file("A.short.txt")
		>>> A = DNAString(*string1)
		>>> string2 = bcbutils.read_fasta_file("B.short.txt")
		>>> B = DNAString(*string2)
		>>> AB_local_alignment = LocalAlignment(A, B, 10, -20, 40, 2)
		>>> AB_local_alignment.length
		22
		"""
		return len(self.optimal_alignment)
		
	def process_alignment(self):
		"""
		Partitions self.optimal_alignment into subsequences of size 70 and saves the chunks into a list of strings.

		Returns formatted sections of the alignment ready for printing to the console and counts of deletions, insertions, matches
		mismatches, and the position in each string where the alignments end.
		
		This one is really ugly and got away from me because I was working on this one last.  
		"""
		insertion_count = 0
		deletion_count = 0
		match_count = 0
		mismatch_count = 0
		
		# Breaks the optimal alignment into subsequences 70 pairs long per assignment instructions.
		optimal_alignment_sections_raw = [section for section in bcbutils.sectioner(self.optimal_alignment)]
		
		optimal_alignment_sections_processed = []
		
		
		align_position_A = self.max_score_row + 1
		align_position_B  = self.max_score_column + 1
		
		# n counts our spot in the alignment and puts a '.' or ':' every 5 spots, alternating, to ease reading.
		n = 0
			
		for section in optimal_alignment_sections_raw:
			spacer_string = '{:5} '.format(' ')
			aligned_string_A = '{:>5} '.format(n + align_position_A)
			aligned_string_mid = '{:5} '.format(' ')
			aligned_string_B = '{:>5} '.format(n + align_position_B)
			for i in section:
				if ((n + 1) % 10)  == 5:
					spacer_string += '.'
				elif ((n + 1) % 10) == 0:
					spacer_string += ':'
				else:
					spacer_string += ' ' 
				aligned_string_A += i[0]
				aligned_string_B += i[1]
				if i[0] == ' ':
					aligned_string_mid += "-"
					insertion_count += 1
					align_position_A -= 1
				elif i[1] == ' ':
					aligned_string_mid += "-"
					deletion_count += 1
					align_position_B -= 1
				elif i[0] != i[1]:
					aligned_string_mid += 'x'
					mismatch_count += 1
				else:
					aligned_string_mid += "|"
					match_count += 1
					
				n += 1
			optimal_alignment_sections_processed.append("{:>7}\n{:>7}\n{:>7}\n{:>7}".format(spacer_string, aligned_string_A, aligned_string_mid, aligned_string_B))
			
			A_end = n + align_position_A - 1
			B_end = n + align_position_B - 1
		
			
		return optimal_alignment_sections_processed, A_end, B_end, insertion_count, deletion_count, mismatch_count, match_count
		
	def __str__(self):
		sectioned_string_alignment, A_end, B_end, insertion_count, deletion_count, mismatch_count, match_count = self.process_alignment()	
		
		
		# Extra blank character at the end of the argument to make the spacing nice.
		score_info_header = "{:^15}{:^20}{:^20}{:^28}".format(
			"Match Score", "Mismatch Score", "Gap Open Penalty", "Gap Extension Penalty"
			)
		score_info_values = "{:^15}{:^20}{:^20}{:^28}".format(
			self.match_score, self.mismatch_score, self.gap_open_penalty, self.gap_extend_penalty
			)
		
		score_info = "\n".join([score_info_header, score_info_values, ''])
		
		
		# Have mercy on me for the awful spacing.
		seq_info = "{:>20}: {}\n{:>20}: {}\n{:>20}: {}\n{:>20}: {}\n".format(
			"Sequence A", self.dnastring1.header, 
			"Length", self.dnastring1.length, 
			"Sequence B", self.dnastring2.header, 
			"Length", self.dnastring2.length
			)
			
		align_info = "{:>20}: {}\n{:>20}: {}\n{:>20}: {}\n{:>20}: {}\n{:>20}: {}\n{:>20}: {}\n".format(
			"Alignment Score", self.max_alignment_score, 
			"Length", self.length, 
			"Start in A", str(self.max_score_row + 1), 
			"Start in B", str(self.max_score_column + 1), 
			"End in A", A_end, 
			"End in B", B_end
			)
			
		match_info = "{:>20}: {}\n{:>20}: {}\n{:>20}: {:.6f}\n\n{:>20}: {}\n{:>20}: {}\n{:>20}: {}\n".format(
			"Number of matches", match_count, 
			"Number of mismatches", mismatch_count, 
			"Match percentage", match_count / self.length, 
			"Number of deletions", deletion_count, 
			"Number of insertions", insertion_count, 
			"Total length of gaps", deletion_count +insertion_count
			)
		
		
		
		return "\n".join([score_info, seq_info, align_info, match_info, *sectioned_string_alignment])
		
if __name__ == "__main__":
	doctest.testmod()
	


	
		