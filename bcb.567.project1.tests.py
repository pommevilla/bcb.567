import numpy as np
import unittest
import scratch as bcb1

class local_alignment_tests(unittest.TestCase):

	def test_delta(self):
		self.assertTrue(bcb1.delta("C", "C", 10, -20) == 10)
		self.assertTrue(bcb1.delta("C", "T", 10, -30) == -30)
		self.assertFalse(bcb1.delta("G", "T", 10, -30) == -75)
		
	def test_read_fasta(self):
		self.assertTrue(bcb1.read_fasta("A.short.txt") == 'GATCGTAGAGTGAGACCTAGTGTTTG')
		self.assertTrue(bcb1.read_fasta("B.short.txt") == 'CTCGTAGGTGAGATTCCTAGTGCC')
	   

	def test_traceback(self):
		A = bcb1.read_fasta("A.short.txt")
		B = bcb1.read_fasta("B.short.txt")
		S, D, I, max_score, max_score_row, max_score_column = bcb1.calculate_matrices(A, B, 10, -20, 40, 2)
		gap_init_penalty = 42
		expected = "\t{}\n\t{}\n\t{}".format("TCGTAGAGTGAGA  CCTAGTG", "||||||-||||||--|||||||", "TCGTAG GTGAGATTCCTAGTG")
		
		self.assertTrue(bcb1.traceback(S, D, I, A, B, gap_init_penalty, max_score, max_score_row, max_score_column) == expected)

		
		
def main():
	unittest.main()
	
if __name__ == '__main__':
	main()
