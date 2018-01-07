'''
    File name: bcb567_project2_tests.py
    Author: Paul Villanueva
    Date created: 10/7/2017
    Date last modified: 10/26/2017
    Python Version: 2.7
'''

import unittest
import bcb567_project2_classes as bcbclasses
import bcb567_project2_utils as bcbutils

class SuperwordTesting(unittest.TestCase):
    def setUp(self):
        word_model = bcbutils.read_word_model('word_model_3.txt')
        dna_string = bcbutils.read_fasta_file('D.short.txt')
        wlcut = 3
        self.superword_array_1 = bcbclasses.SuperwordArray(dna_string, word_model, wlcut)
     
    def test_word_code_array(self):
        self.assertTrue(self.superword_array_1.word_code_array == [4, 3, 13, 6, 11, 14, 8, -1, -1, 3, 13, 6, 11, 14, 11, -1, -1])
        
    def test_largest_block_info(self):
        self.assertTrue(self.superword_array_1.max_block_start == 5)
        self.assertTrue(self.superword_array_1.max_block_end == 6)
        self.assertTrue(self.superword_array_1.max_block_size == 2)
        self.assertTrue(self.superword_array_1.max_block_superword == (3, 6, 14))
        
        
    def test_get_superword_code(self):
        self.assertTrue(self.superword_array_1.get_superword_code(1) == (3, 6, 14))
        self.assertTrue(self.superword_array_1.get_superword_code(16) == (-1, -1, -1))
        self.assertTrue(self.superword_array_1.get_superword_code(9) == (3, 6, 14))
        self.assertTrue(self.superword_array_1.get_superword_code(13) == (14, -1, -1))        
            
if __name__ == '__main__':
    unittest.main(verbosity = 2)