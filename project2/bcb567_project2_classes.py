'''
    File name: bcb567_project2_utils.py
    Author: Paul Villanueva
    Date created: 10/7/2017
    Date last modified: 10/7/2017
    Python Version: 2.7
'''

import doctest
import bcb567_project2_utils as bcbutils

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

class SuperwordCode:

    def __init__(self, pos, swcode):
        self.pos = pos 
        self.swcode = swcode 
        
    def __eq__(self, other):
        """
        >>> a = SuperwordCode(3, (2, 3, 4))
        >>> b = SuperwordCode(3, (2, 3, 4))
        >>> a == b 
        True
        
        >>> x = (3, 2, 1)
        >>> a == x 
        Traceback (most recent call last):
        ...
        AssertionError
        
        >>> c = SuperwordCode(7, (-1, 5, 10))
        >>> d = SuperwordCode(9, (-1, 3, 10))
        >>> c == d
        False
        """
        assert isinstance(other, SuperwordCode)
        return self.pos == other.pos and self.swcode == other.swcode
        
    def __lt__(self, other):
        """
        >>> a = SuperwordCode(3, (2, 3, 4))
        >>> b = SuperwordCode(3, (2, 3, 4))
        >>> a < b 
        False
        
        >>> x = (3, 2, 1)
        >>> a < x 
        Traceback (most recent call last):
        ...
        AssertionError
        
        >>> c = SuperwordCode(7, (-1, 5, 10))
        >>> d = SuperwordCode(9, (-1, 3, 10))
        >>> d < c
        True
        """
        assert isinstance(other, SuperwordCode)
        if -1 not in self.swcode and -1 not in other.swcode:
            if self.swcode == other.swcode:
                return self.pos < other.pos
            else:
                return self.swcode < other.swcode
        
    # def __le__(self, other):
    
    # def __ge__(self, other):
    
    # def __gt__(self, other):
    
    def __ne__(self, other):
        """
        >>> a = SuperwordCode(3, (2, 3, 4))
        >>> b = SuperwordCode(3, (2, 3, 4))
        >>> a != b 
        False
        
        >>> x = (3, 2, 1)
        >>> a != x 
        Traceback (most recent call last):
        ...
        AssertionError
        
        >>> c = SuperwordCode(7, (-1, 5, 10))
        >>> d = SuperwordCode(9, (-1, 3, 10))
        >>> c != d
        True
        """
        assert isinstance(other, SuperwordCode)
        return self.pos != other.pos or self.swcode != other.swcode
        
        
if __name__ == "__main__":
    doctest.testmod()
    # pass