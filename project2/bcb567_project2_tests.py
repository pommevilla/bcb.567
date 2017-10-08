import unittest
import bcb567_project2_classes as bcbclasses
import bcb567_project2_utils as bcbutils

class SuperwordCodeTesting(unittest.TestCase):
    def setUp(self):
        self.a = bcbclasses.SuperwordCode(3, (2, 3, 4))
        self.b = bcbclasses.SuperwordCode(3, (2, 3, 4))
        self.c = bcbclasses.SuperwordCode(7, (-1, 5, 10))
        self.d = bcbclasses.SuperwordCode(9, (-1, 3, 10))
        self.e = bcbclasses.SuperwordCode(6, (-1, -1, 3))
        self.f = bcbclasses.SuperwordCode(10, (-1, 4, 3))
        self.x = (3, 2, 1) 
        
    def test_superwordcode_eq(self):
        self.assertTrue(self.a == self.b)
        self.assertFalse(self.a != self.b)
        self.assertFalse(self.a == self.c)
        with self.assertRaises(AssertionError):
            self.assertTrue(self.a == self.x)
        
    def test_superwordcode_neq(self):
        self.assertFalse(self.a != self.b)
        self.assertTrue(self.c != self.d)
        self.assertTrue(self.a != self.c)
        self.assertTrue(self.b != self.c)
        with self.assertRaises(AssertionError):
            self.assertTrue(self.b != self.x)
            
    # def test_superwordcode_lt(self):
        # self.assertTrue(self.c < self.d)
        # self.assertFalse(self.a < self.f)
        # self.assertTrue(self.e < self.b)
        # self.assertTrue(self.c < self.f)
        # with self.assertRaises(AssertionError):
            # self.assertTrue(self.b < self.x)
            
            
if __name__ == '__main__':
    unittest.main()