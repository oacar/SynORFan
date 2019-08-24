import unittest
from bioconductor import get_frame_mapping
from bioconductor import union_range, intersect_range, count_identical_chars
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
class MyTestCase(unittest.TestCase):
    def test_get_frame_mapping(self):
        seq = 'A--TTTGGATGGGG-G-GGGCC'
              #0111201201201220012012
        expected = {0:0,1:1,2:1,3:1,4:2,5:0,6:1,7:2,8:0,
                       9:1,10:2,11:0,12:1,13:2,14:2,15:0,16:0,
                        17:1,18:2,19:0,20:1,21:2}
        self.assertEqual(get_frame_mapping(seq,8), expected)
    def test_range_union_intersection(self):
        range1 = [4, 10]
        range2 = [2, 8]
        un_range = union_range(range1, range2)
        int_range = intersect_range(range1, range2)

        self.assertEqual(un_range, [2, 10])
        self.assertEqual(int_range, [4, 8])
    def test_count_identical_chars(self):
        self.assertEqual(count_identical_chars('MTG--GCR','M-G-GGRC'),3)
if __name__ == '__main__':
    unittest.main()
