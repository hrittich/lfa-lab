import unittest
from lfa_lab.stencil import *

class StencilTest(unittest.TestCase):

    def test_element_store(self):
        s = SparseStencil([((0,), 1.0)])
        self.assertEqual(s[0], ((0,), 1.0))

        s = SparseStencil([((0,), 1.0j)])
        self.assertEqual(s[0], ((0,), 1.0j))

    def test_multiply(self):

        # test whether to applications of the approximation of the first
        # derivative approximation gives the second one
        S1 = SparseStencil([((0,), -1),
                            ((1,),  1)])

        S2 = SparseStencil([((-1,), -1),
                            (( 0,),  1)])

        es = { o:v for o, v in (S1 * S2) }

        self.assertAlmostEqual( es[(-1,)],  1)
        self.assertAlmostEqual( es[( 0,)], -2)
        self.assertAlmostEqual( es[( 1,)],  1)
        self.assertEqual(len(es), 3)


if __name__ == '__main__':
    unittest.main()

