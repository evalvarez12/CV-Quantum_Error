# -*- coding: utf-8 -*-
"""
Unittest fucntions for covariance_matrix

Created on Thu Mar 14 13:29:56 2019

@author: Eduardo Villasenor
"""

import unittest
import numpy as np
import covariance_matrix2 as cm

class TestCovarianceMatrixMethods(unittest.TestCase):

    def test_kron(self):
        N = 5
        a = np.random.rand(N, N)
        b = np.random.rand(N, N)
        c = np.random.rand(N, N)

        test_assert = np.kron(a, np.kron(b, c))
        test = cm.kron([a, b, c])

        np.testing.assert_array_equal(test, test_assert)


    def test_kron_sum(self):
        N = 2
        a = np.random.rand(N, N)
        b = np.random.rand(N, N)
        c = np.random.rand(N, N)

        eye = np.identity(N)

        test_assert = np.kron(a, np.kron(eye, eye)) + \
                      np.kron(eye, np.kron(b, eye)) + \
                      np.kron(eye, np.kron(eye, c))


        test = cm.kron_sum([a, b, c], [0, 1, 2], 3)

        np.testing.assert_array_equal(test, test_assert)


if __name__ == '__main__':
    unittest.main()
