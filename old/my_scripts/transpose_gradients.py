#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(
        description='Transpose bvec file')

import numpy as np

parser.add_argument('in_bvec', help='Input bvec file')
parser.add_argument('out_bvec', help='Output bvec file')

args = parser.parse_args()

bvecs = np.loadtxt(args.in_bvec)
trans_bvec = bvecs.transpose()
np.savetxt(args.out_bvec, trans_bvec, fmt='%1.16f')
