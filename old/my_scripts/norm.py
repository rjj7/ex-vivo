#!/usr/bin/env python
import argparse

import nibabel as nb
import numpy as np


parser = argparse.ArgumentParser(
        description='Normalizes image intensity to 0-1.')
parser.add_argument('input', metavar='input', help='Input image')
parser.add_argument('output', metavar='output', 
                    help='Normalized output image.')


args = parser.parse_args()

img = nb.load(args.input)
data = img.get_data()
affine = img.get_affine()

data=data.astype('float64')
data_norm = np.zeros (data.shape)

for ii in np.ndindex (data.shape):
    m = data[ii]
    m= m - (np.abs(m).min())
    m = m / (np.abs(m).max() - np.abs(m).min())
    data_norm[ii] = m

nb.save(nb.Nifti1Image(data_norm.astype(np.float32),
                                   img.affine, img.header), args.output)

