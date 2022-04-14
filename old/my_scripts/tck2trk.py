#!/usr/bin/env python

import argparse

import os
import nibabel as nb
import numpy as np

parser = argparse.ArgumentParser(
        description='Converts tractography file from .tck to .trk')
parser.add_argument('in_tract_file', metavar='Input file',
                   help='Tractography file name. Format must be \n'
                        'readable by Nibabel.')
parser.add_argument('out_tract_name', metavar='Output name',
                   help='Output filename. Format must be \n'
                        'writable by Nibabel.')
parser.add_argument('reference', metavar='Reference image',
		   help='Reference file used to create header, must be .nii.gz')

args = parser.parse_args()

input_tract = nb.streamlines.load(args.in_tract_file)
streamlines = input_tract.streamlines

ref_img = nb.load(args.reference)
new_header = nb.streamlines.TrkFile.create_empty_header()
new_header[nb.streamlines.Field.VOXEL_SIZES] = tuple(ref_img.header.
                                                          get_zooms())[:3]
new_header[nb.streamlines.Field.DIMENSIONS] = tuple(ref_img.shape)[:3]
new_header[nb.streamlines.Field.VOXEL_TO_RASMM] = (ref_img.header.
                                                        get_best_affine())
affine = new_header[nb.streamlines.Field.VOXEL_TO_RASMM]
new_header[nb.streamlines.Field.VOXEL_ORDER] = ''.join(
           nb.aff2axcodes(affine))

new_tractogram = nb.streamlines.Tractogram(streamlines,
                                            affine_to_rasmm=np.eye(4))
tractogram_type = nb.streamlines.detect_format(args.out_tract_name)
fileobj = tractogram_type(new_tractogram,
                                  header=new_header)
nb.streamlines.save(fileobj, args.out_tract_name)

