#!/usr/bin/env python
import argparse

import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
import os
import glob
import pickle
from collections import defaultdict
from dipy.core.ndindex import ndindex
from dipy.data import get_sphere, HemiSphere, Sphere
from dipy.io.image import save_nifti
from dipy.io.streamline import load_trk
from dipy.io.streamline import save_trk
from dipy.tracking import utils
from dipy.tracking.streamline import Streamlines
from dipy.tracking.utils import random_seeds_from_mask, move_streamlines, target, density_map, connectivity_matrix
from dipy.tracking._utils import (_mapping_to_voxel, _to_voxel_coordinates)


parser = argparse.ArgumentParser(
        description='computing connectivity matrix')
parser.add_argument('labels_file', help='int')
parser.add_argument('track_file', help='tractography file')
parser.add_argument('output_matrix_py', help='matrix in .npy')
parser.add_argument('output_matrix_txt', help='matrix in .txt')
parser.add_argument('output_matrix_pkl', help='matrix in .pkl')
parser.add_argument('output_matrix_jpg', help='matrix in .jpg')

args = parser.parse_args()

img = nib.load(args.labels)
labels = img.get_data()
labels = labels.clip(min=0)
labels = labels.astype('int32')

tracts, hdr = load_trk(args.track_file)
streamlines = Streamlines(tracts)
streamlines = move_streamlines(streamlines, np.eye(4), img.affine)
matrix, mapping = connectivity_matrix(streamlines, labels, affine=np.eye(4), return_mapping=True, symmetric=True, mapping_as_streamlines=True)
np.save(args.output_matrix_py, matrix)
np.savetxt(args.output_matrix_txt, matrix, delimiter=" ", fmt='%1.3f')
pkl_file = open(args.output_matrix_pkl, "wb")
pickle.dump(mapping, pkl_file)
pkl_file.close()
matrix[:1, :] = 0
matrix[:, :1] = 0
np.fill_diagonal(matrix, 0)
labels = np.loadtxt('lut.txt', dtype = 'string', delimiter=' ', usecols = (1))
myplot = plt.imshow(np.log1p(matrix), interpolation='nearest', cmap='coolwarm')
plt.colorbar(myplot, spacing='proportional')
plt.yticks(np.arange(29), labels, fontsize=8)
plt.xticks(np.arange(29), labels, fontsize=8, rotation=90)
plt.savefig(args.output_matrix_jpg, bbox_inches='tight')
