#!/usr/bin/env python
import argparse


import numpy as np
import nibabel as nib
import os
import glob 

from nibabel.streamlines import save as save_trk
from nibabel.streamlines import Tractogram

from dipy.core.gradients import gradient_table, GradientTable
from dipy.data import get_sphere, Sphere, HemiSphere
from dipy.direction import peaks_from_model, ProbabilisticDirectionGetter
from dipy.io import read_bvals_bvecs
from dipy.io.image import save_nifti
from dipy.reconst.peaks import reshape_peaks_for_visualization
from dipy.segment.mask import median_otsu
from dipy.reconst.shm import sh

parser = argparse.ArgumentParser(
        description='Computing Generalized Q-ball Imaging (GQI)')
parser.add_argument('dwi_data', help='DWI Data')
parser.add_argument('bvecs', help='Bvecs file')
parser.add_argument('bvals', help='Bvals file')
parser.add_argument('--mask', help='Binary mask')
parser.add_argument('output_ODF', help='ODFs')
parser.add_argument('output_peaks', help='Peaks for visualization')
parser.add_argument('output_gfa', help='Generalized Fractional Anisotropy map')

args = parser.parse_args()

# Import DWI data, bvec, bval, and binary mask
fimg =  (args.dwi_data)
img = nib.load(fimg)
data = img.get_data()
print("Data Shape: " + str(data.shape))
 
# Bval and bvec file information import
fbval= (args.bvals)
fbvec = (args.bvecs)
bvals, bvecs = read_bvals_bvecs(fbval, fbvec)
gtab = gradient_table(bvals,bvecs)

# Read the voxel size from the image header.
voxel_size = img.header.get_zooms()[:3]
print ('Voxel Size: ' + str(voxel_size))

# Mask 
#add if mask none
img = nib.load(args.mask)
mask = img.get_data()

sphere = get_sphere('symmetric724')

#GQI
l_values = np.sqrt(gtab.bvals * 0.01506)
tmp=np.tile(l_values, (3,1))
gradsT = gtab.bvecs.T
b_vector = gradsT * tmp
b_vector = b_vector.T
gqi_vector = np.real(np.sinc(np.dot(b_vector, sphere.vertices.T)* 0.3/np.pi))
odfs = np.dot(data, gqi_vector)

#min/max normalization
ijk = np.ascontiguousarray(np.array(np.nonzero(seed_mask)).T)
shape = data.shape[:-1]
odfs_norm = np.zeros ((shape + (len(sphere.vertices),)))

for (k, center) in enumerate(ijk):
    m = odfs[tuple(center.astype(np.int))].copy()
    m = m - (np.abs(m).min())
    m = m / (np.abs(m).max())
    odfs_norm[tuple(center.astype(np.int))] = m

npeaks = 3
relative_peak_threshold=.5
min_separation_angle=25
shape=mask.shape    
peak_dirs = np.zeros ((shape + (npeaks, 3)))
peak_values = np.zeros ((shape + (npeaks,)))
peak_indices = np.zeros((shape + (npeaks,)), dtype='int')

odfs_norm =odfs_gqi.astype(float)

for idx in ndindex(shape):
    if not mask[idx]:
        continue
    direction, pk, indices = peak_directions(odfs_norm[idx], sphere, relative_peak_threshold, min_separation_angle)
    n = min(npeaks, pk.shape[0])
    
    peak_dirs[idx][:n] = direction[:n]
    peak_values[idx][:n] = pk[:n]
    peak_indices[idx][:n] = indices[:n]
    #nib.save(nib.Nifti1Image(reshape_peaks_for_visualization(peak_dirs),
                                     img.affine),os.path.join(sub, 'peaks)gqi.nii.gz'))

nib.save(nib.Nifti1Image(odfs_norm.astype(np.float32),
                                  imgm.affine, imgm.header), 'I43/slab2_03_01/odfs_gqi_l03.nii.gz')
#look for numpy.menmap for memory stuff
