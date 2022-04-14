#!/usr/bin/env python
import argparse
import nibabel as nib
import numpy as np

from dipy.core.gradients import gradient_table, GradientTable
from dipy.data import get_sphere, Sphere, HemiSphere
from dipy.direction import (peaks_from_model, ProbabilisticDirectionGetter)
from dipy.io import read_bvals_bvecs
from dipy.io.image import save_nifti
from dipy.io.streamline import save_trk, load_trk
from dipy.tracking.local import (LocalTracking, ThresholdTissueClassifier, ActTissueClassifier, CmcTissueClassifier)
from dipy.tracking.utils import (random_seeds_from_mask) 
from dipy.tracking.streamline import Streamlines

parser = argparse.ArgumentParser(
        description='Tracking from PMF using tissue classifier')
parser.add_argument('odfs', help= 'odfs sf format')
parser.add_argument('bvecs', help='Bvecs file, in FSL format')
parser.add_argument('bvals', help='Bvals file, in FSL format')
parser.add_argument('output_trk', help='Output filename trk')
parser.add_argument('--fa_mask', dest='fa_mask', default=None, help='FA mask')
parser.add_argument('--seed_mask', dest='seed_mask', help='Binary mask to randomly seed streamlines from')
#parser.add_argument('--max_angle', dest='max_angle', default=3, help='Maximum angle between directions, int. Default=20')
#parser.add_argument('--fa_thr', dest='fa_thr', default=.4, help='FA threshold for tissue classifier, float. Default=.4')
#parser.add_argument('--seeds_count', dest='seeds_count', default='2', help= 'Number of seeds/voxel, int. Default=2')
#parser.add_argument('--step_size', dest='step_size', default='2', help= 'Step Size, int. Default=0.175')

args = parser.parse_args()

fimg =  (args.fa_mask)
img = nib.load(fimg)
fa = img.get_data()
sphere = get_sphere('symmetric724')

fimg = (args.odfs)
img2 = nib.load(fimg)
odfs = img2.get_data()


odfs.clip(min=0)  #in attempt tp remove Gibbs ringingn

#if not args.fa_thr:
#	fa_thr=0.4
#else:
#	fa_thr=args.fa_thr

#if not args.seeds_count:
#	seeds_count=2
#else:
#	seeds_count=args.seeds_count
#	seeds_count=seeds_count.astype('int')
	
#if not args.max_angle:
#	max_angle = 20
#else:
#	max_angle=args.max_angle

#if not args.step_size:
#	step_size=0.175
#else:
#	step_size=args.step_size

if args.seed_mask:
        img =(args.seed_mask)
        imgm = nib.load(img)
        seed_mask = imgm.get_data()
        seeds = random_seeds_from_mask(seed_mask, seeds_count=2)
else:
	seeds = random_seeds_from_mask(fa > fa_thr, seeds_count=seeds_count)

tissue_classifier = ThresholdTissueClassifier(fa, 0.4)
prob_dg = ProbabilisticDirectionGetter.from_pmf(odfs, max_angle=20, sphere=sphere,
						pmf_threshold=0.7)
all_streamlines_tissue_classifier = LocalTracking(prob_dg,
                                               tissue_classifier,
                                               seeds, 
                                               affine=np.eye(4), step_size=0.175,
                                               return_all=True)
streamlines = Streamlines(all_streamlines_tissue_classifier)
save_trk(args.output_trk,
                     streamlines,
                     img2.affine,
                     shape=img2.shape[:3], vox_size=img2.header.get_zooms()[:3])


