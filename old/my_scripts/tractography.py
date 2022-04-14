#!/usr/bin/env python
import argparse

import numpy as np
import nibabel as nib
import os
import glob

from nibabel.streamlines import save as save_trk
from nibabel.streamlines import Tractogram

from dipy.data import get_sphere
from dipy.direction import peaks_from_model, ProbabilisticDirectionGetter, DeterministicMaximumDirectionGetter
from dipy.direction.peaks import peak_directions
from dipy.io.stateful_tractogram import Space, StatefulTractogram
from dipy.tracking.streamline import Streamlines
from dipy.tracking.utils import random_seeds_from_mask
from dipy.tracking.local_tracking import LocalTracking
from dipy.io.streamline import save_trk, save_tck
from dipy.tracking.stopping_criterion import (ActStoppingCriterion, BinaryStoppingCriterion, ThresholdStoppingCriterion)


parser = argparse.ArgumentParser(
        description='Computes probabilistic or deterministic local tractography using DIPY modules from either SH or SF ODFS files')

parser.add_argument('ODFs_file', help='ODFs')
parser.add_argument('mask', help='Binary mask. Streamlines outside this mask will be truncated')
parser.add_argument('seed_mask', help='Binary mask to seed from. Seeds will be randomly placed in each voxel.')
parser.add_argument('output_file', type=str, help='Output track file, in tck or trk format')

parser.add_argument('--stopping_criteria', default='binary', choices=['binary', 'threshold', 'act'], help='Specify classifier to be used to constrain tractography. Can be binary, threshold, or ACT. Threshold classifier requires FA map, ACT requires WM GM CSF segmentation maps from T1.')
parser.add_argument('--fa', help='FA map to be provided with threshold stopping criteria')
parser.add_argument('--fa_thr', default=0.2, type=float, help='FA threshold to be used with threshold classifier')
parser.add_argument('--act_wm', help='White matter PVE segmentation map to be provided with ACT classifier')
parser.add_argument('--act_gm', help='Gray matter PVE segmentation map to be provided with ACT classifier')
parser.add_argument('--act_csf', help='CSF PVE segmentation map to be provided with ACT classifier')

parser.add_argument('--track_algo', default='prob', choices=['det', 'prob'], help='Deterministic or probabilistic tractography. Deterministic: the maxima of the spherical function (SF) the most closely aligned to the previous direction. Probabilistic: a direction drawn from the empirical distribution function defined from the SF. Default is probabilistic')

parser.add_argument('--seeds_n', dest='seeds_n', default=2, type=int,  help='Number of seeds per voxel. Randomly placed')
parser.add_argument('--step_size', dest='step_size', type=float, help='Tracking step-size in mm. Default is half of a voxel')
parser.add_argument('--pmf_treshold', dest='pmf_thr', type=float, help='', default=0.8)
parser.add_argument('--max_angle', dest='max_angle', type=int, help='', default=30)
parser.add_argument('--sphere', dest='sphere', help='Set of directions on the sphere to be used for tractography', default='symmetric362')

args = parser.parse_args()

print ('Loading data...seat back and relax!')
odfs = nib.load(args.ODFs_file).get_data()
mask_img = nib.load(args.seed_mask) 
seed_mask = mask_img.get_data().astype('bool')

sphere = get_sphere(args.sphere)

if not nib.streamlines.is_supported(args.output_file):
        parser.error('Invalid output streamline file format (must be trk or ' +
                     'tck): {0}'.format(args.output_file))

if args.step_size:
	step_size = args.step_size
else: 
	voxel_size = mask_img.header.get_zooms()[0]
	step_size = voxel_size / 2

print ('Voxel size = '+ str(voxel_size))
print ('Step size = ' + str(step_size))
	 

if args.stopping_criteria in ['binary', 'threshold', 'act']:
	if args.stopping_criteria == 'threshold':
		if not args.fa:
			parser.error('Threshold classifier requires FA map')
		fa = nib.load(args.fa).get_data()
		s_criterion = ThresholdStoppingCriterion(fa, args.fa_thr)
	if args.stopping_criteria == 'act':
		if not args.act_wm:
			parser.error('ACT classifier requires WM PVE')
		if not args.act_gm:
			parser.error('ACT classifier requires GM PVE')
		if not args.act_csf:
			parser.error('ACT classifier requires CSF PVE')
		background = np.ones(img_pve_gm.shape)
		background[nib.load(args.gm).get_data() +
            		nib.load(args.wm).get_data() +
            		nib.load(args.csf).get_data() > 0] = 0
		include_map = nib.load(args.gm).get_data()
		include_map[background > 0] = 1
		exclude_map = img_pve_csf.get_data()
		s_criterion = ActStoppingCriterion(include_map, exclude_map)
	if args.stopping_criteria == 'binary':
		mask = nib.load(args.mask).get_data()
		s_criterion = BinaryStoppingCriterion(mask == 1) 	

print ('Stopping critearia is ' + str(args.stopping_criteria))

if odfs.shape[-1] > 45:
	print ('ODFs in SF format')
	pmf = odfs.clip(min=0)
	if args.track_algo == 'Deterministic':
		dg =DeterministicMaximumDirectionGetter.from_pmf(odfs, max_angle, sphere=sphere, pmf_threshold = args.pmf_thr) 	
	dg = ProbabilisticDirectionGetter.from_pmf(pmf, max_angle=args.max_angle, sphere=sphere, pmf_threshold=args.pmf_thr) 
else:
	print ('ODFs in SH format')
	if args.track_algo == 'Deterministic':
		dg = DeterministicMaximumDirectionGetter.from_shcoeff(odfs, max_angle, sphere=sphere, pmf_threshold = args.pmf_thr)
	dg = ProbabilisticDirectionGetter.from_shcoeff(odfs, max_angle=args.max_angle, sphere=sphere, pmf_threshold=args.pmf_thr)

seeds = random_seeds_from_mask(seed_mask, affine=mask_img.affine, 
			seeds_count=args.seeds_n, seed_count_per_voxel=True)

print ('Performing...' + str(args.track_algo) + 'tractography...')
all_streamlines = LocalTracking(dg, s_criterion, seeds, 
                  affine=mask_img.affine, step_size=step_size, return_all=True)

streamlines = Streamlines(all_streamlines)

print ('Saving file...')

sft = StatefulTractogram(streamlines, mask_img, Space.RASMM) 

if str(args.output_file).endswith('trk'):
	save_trk(sft, args.output_file)
elif str(args.output_file).endswith('tck'):
	save_tck(sft, args.output_file)
else:
	parser.error('Tractography file format non available')

print ('Tractography done!')
