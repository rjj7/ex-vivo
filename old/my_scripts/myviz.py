#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import logging
import nibabel
import numpy as np

from dipy.data import get_sphere
from dipy.reconst.shm import sh_to_sf, order_from_ncoef
from dipy.viz import window, actor, ui
from dipy.io.utils import is_header_compatible, get_reference_info
from dipy.tracking.streamline import transform_streamlines
from dipy.io.streamline import load_tractogram

p= argparse.ArgumentParser(
        description='Script to visualize ODFs (SH or SF), tractography (trk or tck), and peaks.')

p.add_argument(
        '--anat_file', help='NIFTI 3D volume' )
p.add_argument(
	'--track_file', help='tractography file in tck or trk formats')
p.add_argument(
	'--odf_file', help='ODFs in either SH or SF format')
p.add_argument(
        '--peaks_dirs_file', help='Peaks directions file')
p.add_argument(
        '--peaks_values_file', help='Peaks values file')

p.add_argument('--index_x', type=int, help='index of the slice')
p.add_argument('--index_y', type=int, help='index of the slice')
p.add_argument('--index_z', type=int, help='index of the slice')

p.add_argument('-x', '--x', action='store_true',
                       help='display yz slice (sagital)[default]')
p.add_argument('-y', '--y', action='store_true',
                       help='display xz slice (coronal)')
p.add_argument('-z', '--z', action='store_true',
                       help='display xy slice (axial)')
p.add_argument('--interactive', '-i', type=bool, default=True, 
		help='Activate or suppress visualization')

p.add_argument(
        '--basis', action='store',
        type=str, choices=["descoteaux07", "tournier07"],
        help="Basis used for the spherical harmonic coefficients.\n" +
        "If ODFs obtained with MRtrix choose 'tournier07'; if obtained \n" +
	"with dipy choose 'descoteaux07'. Default is Descoteux07.", default='descoteaux07')
p.add_argument('--sh_order', type=int, default=6,
		help='Order of the shperical harmonics (lmax). Default = 6')
p.add_argument(
        '--sphere', action='store', default='symmetric362', type=str,
        help="Sphere used for discretization ['%(default)s'].\n" +
	"(must be the same number of your ODFs 4th dimension (if SF)")
p.add_argument(
        '--mask', action='store', default=None,
        help="Binary mask. Only displays ODFs and peaks inside the mask")
p.add_argument(
        '--norm', type=bool, default=True,
        help="normalize fodf size by diffusion mask")
p.add_argument(
	'--peaks_color', type=int, 
	help='RGB color code format', default=(0,0,1))    
p.add_argument(
	'--odf_colormap', type=str, default='jet',
	 help="Color map to use when displaying ODFs.\n" +
       " If None then white color is used. Matplotlib colormaps are\n" +
       " supported (e.g., 'inferno'). Default is ['%(default)s']")
p.add_argument('--scale_odfs', type=float, default=0.6, help='Distance between spheres')
p.add_argument('--odfs_opacity', type=float, default=1.0, help='ODFs Opacity')
p.add_argument('--anat_opacity', type=float, default=1.0, help='Slice Opacity')
p.add_argument('--peaks_opacity', type=float, default=1.0, help='Peaks Opacity')
p.add_argument('--background', type=float, help='Color background. RGB code (e.g. (0.5, 0.5, 0.5). Default is black', default=(0,0,0))
p.add_argument('--linewidth', type=float, default=1.0, help='Peaks linewidth')
p.add_argument('--zoom', type=float, help='Camera zoom', default=1)
p.add_argument('--cam_position', type=str, help='Choice of axis for camera position.', default='x' )
p.add_argument('--window_size', type=int, help='Renderer window size', default=(1024,800))
p.add_argument('--screenshot', '-s', help='Path to .png file')
p.add_argument('--nd', help='No interactive')
args = p.parse_args()

if not args.anat_file:
	if not args.odf_file:
		if not args.track_file:
			p.error("You must specify at least one file input.\n" +
                		" Anat, track, ODFs, peaks, or dyads")

ren = window.Renderer()
if args.background:
	ren.background(args.background)

sphere = get_sphere(args.sphere)

# Load data
if args.odf_file:
	print("Loading ODF file..")
	img = nibabel.load(args.odf_file)
	odfs = img.get_data().astype('float64')
	if odfs.shape[-1] > 45:
		odfs_format = 'sf'
		print ('ODFs in is SF format')
	else: 
		odfs_format = 'sh'
		print ('ODFs in is SH format.')
	shape = odfs.shape
	print(shape)

if args.anat_file:
	img = nibabel.load(args.anat_file)
	anat = img.get_data()
	shape = anat.shape
	
#if args.reference:
#	img = nibabel.load(args.reference)

if args.mask:
	mask = nibabel.load(args.mask).get_data()
	odfs *= mask

if args.track_file:
	if str(args.track_file).endswith('trk'):
		track = load_tractogram(args.track_file, reference = 'same')
	elif str(args.track_file).endswith('tck'):
		if not (args.anat_file or args.odf_file or args.reference):
			p.error("If track format tck you need to\n" +
				"provide reference image")
		track = load_tractogram(args.track_file, img)
	else:
		 p.error("Track_file format not unrecognized")	
	if (args.anat_file or args.odf_file):
		if not is_header_compatible(img, track):
			p.error("Reference file and track_file\n" +
                        "do not have compatible headers!")
		track = transform_streamlines(track, 
						np.linalg.inv(img.affine))
	else:
		track = track.streamlines
	track_actor = actor.line(track)
	ren.add(track)

if args.peaks_dirs_file:
	peaks_dirs = nibabel.load(args.peaks_dirs_file).get_data()
	#if len(peaks_dirs) < 5:
		#if args.peaks_oct:
		#	peaks_dirs = np.expand_dims(peaks_dirs, axis = 3)
		#else:
		#	p.error("Peaks file needs to be in dipy format") 
		#To do: add reshape peaks if not

if args.peaks_values_file:
	peaks_values = nibabel.load(args.peaks_values_file).get_data() 
else:
	peaks_values = None
      
if args.x:
	if args.index_x:
		ind_x  = args.index_x
	else:
		ind_x = shape[0] // 2
	print(ind_x)
	if args.anat_file:
		x_actor = actor.slicer(anat, opacity=args.anat_opacity)
		x_actor.display_extent(ind_x, ind_x,
					0, shape[1] -1,
					0, shape[2] -1)
		ren.add(x_actor)
	if args.odf_file:
		odfs_cut = odfs[ind_x:ind_x+1, :, :]
		if odfs_format == 'sh':
			odfs_cut = sh_to_sf(odfs_cut, sphere, 
				args.sh_order, basis_type = args.basis)
		x_odf_actor = actor.odf_slicer(odfs_cut, sphere=sphere,
			colormap=args.odf_colormap, scale=args.scale_odfs,
			opacity=args.odfs_opacity)
		x_odf_actor.display(x=0)
		ren.add(x_odf_actor)
	if args.peaks_dirs_file:
		x_peaks_actor = actor.peak_slicer(peaks_dirs, peaks_values,
			linewidth=args.linewidth, colors= args.peaks_color,
			opacity=args.peaks_opacity)
		x_peaks_actor.display_extent(ind_x, ind_x,
                                        0, shape[1] -1,
                                        0, shape[2] -1)
		ren.add(x_peaks_actor)
		
if args.y:
	if args.index_y:
		ind_y = args.index_y
	else:
		ind_y = int(np.round(shape[1] // 2))
	if args.anat_file:
                y_actor = actor.slicer(anat, opacity=args.anat_opacity)
                y_actor.display_extent(0, shape[0] - 1,
					ind_y, ind_y,
                                        0, shape[2] -1)
                ren.add(y_actor)
	if args.odf_file:
		odfs_cut = odfs[:, ind_y:ind_y+1, :]
		if odfs_format == 'sh':
                        odfs_cut = sh_to_sf(odfs_cut, sphere,
                                args.sh_order, basis_type = args.basis)
		y_odf_actor = actor.odf_slicer(odfs_cut, sphere=sphere,
                        norm=args.norm, colormap=args.odf_colormap,
                        scale=args.scale_odfs, opacity=args.odfs_opacity)
		y_odf_actor.display(y=0)
		ren.add(y_odf_actor)
	if args.peaks_dirs_file:
                y_peaks_actor = actor.peak_slicer(peaks_dirs, peaks_values,
                        linewidth=args.linewidth, colors= args.peaks_color,
                        opacity=args.peaks_opacity)
                y_peaks_actor.display_extent(0, shape[0] - 1,
					ind_y, ind_y,
                                        0, shape[2] -1)
                ren.add(y_peaks_actor)

if args.z:
	if args.index_z:
		ind_z = args.index_z
	else:
		ind_z = int(np.round(shape[2] // 2))
	if args.ant_file:
                z_actor = actor.slicer(anat, opacity=args.anat_opacity)
                z_actor.display_extent(0, shape[0] - 1,
					0, shape[1] - 1,
                                        ind_z, ind_z)
                ren.add(z_actor)
	if args.odf_file:
		odfs_cut = odfs[:, :, ind_x:ind_x+1]
		if odfs_format == 'sh':
                        odfs_cut = sh_to_sf(odfs_cut, sphere,
                                args.sh_order, basis_type = args.basis)
		z_odf_actor = actor.odf_slicer(odfs_cut, sphere=sphere,
                        norm=args.norm, colormap=args.odf_colormap,
                        scale=args.scale_odfs, opacity=args.odfs_opacity)
		z_odf_actor.display(z=0)
		ren.add(z_odf_actor)
	if args.peaks_dirs_file:
		z_peaks_actor = actor.peak_slicer(peaks_dirs, peaks_values,
                        linewidth=args.linewidth, colors= args.peaks_color,
                        opacity=args.peaks_opacity)
		z_peaks_actor.display_extent(0, shape[0] - 1,
					0, shape[1] - 1,
                                        ind_z, ind_z)
		ren.add(z_peaks_actor)

if not args.x:
	if not args.y:
		if not args.z:
			ind_x = int(np.round(shape[0] // 2))
			print(ind_x)
			if args.anat_file:
				x_actor = actor.slicer(anat, opacity=args.anat_opacity)
				x_actor.display_extent(index_x, index_x,
                                        		0, shape[1] -1,
                                        		0, shape[2] -1)
				ren.add(x_actor)
			if args.odf_file:
				odfs_cut = odfs[ind_x:ind_x+1, :, :]
				if odfs_format == 'sh':
					odfs_cut = sh_to_sf(odfs_cut, sphere,
                                		args.sh_order, basis_type = args.basis)
				x_odf_actor = actor.odf_slicer(odfs_cut, sphere=sphere,
                        		norm=args.norm, colormap=args.odf_colormap,
                        		mask = mask, scale=args.scale_odfs,
                        		opacity=args.odfs_opacity)
				x_odf_actor.display(x=0)
				ren.add(x_odf_actor)
			if args.peaks_dirs_file:
				peaks_actor = actor.peak_slicer(peaks_dirs, peaks_values,
                        		linewidth=args.linewidth, colors= args.peaks_color,
                        		opacity=peaks_opacity)
				peaks_actor.display_extent(index_x, index_x,
                                        		0, shape[1] -1,
                                        		0, shape[2] -1)

# Rendering
if args.cam_position == 'x':
        camera_position = (10, 0, 0)
elif args.cam_position == 'y':
        camera_position = (0, 0, 10)
elif args.cam_position == 'z':
	camera_position = (0, 10, 0)

ren.set_camera(position=camera_position, view_up=(0, 0, 1))
ren.zoom(args.zoom)

if not args.nd:
	window.show(ren, size=(args.window_size))
if args.screenshot:
	window.record(ren, out_path=args.screenshot, size=(args.window_size))

