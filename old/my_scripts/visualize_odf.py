#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import argparse
import logging
import nibabel

from scilpy.io.utils import assert_inputs_exist

from dipy.data import get_sphere
from dipy.reconst.shm import sh_to_sf, order_from_ncoef
from dipy.viz import window, actor

DESCRIPTION = """
    Script to visualize SH with dipy.
    """

def buildArgsParser():
    p = argparse.ArgumentParser(description=DESCRIPTION,
                                formatter_class=argparse.RawTextHelpFormatter)

    p.add_argument(
        'sh_file', type=str,
        help='SH file (odf, fodf)')

    group = p.add_mutually_exclusive_group()
    group.add_argument('-x', '--x', action='store_true',
                       help='display yz slice (sagital)')
    group.add_argument('-y', '--y', action='store_true',
                       help='display xz slice (coronal) [default]')
    group.add_argument('-z', '--z', action='store_true',
                       help='display xy slice (axial)')

    p.add_argument(
        '--index', '-i' , type=int,
        help='index of the slice')

    p.add_argument(
        '--basis', action='store',
        default='dipy', type=str, choices=["mrtrix", "dipy"],
        help="Basis used for the spherical harmonic coefficients.\n" +
        "(must be 'mrtrix' or 'dipy'). [%(default)s]")

    p.add_argument(
        '--sphere', action='store', default='repulsion724', type=str,
        help="Sphere used for discretization [%(default)s]")

    p.add_argument(
        '--mask', action='store',
        help="voxel mask")

    p.add_argument(
        '--normalize', action='store_true',
        help="normalize fodf size by diffusion mask")

    return p


def main():
    parser = buildArgsParser()
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    assert_inputs_exist(parser, required=[args.sh_file], optional=[args.mask])

    # Load data
    sh_data = nibabel.load(args.sh_file).get_data().astype('float64')
    if args.mask:
        mask = nibabel.load(args.sh_file).get_data()
        sh_data *= mask

    # Index for the chosen slice (or mid section)
    if args.index:
        index = args.index
    else:
        if args.x:
            index  = sh_data.shape[0] // 2
        elif args.z:
            index  = sh_data.shape[2] // 2
        else:
            index  = sh_data.shape[1] // 2

    # Slicing of the data
    if args.x:
        sh_data_cut = sh_data[index:index+1, :, :]
    elif args.z:
        sh_data_cut = sh_data[:, :, index:index+1]
    else:
        sh_data_cut = sh_data[:, index:index+1, :]

    # Basis and SH coefficient to SF
    if args.basis == 'dipy':
        basis_name = None
    else:
        basis_name = args.basis

    sphere = get_sphere(args.sphere)
    sh_order = order_from_ncoef(sh_data.shape[-1])
    sf_data = sh_to_sf(sh_data_cut, sphere, sh_order, basis_type=basis_name)

    model_kernel = actor.odf_slicer(
        sf_data, sphere=sphere, norm=args.normalize, colormap='jet', scale=0.4)

    # Rendering
    ren = window.Renderer()

    if args.x:
        model_kernel.display(x=0)
        camera_position = (10, 0, 0)
    elif args.z:
        model_kernel.display(z=0)
        camera_position = (0, 0, 10)
    else:
        model_kernel.display(y=0)
        camera_position = (0, 10, 0)

    ren.set_camera(position=camera_position, view_up=(0, 0, 1))
    ren.add(model_kernel)
    window.show(ren,size=(1024, 800))


if __name__ == "__main__":
    main()
