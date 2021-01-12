#!/usr/bin/env python
"""Read netCDF files to compare dark levels with and without the spinning telescope."""

__version__ = "0.0.1"
__author__ = "Christopher Field <Christopher.T.Field@nasa.gov>"

import argparse
from copy import deepcopy
import datetime
from math import ceil, floor
from pprint import pprint
import os
from pathlib import Path
import sys

from matplotlib import pyplot as plt
import matplotlib as mpl                # Used to get color map.
from mpl_toolkits.mplot3d import axes3d
from netCDF4 import Dataset         #https://unidata.github.io/netcdf4-python/netCDF4/index.html
import numpy as np
import numpy.ma as ma
from pptx import Presentation
from pptx.util import Inches
from scipy.fftpack import fft
from scipy.optimize import curve_fit

# sys.path.append(os.path.abspath(os.curdir))
# import .context

from python_pptx_helper.create_ppt import create_ppt
import read_oci_1a_class as oci_ss
# from read_oci_1a_class import READ_OCI_SPA_SPE_L1A as oci_ss


def commandline():
    """Define the command line arguments. Return parser.

    Usage: args = commandline().parse_args()

    This structure permits script documentation by sphnix. Use  .. argparse::.
    See https://readthedocs.org/projects/sphinx-argparse/ for more info.
    """

    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=" "
    )
    parser.add_argument('data_directory',
                        help='Path to source data folder.')
    parser.add_argument('ppt_input',
                        help='Path to source Power Point Template file')
    parser.add_argument('ppt_output',
                        help='Path to destination Power Point Template file')
    # parser.add_argument('ppt_output',
    #                     help='Path to destination Power Point Template file')

    return parser


def plot(x, rotating, static, figsize, gridspec_kw):

    stat_rot = np.zeros((rotating.shape[0], 2), dtype=np.float64)
    stat_rot[:, 0] = rotating.mean(axis=1)
    stat_rot[:, 1] = rotating.std(axis=1)

    stat_static = np.zeros((static.shape[0], 2), dtype=np.float64)
    stat_static[:, 0] = static.mean(axis=1)
    stat_static[:, 1] = static.std(axis=1)
    fig, ax = plt.subplots(2, 2, sharex=True, figsize=figsize, gridspec_kw=gridspec_kw)
    fig.suptitle("Comparison of static and rotating dark counts from XINA markers 4008 and 4009 (DN16)")
    ax[0, 0].set_title("Mean counts across all scans and spatial pixels")
    ax[0, 1].set_title("std of counts across all scans and spatial pixels")
    ax[0, 0].set_ylabel("Dark counts")
    ax[1, 0].set_ylabel("Differnce in dark counts")
    [a.set_xlabel("Wavelength (nm)") for a in ax[1, :]]
    ax[0, 0].plot(x, stat_rot[:, 0], label="Rotating")
    ax[0, 0].plot(x, stat_static[:, 0], label="Static")
    ax[1, 0].plot(x, stat_rot[:, 0] - stat_static[:, 0], label="Rotating - Static")
    ax[0, 1].plot(x, stat_rot[:, 1], label="Rotating")
    ax[0, 1].plot(x, stat_static[:, 1], label="Static")
    ax[1, 1].plot(x, stat_rot[:, 1] - stat_static[:, 1], label="Rotating - Static")
    return fig, ax


def main(args):
    """Main routine. Could be imported and called by another Python program

    :param args:    Argparse generated command line arguments
    :return:        None
    """
    pattern = ""

    ppx_input = "/Users/cfield/Documents/PACE/python_code/Presentation_template_in2.pptx"
    ppx_output = "/Users/cfield/Documents/PACE/micro.pptx"
    # ppx_output = None
    # ppx_input = ppx_output
    midnight = datetime.datetime(2020, 11, 9, 0, 0, 0)

    root = Path("/Users/cfield/Documents/PACE/Data/microphonics_test/")
    print("Rotating Files")
    data_path_rot = root / "rotating/"
    start_rot = datetime.datetime(2020, 11, 9, 22, 17, 0)
    end_rot = datetime.datetime(2020, 11, 9, 22, 18, 26)
    files_rot = oci_ss.get_pace_oci_l1a_files(data_path_rot, pattern)
    pprint(files_rot)

    print("Static Files")
    data_path_static = root / "static/"
    start_static = datetime.datetime(2020, 11, 9, 22, 21, 37)
    end_static = datetime.datetime(2020, 11, 9, 22, 23, 11)
    files_static = oci_ss.get_pace_oci_l1a_files(data_path_static, pattern)
    pprint(files_static)
    # ppx_input = ppx_output  # Append new slides to old file

    # ppx_output = None
    pattern = ""
    dark_zero_limit = 40
    zero_limit = 41

    figsize = (11.5, 5.75)
    gridspec_kw = {"left": 0.07,
                   "right": 0.95,
                   "top": 0.90,
                   "bottom": 0.15,
                   "wspace": 0.10,
                   "hspace": 0.20
                   }

    oci_rot = oci_ss.READ_OCI_SPA_SPE_L1A(None)
    for file in files_rot:
        oci_rot.append(oci_ss.READ_OCI_SPA_SPE_L1A(file))
    oci_rot.append_conclude()

    oci_static = oci_ss.READ_OCI_SPA_SPE_L1A(None)
    for file in files_static:
        oci_static.append(oci_ss.READ_OCI_SPA_SPE_L1A(file))
    oci_static.append_conclude()

    DC_red_data_flat_rot = np.swapaxes(oci_rot.DC_red_data, 0, 1).reshape(oci_rot.DC_red_data.shape[1], -1)
    DC_red_data_flat_static = np.swapaxes(oci_static.DC_red_data, 0, 1).reshape(oci_static.DC_red_data.shape[1], -1)
    fig_red, ax_red = plot(oci_rot.red_wavelengths, DC_red_data_flat_rot, DC_red_data_flat_static, figsize, gridspec_kw)
    ax_red[0, 0].legend(loc=(0.5, 0.7))
    ax_red[0, 1].legend(loc=(0.5, 0.7))
    ax_red[1, 0].legend(loc='upper center')
    ax_red[1, 1].legend(loc='upper center')

    fig_filename_red = "/Users/cfield/Desktop/temp/microphonics_red.png"
    fig_red.savefig(fig_filename_red, dpi=600)

    DC_swir_data_flat_rot = np.swapaxes(oci_rot.DC_swir_data, 0, 1).reshape(oci_rot.DC_swir_data.shape[1], -1)
    DC_swir_data_flat_static = np.swapaxes(oci_static.DC_swir_data, 0, 1).reshape(oci_static.DC_swir_data.shape[1], -1)
    fig_swir, ax_swir = plot(range(5), DC_swir_data_flat_rot[:5, :], DC_swir_data_flat_static[:5, :], figsize, gridspec_kw)
    ax_swir[0, 0].legend(loc='upper right')
    ax_swir[0, 1].legend(loc='upper left')
    ax_swir[1, 0].legend(loc='upper left')
    ax_swir[1, 1].legend(loc='upper left')
    ax_swir[1, 0].set_ylim(-0.15, 0.10)
    ax_swir[1, 0].set_xticks(range(5))
    ax_swir[1, 0].set_xticklabels(["940", "1038", "1250", "1250 HG", "1378"])
    fig_filename_swir = "/Users/cfield/Desktop/temp/microphonics_swir.png"
    fig_swir.savefig(fig_filename_swir, dpi=600)

    # Now put the plot in a Power Point file
    print("Processing is all done. Now make the Power Point file.")
    print("Defining the slides")
    slides = [{"layout" : 0,    # Title slide
               "text" : [(0, "Microphonics Investigation"),
                         (1, "Generated on {:%m-%d-%Y}\n by Christopher Field".format(datetime.date.today()))]
               },
              {"layout" : 11,   # Picture slide
               "text" : [(0, "Red Dark Level (DN16)"),
                         # (10, "Date"),
                         # (11, "Footer"),
                         # (12, "Slide #")
                         (13, "Source data from {:}".format(root)),
                        ],
               "picture" : [(1, fig_filename_red)]
               },
              {"layout" : 11,   # Picture slide
               "text" : [(0, "SWIR Dark Level"),
                         # (10, "Date"),
                         # (11, "Footer"),
                         # (12, "Slide #")
                         (13, "Source data from {:}".format(root)),
                        ],
               "picture" : [(1, fig_filename_swir)]
               },
              ]
    print("Creating the file.")
    create_ppt(ppx_input, ppx_output, slides, show_placeholders=False)

    plt.show()


if __name__ == "__main__":
    # args = commandline().parse_args()
    args = None
    print(args)
    main(args)
