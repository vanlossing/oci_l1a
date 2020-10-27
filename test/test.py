#!/usr/bin/env python
"""test read_oci_1a_class"""

__version__ = "0.0.1"
__author__ = "Christopher Field <Christopher.T.Field@nasa.gov>"

import argparse
import os
import sys

import numpy as np
from matplotlib import pyplot as plt

# sys.path.append(os.path.abspath(os.curdir))
# import .context

from read_oci_1a_class import READ_OCI_SPA_SPE_L1A as OCI_SPA_SPE

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
    parser.add_argument('input',
                        help='Input argument')

    return parser


def main(args):
    """Main routine. Could be imported and called by another Python program

    :param args:    Argparse generated command line arguments
    :return:        None
    """
    pass


if __name__ == "__main__":

    print("os.curdir()")
    print(os.path.abspath(os.curdir))
    path = ['/Users/cfield/Downloads/PACE_OCI_SPEC_2020076013.20200316T132700.L1A.nc',
            '/Users/cfield/Downloads/PACE_OCI_SPEC_2020076013.20200316T132800.L1A.nc',
            '/Users/cfield/Downloads/PACE_OCI_SPEC_2020076013.20200316T133100.L1A.nc']

    oci = OCI_SPA_SPE(None)
    for file in path:
        oci.append(OCI_SPA_SPE(file))
    oci.append_conclude()

    fig, ax = plt.subplots(2, 2, sharex=True)
    fig.suptitle("Red CCD counts vs time from all files.")
    ax[0, 0].plot(oci.scan_start_time, oci.Sci_red_data[:, :, 0], 'x')
    ax[0, 0].set_title("Sci_red[:, :, 0] vs time")
    ax[1, 0].plot(oci.scan_start_time, oci.Sci_red_data[:, :, 400], 'x')
    ax[1, 0].set_title("Sci_red[:, :, 400] vs time")
    ax[0, 1].plot(oci.scan_start_time, oci.Sci_red_data[:, 0, :], 'x')
    ax[0, 1].set_title("Sci_red[:, 0, :] vs time")
    ax[1, 1].plot(oci.scan_start_time, oci.Sci_red_data[:, 400, :], 'x')
    ax[1, 1].set_title("Sci_red[:, 400, :] vs time")

    Sci_red_data_mean = oci.Sci_red_data.mean(axis=0)
    Sci_red_data_std = oci.Sci_red_data.std(axis=0)
    # print(oci.Sci_red_data.shape, oci.Sci_red_data_mean.shape)

    fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Scan Angle')
    ax.set_zlabel('CCD Counts')
    # angle, lam = np.meshgrid(range(Sci_red_data_mean.shape[1]), range(Sci_red_data_mean.shape[0]))
    oci.pixel_angle[oci.pixel_angle > 200] = np.nan
    angle, lam = np.meshgrid(oci.pixel_angle, oci.wavelength)
    ax.plot_surface(lam, angle, ma.filled(Sci_red_data_mean, np.nan))
    plt.show()
