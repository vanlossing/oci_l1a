#!/usr/bin/env python
"""test read_oci_1a_class"""

__version__ = "0.0.1"
__author__ = "Christopher Field <Christopher.T.Field@nasa.gov>"

import argparse
from copy import deepcopy
from pprint import pprint
import os
from pathlib import Path
import sys

from matplotlib import pyplot as plt
import matplotlib as mpl                # Used to get color map.
from netCDF4 import Dataset         #https://unidata.github.io/netcdf4-python/netCDF4/index.html
import numpy as np
import numpy.ma as ma

# sys.path.append(os.path.abspath(os.curdir))
# import .context

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

    data_path = Path("/Users/cfield/Documents/PACE/python_code/oci_l1a/test/dark_view_correction_examples/")
    if 1:
        pattern = "2020055138.20200224T183300"
        # Read Nick's background subtracted data.
        pattern_reference = "/ncollins_PACE_OCI_2020055138.20200224T183300.L1A.nc"
        wavelength_file = "red_waves_8x8.csv"                   # Used to verify internal wavelength calculation.
    else:
        pattern = "2020055140.20200224T184932"
        # Read Nick's background subtracted data.
        pattern_reference = "/ncollins_PACE_OCI_DIAG_2020055140.20200224T184932.L1A.nc"
        wavelength_file = "red_waves_1x1.csv"                  # Used to verify internal wavelength calculation.

    files = oci_ss.get_pace_oci_l1a_files(data_path, pattern)
    # print("os.curdir()")
    # print(os.path.abspath(os.curdir))
    # pprint(files)

    # oci = oci_ss.read_oci_1a_class(None)
    # for file in path:
    #     oci.append(oci_ss.read_oci_1a_class(file))
    # oci.append_conclude()

    oci__ = oci_ss.READ_OCI_SPA_SPE_L1A(files[0])
    oci = deepcopy(oci__)
    # oci.subtract_dark()
    # oci2 = deepcopy(oci)


    Sci_red_data_global_mean = oci.Sci_red_data[1:-1, :, :].mean()
    Sci_red_data_global_std = oci.Sci_red_data[1:-1, :, :].std()
    # print("Sci_red_data")
    # print("Pre-mask:   Sci_red_data_global_mean,  Sci_red_data_global_std", Sci_red_data_global_mean, Sci_red_data_global_std)
    # print("Pre-mask:   Sci_red_data_global_max(), Sci_red_data_global_min()", oci.Sci_red_data[1:-1, :, :].max(), oci.Sci_red_data[1:-1, :, :].min())
    ma.masked_outside(oci.Sci_red_data,
                      Sci_red_data_global_mean - 4*Sci_red_data_global_std,
                      Sci_red_data_global_mean + 4*Sci_red_data_global_std,
                      copy=False)
    Sci_red_data_global_mean = oci.Sci_red_data[1:-1, :, :].mean()
    Sci_red_data_global_std = oci.Sci_red_data[1:-1, :, :].std()
    # print("Post-mask:  Sci_red_data_global_mean,  Sci_red_data_global_std", Sci_red_data_global_mean, Sci_red_data_global_std)
    # print("Post-mask:  Sci_red_data_global_max(), Sci_red_data_global_min()", oci.Sci_red_data[1:-1, :, :].max(), oci.Sci_red_data[1:-1, :, :].min())
    # print(oci.Sci_red_data.flags)

    fig, ax = plt.subplots(2, 2, sharex=True)
    fig.suptitle("Red CCD counts vs time from {:}.".format(oci.path.name))
    ax[0, 0].plot(oci.scan_start_time, oci.Sci_red_data[:, :, 0], 'x')
    ax[0, 0].set_title("Sci_red[:, :, 0] vs time")
    ax[1, 0].plot(oci.scan_start_time, oci.Sci_red_data[:, :, 400], 'x')
    ax[1, 0].set_title("Sci_red[:, :, 400] vs time")
    ax[0, 1].plot(oci.scan_start_time, oci.Sci_red_data[:, 0, :], 'x')
    ax[0, 1].set_title("Sci_red[:, 0, :] vs time")
    ax[1, 1].plot(oci.scan_start_time, oci.Sci_red_data[:, 30, :], 'x')
    ax[1, 1].set_title("Sci_red[:, 400, :] vs time")


    Sci_red_data_mean = oci.Sci_red_data[1:-1, :, :].mean(axis=0)
    Sci_red_data_std = oci.Sci_red_data[1:-1, :, :].std(axis=0)

    fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
    fig.suptitle("Mean Red CCD counts vs $\lambda$ and Angle \n{:}.".format(oci.path.name))
    ax.set_xlabel('Wavelength (ns)')
    ax.set_ylabel('Scan Angle (deg)')
    ax.set_zlabel('File CCD Counts')
    # angle, lam = np.meshgrid(range(Sci_red_data_mean.shape[1]), range(Sci_red_data_mean.shape[0]))
    oci.pixel_angle[oci.pixel_angle > 200] = np.nan
    angle, lam = np.meshgrid(oci.pixel_angle, oci.red_wavelengths)
    surf = ax.plot_surface(lam, angle, ma.filled(Sci_red_data_mean, np.nan),
                           cmap='viridis', edgecolor='none',
                           vmin=0, vmax=1000) #cmap=mpl.rcParams['image.cmap'])
    fig.colorbar(surf, ax=ax,
                 shrink=0.5, aspect=5)

    oci.subtract_dark(start=0)

    Sci_red_data_global_mean = oci.Sci_red_data[1:-1, :, :].mean()
    Sci_red_data_global_std = oci.Sci_red_data[1:-1, :, :].std()
    # print("Pre-mask:   Sci_red_data_global_mean,  Sci_red_data_global_std", Sci_red_data_global_mean, Sci_red_data_global_std)
    # print("Pre-mask:   Sci_red_data_global_max(), Sci_red_data_global_min()", oci.Sci_red_data[1:-1, :, :].max(), oci.Sci_red_data[1:-1, :, :].min())
    ma.masked_outside(oci.Sci_red_data,
                      Sci_red_data_global_mean - 1 * Sci_red_data_global_std,
                      Sci_red_data_global_mean + 1 * Sci_red_data_global_std,
                      copy=False)
    Sci_red_data_global_mean = oci.Sci_red_data[1:-1, :, :].mean()
    Sci_red_data_global_std = oci.Sci_red_data[1:-1, :, :].std()
    # print("Post-mask:  Sci_red_data_global_mean,  Sci_red_data_global_std", Sci_red_data_global_mean, Sci_red_data_global_std)
    # print("Pre-mask:   Sci_red_data_global_max(), Sci_red_data_global_min()", oci.Sci_red_data[1:-1, :, :].max(), oci.Sci_red_data[1:-1, :, :].min())
    # print(oci.Sci_red_data.flags)

    Sci_red_data_mean = oci.Sci_red_data.mean(axis=0)
    Sci_red_data_std = oci.Sci_red_data.std(axis=0)
    fig2, ax2 = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
    fig2.suptitle("Mean Red CCD counts vs $\lambda$ and Angle \n{:}.".format(oci.path.name))
    ax2.set_xlabel('Wavelength (ns)')
    ax2.set_ylabel('Scan Angle (deg)')
    ax2.set_zlabel('Background subtracted Counts')
    angle, lam = np.meshgrid(oci.pixel_angle, oci.red_wavelengths)
    surf2 = ax2.plot_surface(lam, angle, ma.filled(Sci_red_data_mean, np.nan), cmap='viridis', edgecolor='none',
                             vmin=-1, vmax=1) #cmap=mpl.rcParams['image.cmap'])
    fig2.colorbar(surf2, ax=ax2,
                  shrink=0.5, aspect=5)


    files = oci_ss.get_pace_oci_l1a_files(data_path, pattern_reference)
    with Dataset(files[0], 'r') as fid_cdf:
        science_group = fid_cdf.groups['science_data']  # Read the CCD raw data
        Sci_red_ds = science_group.variables['sci_red_ds'][1:-1, :, :]

    # print("Sci_red_ds")
    # print(Sci_red_ds.flags)
    Sci_red_ds_global_mean = Sci_red_ds.mean()
    Sci_red_ds_global_std = Sci_red_ds.std()
    # print("Pre-mask:   Sci_red_ds_global_mean,    Sci_red_ds_global_std", Sci_red_ds_global_mean, Sci_red_ds_global_std)
    # print("Pre-mask:   Sci_red_ds.max(),          Sci_red_ds.min()     ", Sci_red_ds[1:-1, :, :].max(), Sci_red_ds[1:-1, :, :].min())
    # print(type(Sci_red_ds))
    Sci_red_ds = ma.masked_outside(Sci_red_ds,
                                   0.1, #Sci_red_ds_global_mean - 1,
                                   Sci_red_ds_global_mean + 1,
                                   copy=False)
    Sci_red_ds_global_mean = Sci_red_ds.mean()
    Sci_red_ds_global_std = Sci_red_ds.std()
    # print("Post-mask:  Sci_red_ds_global_mean,    Sci_red_ds_global_std", Sci_red_ds_global_mean, Sci_red_ds_global_std)
    # print("Post-mask:  Sci_red_ds.max(),          Sci_red_ds.min()     ", Sci_red_ds[1:-1, :, :].max(), Sci_red_ds[1:-1, :, :].min())

    Sci_red_ds_mean = Sci_red_ds.mean(axis=0)
    Sci_red_ds_std = Sci_red_ds.std(axis=0)
    # print("Sci_red_ds.mean(), Sci_red_ds.std() ", Sci_red_ds.mean(), Sci_red_ds.std())
    fig3, ax3 = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
    fig3.suptitle("Mean Red CCD counts vs $\lambda$ and Angle \n{:}.".format(files[0].name))
    ax3.set_xlabel('Wavelength (ns)')
    ax3.set_ylabel('Scan Angle (deg)')
    ax3.set_zlabel("Nick's Background subtracted Counts")
    # angle, lam = np.meshgrid(range(Sci_red_data_mean.shape[1]), range(Sci_red_data_mean.shape[0]))
    # oci.pixel_angle[oci.pixel_angle > 200] = np.nan

    angle, lam = np.meshgrid(oci.pixel_angle, oci.red_wavelengths)
    surf3 = ax3.plot_surface(lam, angle, ma.filled(Sci_red_ds_mean, np.nan), cmap='viridis',
                             edgecolor='none',
                             vmin=-1, vmax=1) #cmap=mpl.rcParams['image.cmap'])
    fig3.colorbar(surf3, ax=ax3,
                  shrink=0.5, aspect=5)

    wavelengths = np.loadtxt(data_path/wavelength_file)
    print("Wavelengths from csv file")
    print(wavelengths)
    print("Wavelengths from nc file")
    print(oci.red_wavelengths)


    plt.show()
