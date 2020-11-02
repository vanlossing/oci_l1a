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
from mpl_toolkits.mplot3d import axes3d
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
        pattern_reference_n = "/ncollins_PACE_OCI_2020055138.20200224T183300.L1A.nc"
        pattern_reference_h = "/hchoi_PACE_OCI_2020055138.20200224T183300.L1A.nc"
        dark_zero_limit = 40
        zero_limit = 41
        wavelength_file = "red_waves_8x8.csv"                   # Used to verify internal wavelength calculation.
    else:
        pattern = "2020055140.20200224T184932"
        # Read Nick's background subtracted data.
        pattern_reference_n = "/ncollins_PACE_OCI_DIAG_2020055140.20200224T184932.L1A.nc"
        pattern_reference_h = "/hchoi_PACE_OCI_DIAG_2020055140.20200224T184932.L1A.nc"
        dark_zero_limit = 328
        zero_limit = 320
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

    # Clip the far out points so that plots have reasonable scale.
    Sci_red_data_global_mean_a = oci.Sci_red_data.mean()
    Sci_red_data_global_std_a = oci.Sci_red_data.std()
    ma.masked_outside(oci.Sci_red_data,
                      Sci_red_data_global_mean_a - 4*Sci_red_data_global_std_a,
                      Sci_red_data_global_mean_a + 4*Sci_red_data_global_std_a,
                      copy=False)
    Sci_red_data_global_mean = oci.Sci_red_data[1:-1, :, :].mean()
    Sci_red_data_global_std = oci.Sci_red_data[1:-1, :, :].std()
    # print("Post-mask:  Sci_red_data_global_mean,  Sci_red_data_global_std", Sci_red_data_global_mean, Sci_red_data_global_std)
    # print("Post-mask:  Sci_red_data_global_max(), Sci_red_data_global_min()", oci.Sci_red_data[1:-1, :, :].max(), oci.Sci_red_data[1:-1, :, :].min())
    # print(oci.Sci_red_data.flags)

    fig_2d, ax_2d = plt.subplots(2, 2, sharex=True)
    fig_2d.suptitle("Red CCD counts vs time from {:}.".format(oci.path.name))
    ax_2d[0, 0].plot(oci.scan_start_time, oci.Sci_red_data[:, :, 0], 'x')
    ax_2d[0, 0].set_title("Sci_red[:, :, 0] vs time")
    ax_2d[1, 0].plot(oci.scan_start_time, oci.Sci_red_data[:, :, 400], 'x')
    ax_2d[1, 0].set_title("Sci_red[:, :, 400] vs time")
    ax_2d[0, 1].plot(oci.scan_start_time, oci.Sci_red_data[:, 0, :], 'x')
    ax_2d[0, 1].set_title("Sci_red[:, 0, :] vs time")
    ax_2d[1, 1].plot(oci.scan_start_time, oci.Sci_red_data[:, 30, :], 'x')
    ax_2d[1, 1].set_title("Sci_red[:, 400, :] vs time")

    Sci_red_data_mean = oci.Sci_red_data[1:-1, :, :].mean(axis=0)
    Sci_red_data_std = oci.Sci_red_data[1:-1, :, :].std(axis=0)

    fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
    fig.suptitle("Mean Red CCD counts vs $\lambda$ and Angle \n{:}.".format(oci.path.name))
    fig.text(0.5, 0.05,
             "Global mean {:.0f}, std {:.0f}. Clipped to 4 $\sigma$".format(Sci_red_data_global_mean_a,
                                                                            Sci_red_data_global_std_a))
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Scan Angle (deg)')
    ax.set_zlabel('File CCD Counts')
    # angle, lam = np.meshgrid(range(Sci_red_data_mean.shape[1]), range(Sci_red_data_mean.shape[0]))
    oci.pixel_angle[oci.pixel_angle > 200] = np.nan
    angle, lam = np.meshgrid(oci.pixel_angle, oci.red_wavelengths)
    surf = ax.plot_surface(lam, angle, ma.filled(Sci_red_data_mean, np.nan),
                           cmap='viridis', edgecolor='none',
                           vmin=0, vmax=1000) #cmap=mpl.rcParams['image.cmap'])

    # angle, lam = np.meshgrid(np.arange(97)*0.5, oci.red_wavelengths)
    # surf = ax.plot_surface(lam, angle, oci.DC_red_data.mean(axis=0),
    #                        cmap='jet', edgecolor='none',
    #                        vmin=0, vmax=1000) #cmap=mpl.rcParams['image.cmap'])
    #
    # angle, lam = np.meshgrid(oci.pixel_angle, oci.red_wavelengths)
    #
    fig.colorbar(surf, ax=ax,
                 shrink=0.5, aspect=5)


    # plt.show()
    # exit(1)
    oci.subtract_dark(start=dark_zero_limit)


    zlim = [-0.5, 0.5]

    # oci.Sci_red_data[:, :, :40] = 0
    oci.Sci_red_data[:, :, :zero_limit] = 0
    # oci.Sci_red_data[:, :, :41].mask = True
    # oci.Sci_red_data[:, :, :330].mask = 0


    # Sci_red_data_global_mean = oci.Sci_red_data[1:-1, :, :].mean()
    # Sci_red_data_global_std = oci.Sci_red_data[1:-1, :, :].std()
    # # print("Pre-mask:   Sci_red_data_global_mean,  Sci_red_data_global_std", Sci_red_data_global_mean, Sci_red_data_global_std)
    # # print("Pre-mask:   Sci_red_data_global_max(), Sci_red_data_global_min()", oci.Sci_red_data[1:-1, :, :].max(), oci.Sci_red_data[1:-1, :, :].min())
    # ma.masked_outside(oci.Sci_red_data,
    #                   Sci_red_data_global_mean - 1 * Sci_red_data_global_std,
    #                   Sci_red_data_global_mean + 1 * Sci_red_data_global_std,
    #                   copy=False)
    Sci_red_data_global_mean = oci.Sci_red_data[1:-1, :, :].mean()
    Sci_red_data_global_std = oci.Sci_red_data[1:-1, :, :].std()
    # print("Post-mask:  Sci_red_data_global_mean,  Sci_red_data_global_std", Sci_red_data_global_mean, Sci_red_data_global_std)
    # print("Pre-mask:   Sci_red_data_global_max(), Sci_red_data_global_min()", oci.Sci_red_data[1:-1, :, :].max(), oci.Sci_red_data[1:-1, :, :].min())
    # print(oci.Sci_red_data.flags)

    Sci_red_data_mean = oci.Sci_red_data[1:-1, :, :].mean(axis=0)
    Sci_red_data_std = oci.Sci_red_data[1:-1, :, :].std(axis=0)
    fig2, ax2 = plt.subplots(1, 2, subplot_kw={"projection": "3d"})
    fig2.suptitle("Mean Red CCD counts vs $\lambda$ and Angle \n{:}.".format(oci.path.name))
    [a.set_xlabel('Wavelength (nm)') for a in ax2]
    [a.set_ylabel('Scan Angle (deg)') for a in ax2]
    ax2[0].set_zlabel('Background subtracted Counts')
    ax2[1].set_zlabel('Background subtracted std Counts')
    ax2[0].set_title("Mean CCD Counts over {:} scans.".format(oci.Sci_red_data[1:-1, :, :].shape[0]))
    ax2[1].set_title("STD CCD Counts over {:} scans.".format(oci.Sci_red_data[1:-1, :, :].shape[0]))
    angle, lam = np.meshgrid(oci.pixel_angle, oci.red_wavelengths)
    if 1:
        # surf2 = ax2.plot_surface(lam, angle, oci.Sci_red_data[10, :, :], cmap='viridis', edgecolor='none',
        #                          vmin=-1, vmax=1) #cmap=mpl.rcParams['image.cmap'])
        surf2 = ax2[0].plot_surface(lam, angle, Sci_red_data_mean, cmap='viridis', edgecolor='none',
                                 vmin=-0.1, vmax=0.1) #cmap=mpl.rcParams['image.cmap'])
        surf2b = ax2[1].plot_surface(lam, angle, Sci_red_data_std, cmap='viridis', edgecolor='none',
                                 vmin=1.2, vmax=1.5) #cmap=mpl.rcParams['image.cmap'])
        fig2.colorbar(surf2, ax=ax2[0],
                  shrink=0.5, aspect=5)
        fig2.colorbar(surf2b, ax=ax2[1],
                  shrink=0.5, aspect=5)
    ax2[0].set_zlim(*zlim)
    ax2[1].set_zlim(1, 1.6)

    ####################################################
    # Now read Nick's dark adjusted data for comparison.
    ####################################################
    files = oci_ss.get_pace_oci_l1a_files(data_path, pattern_reference_n)
    with Dataset(files[0], 'r') as fid_cdf:
        science_group = fid_cdf.groups['science_data']  # Read the CCD raw data
        Sci_red_ds_n = science_group.variables['sci_red_ds'][:]
        Sci_red_ds_mean_n = Sci_red_ds_n.mean(axis=0)
        Sci_red_ds_std_n = Sci_red_ds_n.std(axis=0)

    # Sci_red_ds_global_mean = Sci_red_ds.mean()
    # Sci_red_ds_global_std = Sci_red_ds.std()

    # print("Sci_red_ds.mean(), Sci_red_ds.std() ", Sci_red_ds.mean(), Sci_red_ds.std())
    fig3, ax3 = plt.subplots(1, 2, subplot_kw={"projection": "3d"})
    fig3.suptitle("Mean Red CCD counts vs $\lambda$ and Angle \n{:}.".format(files[0].name))
    [a.set_xlabel('Wavelength (nm)') for a in ax3]
    [a.set_ylabel('Scan Angle (deg)') for a in ax3]
    ax3[0].set_title("Nick's Mean CCD Counts over {:} scans.".format(Sci_red_ds_n.shape[0]))
    ax3[1].set_title("Nick's STD CCD Counts over {:} scans.".format(Sci_red_ds_n.shape[0]))

    angle, lam = np.meshgrid(oci.pixel_angle, oci.red_wavelengths)
    if 1:
        # surf3 = ax3.plot_surface(lam, angle, Sci_red_ds_n[10, :, :], cmap='viridis',
        #                          edgecolor='none',
        #                          vmin=-1, vmax=1) #cmap=mpl.rcParams['image.cmap'])
    # else:
        surf3 = ax3[0].plot_surface(lam, angle, Sci_red_ds_mean_n, cmap='viridis',
                                    edgecolor='none',
                                    vmin=-0.1, vmax=0.1) #cmap=mpl.rcParams['image.cmap'])
        surf3b = ax3[1].plot_surface(lam, angle, Sci_red_ds_std_n, cmap='viridis',
                                    edgecolor='none',
                                    vmin=1.2, vmax=1.5) #cmap=mpl.rcParams['image.cmap'])
    fig3.colorbar(surf3, ax=ax3[0],
                  shrink=0.5, aspect=5)
    fig3.colorbar(surf3, ax=ax3[1],
                  shrink=0.5, aspect=5)
    ax3[0].set_zlim(*zlim)
    ax3[1].set_zlim(1, 1.6)

    ####################################################
    # Now read Nyeungu's dark adjusted data for comparison.
    ####################################################
    files = oci_ss.get_pace_oci_l1a_files(data_path, pattern_reference_h)
    with Dataset(files[0], 'r') as fid_cdf:
        science_group = fid_cdf.groups['science_data']  # Read the CCD raw data
        Sci_red_ds_h = science_group.variables['sci_red_ds'][:]
        Sci_red_ds_mean_h = Sci_red_ds_h.mean(axis=0)
        Sci_red_ds_std_h = Sci_red_ds_h.std(axis=0)

    # print("Sci_red_ds.mean(), Sci_red_ds.std() ", Sci_red_ds.mean(), Sci_red_ds.std())
    fig3, ax3 = plt.subplots(1, 2, subplot_kw={"projection": "3d"})
    fig3.suptitle("Mean Red CCD counts vs $\lambda$ and Angle \n{:}.".format(files[0].name))
    [a.set_xlabel('Wavelength (nm)') for a in ax3]
    [a.set_ylabel('Scan Angle (deg)') for a in ax3]
    ax3[0].set_title("Hyeungu's Mean CCD Counts over {:} scans.".format(Sci_red_ds_n.shape[0]))
    ax3[1].set_title("Hyeungu's STD CCD Counts over {:} scans.".format(Sci_red_ds_n.shape[0]))

    angle, lam = np.meshgrid(oci.pixel_angle, oci.red_wavelengths)
    if 1:
        # surf3 = ax3.plot_surface(lam, angle, Sci_red_ds_h[10, :, :], cmap='viridis',
        #                          edgecolor='none',
        #                          vmin=-1, vmax=1) #cmap=mpl.rcParams['image.cmap'])
    # else:
        surf3 = ax3[0].plot_surface(lam, angle, Sci_red_ds_mean_h, cmap='viridis',
                                    edgecolor='none',
                                    vmin=-0.1, vmax=0.1) #cmap=mpl.rcParams['image.cmap'])
        surf3b = ax3[1].plot_surface(lam, angle, Sci_red_ds_std_h, cmap='viridis',
                                    edgecolor='none',
                                    vmin=1.2, vmax=1.5) #cmap=mpl.rcParams['image.cmap'])
    fig3.colorbar(surf3, ax=ax3[0],
                  shrink=0.5, aspect=5)
    fig3.colorbar(surf3, ax=ax3[1],
                  shrink=0.5, aspect=5)
    ax3[0].set_zlim(*zlim)
    ax3[1].set_zlim(1, 1.6)

    # wavelengths = np.loadtxt(data_path/wavelength_file)
    # print("Wavelengths from csv file")
    # print(wavelengths)
    # print("Wavelengths from nc file")
    # print(oci.red_wavelengths)

    #####################################################
    # Now plot the difference between the three solutions
    #####################################################
    zlim = [-0.002, 0.006]
    # zlim = [-0.2, 0.06]
    # zlim = [-0.5, 0.5]

    fig4, ax4 = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
    fig4.suptitle("Difference between CTF and NRC's results vs $\lambda$ and Angle \n{:}.".format(files[0].name))
    ax4.set_xlabel('Wavelength (nm)')
    ax4.set_ylabel('Scan Angle (deg)')
    # ax4.set_zlabel("Nick's Background subtracted Counts")

    surf4 = ax4.plot_surface(lam, angle,  Sci_red_data_mean - Sci_red_ds_mean_n, cmap='viridis',
                             edgecolor='none',
                             vmin=-0.002, vmax=0.003) #cmap=mpl.rcParams['image.cmap'])
    fig4.colorbar(surf4, ax=ax4,
                  shrink=0.5, aspect=5)
    ax4.set_zlim(*zlim)

    fig4, ax4 = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
    fig4.suptitle("Difference between CTF and HChoi's results vs $\lambda$ and Angle \n{:}.".format(files[0].name))
    ax4.set_xlabel('Wavelength (nm)')
    ax4.set_ylabel('Scan Angle (deg)')
    # ax4.set_zlabel("Nick's Background subtracted Counts")

    surf4 = ax4.plot_surface(lam, angle,  Sci_red_data_mean - Sci_red_ds_mean_h, cmap='viridis',
                             edgecolor='none',
                             vmin=-0.002, vmax=0.003) #cmap=mpl.rcParams['image.cmap'])
    fig4.colorbar(surf4, ax=ax4,
                  shrink=0.5, aspect=5)
    ax4.set_zlim(*zlim)

    fig4, ax4 = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
    fig4.suptitle("Difference between NRC and HChoi's results vs $\lambda$ and Angle \n{:}.".format(files[0].name))
    ax4.set_xlabel('Wavelength (nm)')
    ax4.set_ylabel('Scan Angle (deg)')
    # ax4.set_zlabel("Nick's Background subtracted Counts")

    surf4 = ax4.plot_surface(lam, angle,  Sci_red_ds_mean_n - Sci_red_ds_mean_h, cmap='viridis',
                             edgecolor='none',
                             vmin=-0.002, vmax=0.003) #cmap=mpl.rcParams['image.cmap'])
    fig4.colorbar(surf4, ax=ax4,
                  shrink=0.5, aspect=5)
    ax4.set_zlim(*zlim)


    # path = "/Users/cfield/Documents/PACE/python_code/oci_l1a/test/dark_view_correction_examples/ncollins_cfield_PACE_OCI_DIAG_2020055140.20200224T184932.L1A.nc"
    # path = "/Users/cfield/Documents/PACE/python_code/oci_l1a/test/dark_view_correction_examples/ncollins_cfield_PACE_OCI_2020055138.20200224T183300.L1A.nc"
    # with Dataset(path, 'a') as fid_cdf:
    #     outVar = fid_cdf.createVariable('/science_data/sci_red_ds_ctf',
    #                                     fid_cdf['/science_data/sci_red_ds'].datatype,
    #                                     fid_cdf['/science_data/sci_red_ds'].dimensions)
    #     outVar[:] = oci.Sci_red_data
    plt.show()

