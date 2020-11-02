#!/usr/bin/env python
"""This module provides a class to read PACE OCI L1a spatial/spectral data
from one or more netCDF4 files. It also provides a method to return the data
spectral bands.

For a description of the L1A file format see:
     **OCI L1A  data format** TDMS Document number **OCI-SCI-SPEC-0079**.

A more descriptive document on the L1A format is:
    **PACE OCI Level-1A Data Product User’s Guide** by Fred Patt

For information on the thermal variables in the L1A files, see:
    **OCI-THRM-SPEC-0108**    and    **OCI-SYS-TN-0270**

"""

# To read/write Excel files see http://www.python-excel.org
# To create Power Point see
# https://python-pptx.readthedocs.io/en/latest/
# https://yoroshikune.com/create-powerpoint-presentation-python/
# https://pbpython.com/creating-powerpoint.html

import argparse
from copy import deepcopy
from pathlib import Path
from pprint import pprint
import sys

from matplotlib import pyplot as plt
import matplotlib as mpl            # Used to get color map.
# from mpl_toolkits.mplot3d import axes3d
from netCDF4 import Dataset         # https://unidata.github.io/netcdf4-python/netCDF4/index.html
import numpy as np
import numpy.ma as ma

import etu_thermal as therm

CCD_DIM = 512
"""Number of pixels in each direction of the CCD. Constant."""
NUM_TAPS = 16
# """Number of CCD readout taps. Constant."""
RED_S_ORIGIN = 597.5
# """Lowest wavelength for the red CCD (nm). Constant."""
BLUE_S_ORIGIN = 320
# """Lowest wavelength for the blue CCD (nm). Constant."""
DELTA_LAMBDA = 0.625
# """Spectral bandwidth of each CCD pixel (nm). Constant."""

def get_pace_oci_l1a_files(path, pattern=""):
    """Return from path a sorted list of PACE OCI L1A files with that match the pattern.

    :param path: Path to a folder with OCI L1A files.
    :param pattern: Center of pattern to match
    :return: Sorted list of files that match pattern.

    If pattern == "" (default), the file match is made against "PACE_OCI_*.L1A.nc".
    If pattern starts with "/", the file match is made against pattern[1:].
    Otherwise, the file match is made against "PACE_OCI_*" + pattern + "*.L1A.nc"
    """
    if not isinstance(path, Path):
        path = Path(path)
    if not path.exists():
        raise ValueError("Target path doesn't exit. {:}".format(path))
    if not path.is_dir():
        raise ValueError("Target path isn't a directory.. {:}".format(path))
    if pattern == "":
        pattern = "PACE_OCI_*.L1A.nc"
    elif pattern[0] == "/":
        pattern = pattern[1:]
    else:
        pattern = "PACE_OCI_*" + pattern + "*.L1A.nc"
    files = list(path.expanduser().glob(pattern))
    files.sort()
    return files


def getCCDbands3(spectral_mode, band):
    """compute PACE OCI CCD wavebands for given input aggregation spectral_mode.

    :param spectral_mode: vector of spectral aggregation factors, one for each tap
    :param band: "blue" or "red"
    :return: output vector of aggregated pixel center wavelengths, max(spectral_mode)

    :raises: ValueError **band** is not one of the valid strings.
    """
    # Inspired by IDL code WRITTEN BY Sam Kitchen-McKinley, SSAI

    CCDperTAP = CCD_DIM / NUM_TAPS
    if band == "red":
        origin = RED_S_ORIGIN
    elif band == "blue":
        origin = BLUE_S_ORIGIN
    else:
        raise ValueError("{:} band value invalid. Must be 'blue' or 'red'.".format(band))

    tap_starts = origin + CCDperTAP * DELTA_LAMBDA * np.arange(NUM_TAPS)
    wavelengths = [s + 0.5 * DELTA_LAMBDA * mode + DELTA_LAMBDA * np.arange(0, CCDperTAP, mode)
                   for s, mode in zip(tap_starts, spectral_mode) if mode != 0]
    wavelengths = np.concatenate(wavelengths, axis=0)
    if not np.all(spectral_mode == spectral_mode[0]):
        print("Warning. spectral_mode not all the same.",
              file=sys.stderr)
    return wavelengths, spectral_mode.max()


def get_CCD_pixel_type(spatial_zone_lines,
                       spatial_aggregation,
                       spatial_zone_data_type,
                       num_spatial_pixels):
    """Use spatial zone line, aggregation, and type to generate scanner angle for each pixel.

    :param spatial_zone_lines: Values from file "spatial_spectral_modes"
    :param spatial_aggregation: Values from file "spatial_spectral_modes"
    :param spatial_zone_data_type: Values from file "spatial_spectral_modes"
    :param num_spatial_pixels: Number of spatial pixels in the scan (CCD pixels).
    :return: 1darray pixel_angle, 1darray pixel_data_type, spatial aggregation
    """
    zone_span = spatial_zone_lines * 360. / spatial_zone_lines.sum()  # Angular span of each zone
    zone_start = np.zeros(len(zone_span) + 1, dtype=zone_span.dtype)
    zone_start[1:] = zone_span.cumsum()                         # Start angle of each zone
    # Compute pixels in each zone. This is a masked array with zones not represented in the data masked.
    pixels_zone = spatial_zone_lines // spatial_aggregation
    pixel_end = pixels_zone.cumsum()  # Ending pixel of each recorded zone.

    # From the line start index, assign a scan angle to each aggregated spatial pixel in the CCD image.
    # Identify each pixel's zone data type.
    pixel_angle = np.zeros((num_spatial_pixels,), dtype=np.float32)
    pixel_data_type = np.zeros((1, 1, num_spatial_pixels,), dtype=np.int16)
    for num, i_end, angle_range, angle_start, zone_type in zip(pixels_zone,
                                                               pixel_end,
                                                               zone_span,
                                                               zone_start,
                                                               spatial_zone_data_type):
        if num is ma.masked:                                    # If # pixels is masked, then no pixels
            continue                                            # Go to the next zone.
        # pixel_angle[i_end - num:i_end] = angle_start + angle_range / num * range(num)
        # pixel_data_type[:, :, i_end - num:i_end] = zone_type
        try:
            pixel_angle[i_end-num:i_end] = angle_start + angle_range/num * range(num)
            pixel_data_type[:, :, i_end-num:i_end] = zone_type
        except ValueError:
            if zone_type == 2:
                print("Warning. File appears to be old format without dark pixels in CCD image array.",
                      file=sys.stderr)
                continue
            raise
    # Get the spatial aggregation mode. Test that all the same
    a = spatial_aggregation[spatial_zone_data_type > 0]
    if not np.all(a == a[0]):
        print("Warning. Spatial aggregation not all the same in captured pixels.",
              file=sys.stderr)
    return pixel_angle, pixel_data_type, a[0]


class READ_OCI_SPA_SPE_L1A(object):
    """Class to read the OCI l1a spatial/spectral data and assign spectral and pixel labels.

        :param path: File path, either a Path object or string. Default (None)
        :param bool mask_dark: Set True to mask dark data in returned CCD image. Default is False.

        If **path** is None, then an empty object is created to which files can be appended
        by using the **append** method.

        The goal of class is to create data attributes of the OCI RED CCD output and BLUE CCD output.
        The CCD outputs are **N_s** x **N_w** x **N_p** arrays, where **N_s** is the number of
        scans, **N_w** is the number of spectral bands (**N_wr** for the red CCD and **N_wb**
        for the blue CCD), and **N_p** is the number of spatial pixels.

        Data attributes created by this class fall into two groups (explained below) and are:

        Group 1:
           #. **aux_param_table**: The aux_param_table read from the file.
           #. **blue_spectral_mode**: The blue_spectral_mode values read from the file.
           #. **red_spectral_mode**: The red_spectral_mode values read from the file.
           #. **spatial_aggregation**: The spatial_aggregation values read from the file.
           #. **spatial_zone_data_type**: The spatial_zone_data_type values read from the file.
           #. **spatial_zone_lines**: The spatial_zone_lines values read from the file.

           #. **red_wavelengths**: The center wavelength of each red CCD spectral band (nm) (Length N_wr).
           #. **blue_wavelengths**: The center wavelength of each blue CCD spectral band (nm) (Length N_wb).
           #. **pixel_angle**: The scan angle from the start marker of each spatial pixel (deg) (Length N_p).
           #. **pixel_data_type**: The spatial zone type of each spatial pixel (unitless) (Length N_p).

        Group 2:
           #. **pixel_data_type_red**: Same as pixel_data_type broadcast to shape
              **N_s** x **N_wr** x **N_p**. Note that due to the magic of numpy indexing this array
              uses only the memory pixel_data_type does.
           #. **pixel_data_type_blue**: Same as pixel_data_type broadcast to shape
              **N_s** x **N_wr** x **N_p**. Note that due to the magic of numpy indexing this array
              uses only the memory as pixel_data_type does.
           #. **DC_swir_data**: The dark count data read from DC_swir (**N_s** x **N_wr** x **N_d**)?,
              where **N_d** is the number of dark pixels recorded.
           #. **DC_red_data**: The dark count data read from DC_red (**N_s** x **N_wr** x **N_d**),
              where **N_d** is the number of dark pixels recorded. Values have been divided by 16.
           #. **DC_blue_data**: The dark count data read from DC_blue (**N_s** x **N_wb** x **N_d**).
              Values have been divided by 16.
           #. **Sci_swir_data**: The raw counts read from sci_swir (**N_s** x **N_wr** x **N_d**)?.
              If **mask_dark** was True, then any pixels not in zone 9 are masked. (NOT YET.)
           #. **Sci_red_data**: The raw counts read from sci_red (**N_s** x **N_wr** x **N_d**).
              If **mask_dark** was True, then any pixels not in zone 9 are masked.
           #. **Sci_blue_data**: The raw counts read from sci_blue (**N_s** x **N_wb** x **N_d**).
              If **mask_dark** was True, then any pixels not in zone 9 are masked.
           #. **scan_start_time**: The values read from scan_start_time (**N_s**).
           #. **scan_start_CCSDS_sec**: The values read from scan_start_CCSDS_sec (**N_s**).
           #. **scan_start_CCSDS_usec**: The values read from scan_start_CCSDS_usec (**N_s**).
           #. **scan_start_TIA**: Scan line start time in TAI (int64 in us) (**N_s**).

        In addition, **path** is defined, the path to a single input file or a list of paths to
        each input file when **append()** is used.

        When data files are accumulated by the append method, the first 6 data attributes in
        group 1 are tested for a match between the previously read file and the current one.
        Any mismatch will raise a ValueError indicating which was the first one to not match.
        Data attributes in Group 1 are all taken from the current file and are always valid.
        However, when appending, data attributes in Group 2 are accumulated in lists
        and don't take the form of numpy masked arrays until after the **append_conclude()**
        method has been called.

        .. note::
            Some methods permanently modify the data attributes in place. To preserve a version of
            the original use **copy.deepcopy()**. For example::

                from copy import deepcopy
                oci.append(read_oci_1a_class(file))
                oci2 = deepcopy(oci)
                oci2.subtract_dark()
                oci2 = deepcopy(oci)    # Reset oci2 to the original values.
    """
    def __init__(self, path=None, mask_dark=False):
        """Constructor for class READ_OCI_SPA_SPE_L1A"""

        self.path = path
        if path is None:
            return
        with Dataset(path,'r') as fid_cdf:
            dark_group = fid_cdf.groups['onboard_calibration_data']             # Read dark data
            self.DC_swir_data = dark_group.variables['DC_SWIR'][:]
            self.DC_red_data = dark_group.variables['DC_red'][:]
            self.DC_blue_data = dark_group.variables['DC_blue'][:]

            science_group = fid_cdf.groups['science_data']                          # Read the CCD raw data
            self.Sci_swir_data =  science_group.variables['sci_SWIR'][()]
            self.Sci_red_data =  science_group.variables['sci_red'][()]
            self.Sci_blue_data =  science_group.variables['sci_blue'][()]

            time_group = fid_cdf.groups['scan_line_attributes']                     # Get scan line attributes
            self.scan_start_time = time_group.variables['scan_start_time'][()]
            self.scan_start_CCSDS_sec = time_group.variables['scan_start_CCSDS_sec'][()]
            self.scan_start_CCSDS_usec = time_group.variables['scan_start_CCSDS_usec'][()]

            # Need the spatial/spectral aggregation first so can calculate dark shift.
            spatial_spectral_group = fid_cdf.groups['spatial_spectral_modes']  # Get aggregation attributes
            self.aux_param_table = spatial_spectral_group.variables['aux_param_table'][()]
            self.blue_spectral_mode = spatial_spectral_group.variables['blue_spectral_mode'][()]
            self.red_spectral_mode = spatial_spectral_group.variables['red_spectral_mode'][()]
            self.spatial_aggregation = spatial_spectral_group.variables['spatial_aggregation'][()]
            self.spatial_zone_data_type = spatial_spectral_group.variables['spatial_zone_data_type'][()]
            self.spatial_zone_lines = spatial_spectral_group.variables['spatial_zone_lines'][()]

        # Compute the wavelengths from the spectral modes.
        self.red_wavelengths, red_agg = getCCDbands3(self.red_spectral_mode, 'red')
        self.blue_wavelengths, blue_agg = getCCDbands3(self.blue_spectral_mode, 'blue')
        # Map the Sci_?_data spatial (third index) to scan angle
        self.pixel_angle, self.pixel_data_type, spat_agg = get_CCD_pixel_type(self.spatial_zone_lines,
                                                                              self.spatial_aggregation,
                                                                              self.spatial_zone_data_type,
                                                                              self.Sci_red_data.shape[2])
        self.pixel_data_type_red = np.broadcast_to(self.pixel_data_type, self.Sci_red_data.shape)
        self.pixel_data_type_blue = np.broadcast_to(self.pixel_data_type, self.Sci_blue_data.shape)
        if mask_dark:
            self.Sci_red_data = ma.masked_where(self.pixel_data_type_red != 9, self.Sci_red_data)
            self.Sci_blue_data = ma.masked_where(self.pixel_data_type_blue != 9, self.Sci_blue_data)
        # Now shift the dark counts by the required amount.
        # The number of bits to shift based on the pixels aggregated.
        dark_shift = {1: 0, 2: 0, 4: 0, 8: 1, 16: 2, 32: 3, 64: 4}
        self.DC_red_data = np.right_shift(self.DC_red_data, dark_shift[red_agg * spat_agg])
        self.DC_blue_data = np.right_shift(self.DC_blue_data, dark_shift[blue_agg * spat_agg])

        # UNIX - double, seconds     , epoch 0h Jan 1 1970, adjusted for leapseconds
        # TAI  - double, seconds     , epoch 0h Jan 1 1958, ignores leapseconds
        # Compute the scan start time in TAI with microsecond resolution
        self.scan_start_CCSDS_sec = self.scan_start_CCSDS_sec.astype(np.int64)
        self.scan_start_TIA = self.scan_start_CCSDS_sec * 1000000 + self.scan_start_CCSDS_usec  # usec since TAI


    def append(self, other):
        """Append the CCD data from "other" to this object.

        :param other: A READ_OCI_SPA_SPE_L1A object
        :return: None

        :raises: ValueError if any of the spatial/spectral aggregation are different.

        Typical use::

            oci_ccd = READ_OCI_SPA_SPE_L1A()
            for file in files:
                oci_ccd.append(READ_OCI_SPA_SPE_L1A(file))
            oci_ccd.append_conclude()
        """
        if self.path is None:
            self.path = []                              # Create all of the lists.
            self.DC_swirdata_lst, self.Sci_swir_data_lst = [], []
            self.DC_red_data_lst, self.Sci_red_data_lst = [], []
            self.DC_blue_data_lst, self.Sci_blue_data_lst = [], []
            self.scan_start_time_lst, self.scan_start_CCSDS_sec_lst, self.scan_start_CCSDS_usec_lst = [], [], []
            self.aux_param_table, self.blue_spectral_mode, self.red_spectral_mode = [], [], []
            self.spatial_aggregation, self.spatial_zone_data_type, self.spatial_zone_lines = [], [], []
            self.scan_start_CCSDS_sec_lst, self.scan_start_TIA_lst = [], []
        else:                                           # Test match between old and new parameters.
            if not np.array_equal(self.aux_param_table, other.aux_param_table): #, equal_nan=True):
                raise ValueError("aux_param_table values are different.")
            if not np.array_equal(self.blue_spectral_mode, other.blue_spectral_mode): #, equal_nan=True):
                raise ValueError("blue_spectral_mode values are different.")
            if not np.array_equal(self.red_spectral_mode, other.red_spectral_mode): #, equal_nan=True):
                raise ValueError("red_spectral_mode values are different.")
            if not np.array_equal(self.spatial_aggregation, other.spatial_aggregation): #, equal_nan=True):
                raise ValueError("spatial_aggregation values are different.")
            if not np.array_equal(self.spatial_zone_data_type, other.spatial_zone_data_type): #, equal_nan=True):
                raise ValueError("spatial_zone_data_type values are different.")
            if not np.array_equal(self.spatial_zone_lines, other.spatial_zone_lines): #, equal_nan=True):
                raise ValueError("spatial_zone_lines values are different.")

        # Assign those values of "other" that (should) be the same as "self" to self.
        self.path.append(other.path)
        self.aux_param_table = other.aux_param_table

        self.red_spectral_mode =  other.red_spectral_mode
        self.blue_spectral_mode = other.blue_spectral_mode
        self.red_wavelengths = other.red_wavelengths
        self.blue_wavelengths = other.blue_wavelengths

        self.spatial_aggregation = other.spatial_aggregation
        self.spatial_zone_data_type = other.spatial_zone_data_type
        self.spatial_zone_lines = other.spatial_zone_lines
        self.pixel_angle = other.pixel_angle
        self.pixel_data_type = other.pixel_data_type

        # Append each data element to the corresponding lists.
        self.DC_swirdata_lstdata_lst.append(other.DC_swir_data_data)
        self.DC_red_data_lst.append(other.DC_red_data)
        self.DC_blue_data_lst.append(other.DC_blue_data)
        self.Sci_swir_data_lst.append(other.Sci_swir_data)
        self.Sci_red_data_lst.append(other.Sci_red_data)
        self.Sci_blue_data_lst.append(other.Sci_blue_data)
        self.scan_start_time_lst.append(other.scan_start_time)
        self.scan_start_CCSDS_sec_lst.append(other.scan_start_CCSDS_sec)
        self.scan_start_CCSDS_usec_lst.append(other.scan_start_CCSDS_usec)
        self.scan_start_TIA_lst.append(other.scan_start_TIA)


    def append_conclude(self):
        """Convert lists accumulated by **append** into masked arrays in group 2.

        :return: None

        Until this method is called, the data attributes in Group 2 are not valid.
        """
        # Concatenate the lists into arrays, map pixel_data_type to the CCD arrays,
        # and delete the lists that were used to accumulate to de-clutter the namespace.
        self.DC_swir_data = ma.concatenate(self.DC_swirdata_lst, axis=0)
        self.DC_red_data = ma.concatenate(self.DC_red_data_lst, axis=0)
        self.DC_blue_data = ma.concatenate(self.DC_blue_data_lst, axis=0)
        self.Sci_swirdata = ma.concatenate(self.Sci_swir_data_lst, axis=0)
        self.Sci_red_data = ma.concatenate(self.Sci_red_data_lst, axis=0)
        self.Sci_blue_data = ma.concatenate(self.Sci_blue_data_lst, axis=0)
        self.scan_start_time = ma.concatenate(self.scan_start_time_lst, axis=0)
        self.scan_start_CCSDS_sec = ma.concatenate(self.scan_start_CCSDS_sec_lst, axis=0)
        self.scan_start_CCSDS_usec = ma.concatenate(self.scan_start_CCSDS_usec_lst, axis=0)
        self.scan_start_TIA = ma.concatenate(self.scan_start_TIA_lst, axis=0)

        self.pixel_data_type_red = np.broadcast_to(self.pixel_data_type,
                                                   self.Sci_red_data.shape)
        self.pixel_data_type_blue = np.broadcast_to(self.pixel_data_type,
                                                    self.Sci_blue_data.shape)
        del self.DC_swir_data_lst,  self.DC_red_data_lst,  self.DC_blue_data_lst
        del self.Sci_swir_data_lst,  self.Sci_red_data_lst,  self.Sci_blue_data_lst
        del self.scan_start_time_lst,  self.scan_start_TIA_lst
        del self.scan_start_CCSDS_sec_lst, self.scan_start_CCSDS_usec_lst


    def subtract_dark(self, start=40, dev=3):
        """Average the dark counts for each scan and subtract from each scan

        :param start: First index to use when computing the dark average (Default=40)
        :param dev: Multiple of sigma to mask when computing dark level (Default=3)
        :return:

        .. note::
            This method permanently modifies **Sci_red_data** and **Sci_blue_data** in memory.
        """
        DC_swir_data = self.DC_swir_data[:, :, start:].astype(np.float64)
        DC_swir_mean = DC_swir_data.mean(axis=2, keepdims=True)
        DC_swir_std_ = DC_swir_data.std(axis=2, keepdims=True)
        DC_swir_mean = ma.masked_where((DC_swir_data < DC_swir_mean - dev * DC_swir_std_)
                                    | (DC_swir_data > DC_swir_mean + dev * DC_swir_std_),
                                      DC_swir_data).mean(axis=2, keepdims=True)
        self.Sci_swir_data = self.Sci_swir_data.astype(np.float64) - DC_swir_mean

        DC_red_data = self.DC_red_data[:, :, start:].astype(np.float64)
        DC_red_mean = DC_red_data.mean(axis=2, keepdims=True)
        DC_red_std_ = DC_red_data.std(axis=2, keepdims=True)
        DC_red_mean = ma.masked_where((DC_red_data < DC_red_mean - dev * DC_red_std_)
                                    | (DC_red_data > DC_red_mean + dev * DC_red_std_),
                                      DC_red_data).mean(axis=2, keepdims=True)
        self.Sci_red_data = self.Sci_red_data.astype(np.float64) - DC_red_mean

        DC_blue_data = self.DC_blue_data[:, :, start:].astype(np.float64)
        DC_blue_mean = DC_blue_data.mean(axis=2, keepdims=True)
        DC_blue_std_ = DC_blue_data.std(axis=2, keepdims=True)
        DC_blue_mean = ma.masked_where((DC_blue_data < DC_blue_mean - dev * DC_blue_std_)
                                    | (DC_blue_data > DC_blue_mean + dev * DC_blue_std_),
                                      DC_blue_data).mean(axis=2, keepdims=True)
        self.Sci_blue_data = self.Sci_blue_data.astype(np.float64) - DC_red_mean


    def select__data_type(self, data_type):
        """Trim data attributes to those spatial pixels with the type **data_type**.

        :param data_type: Scalar data type to keep in the trimmed data.
        :return: None

        .. note::
            This method permanently modifies the data attributes in place.
        """
        index = self.pixel_data_type[0, 0, :] == data_type
        self.Sci_red_data = self.Sci_red_data[:, :, index]
        self.Sci_blue_data = self.Sci_blue_data[:, :, index]
        self.pixel_data_type = self.pixel_data_type[:, :, index]    # shape = (1, 1, n)
        self.pixel_angle = self.pixel_angle[index]
        self.pixel_data_type_red = np.broadcast_to(self.pixel_data_type,
                                                   self.Sci_red_data.shape)
        self.pixel_data_type_blue = np.broadcast_to(self.pixel_data_type,
                                                    self.Sci_blue_data.shape)


class READ_OCI_THERMAL_L1A(object):
    """Class to read the OCI l1a thermal data.

        :param path: File path, either a Path object or string. Default (None)
        :param bool mask_dark: Set True to mask dark data in returned CCD image. Default is False.

        If **path** is None, then an empty object is created to which files can be appended
        by using the **append** method.

        Data attributes created by this class are:

           #. **aob_temperatures**
           #. **blue_fpa_temperatures**
           #. **dau_temperatures**
           #. **dau_tlm_time**
           #. **icdu_temperatures**.  Not implemented.
           #. **red_fpa_temperatures**
           #. **sds_det_temperatures**
           #. **tc_tlm_time**

        In addition, **path** is defined, the path to a single input file or a list of paths to
        each input file when **append()** is used.
    """


    def __init__(self, path=None, mask_dark=False):
        """Constructor for class READ_OCI_THERMAL_L1A"""

        self.path = path
        if path is None:
            return
        with Dataset(path,'r') as fid_cdf:
            # print("fid_cdf.groups")
            # pprint(fid_cdf.groups)
            # pprint(fid_cdf.dimensions)
            # pprint(fid_cdf.variables)
            #
            # print("fid_cdf.groups['engineering_data'].groups")
            # pprint(fid_cdf.groups['engineering_data'].groups)
            # pprint(fid_cdf.groups['engineering_data'].dimensions)
            # pprint(fid_cdf.groups['engineering_data'].variables)

            engineering_data = fid_cdf.groups['engineering_data']
            aob_cnts = engineering_data.variables['AOB_temperatures'][()]
            blue_fpa_cnts = engineering_data.variables['blue_FPA_temperatures'][()]
            dau_cnts = engineering_data.variables['DAU_temperatures'][()]
            self.dau_tlm_time = engineering_data.variables['DAU_tlm_time'][()]      # Time axis
            # icdu_cnts = engineering_data.variables['ICDU_thermisters'][()]
            red_fpa_cnts = engineering_data.variables['red_FPA_temperatures'][()]
            sds_det_cnts = engineering_data.variables['SDS_det_temperatures'][()]
            self.tc_tlm_time = engineering_data.variables['TC_tlm_time'][()]        # Time axis

            SCALE = 2.5 / 32768.                            # Convert all counts to volts
            aob_volts =      SCALE * aob_cnts
            blue_fpa_volts = SCALE * blue_fpa_cnts
            dau_volts =      SCALE * dau_cnts
            # icdu_volts =     SCALE * icdu_cnts
            red_fpa_volts =  SCALE * red_fpa_cnts
            sds_det_volts =  SCALE * sds_det_cnts

            # Use Herner's method to evaluate the polynomials to convert to deg C
            self.aob_temperatures = np.zeros_like(aob_volts, dtype=np.float64)
            coeffs = np.asarray(therm.DAU_AOB_temps["convert"], dtype=np.float64).T
            self.aob_temperatures[:] = coeffs[0, :].copy()
            for c in coeffs[1:, :]:
                self.aob_temperatures *= aob_volts
                self.aob_temperatures += c
            # print("Converting aob_volts")
            # print(self.aob_temperatures)

            self.blue_fpa_temperatures = np.zeros_like(blue_fpa_cnts, dtype=np.float64)
            self.red_fpa_temperatures = np.zeros_like(red_fpa_volts, dtype=np.float64)
            coeffs = np.asarray(therm.DAU_FPA_temps["convert"], dtype=np.float64).T
            self.blue_fpa_temperatures[:] = coeffs[0, :].copy()
            self.red_fpa_temperatures[:] = coeffs[0, :].copy()
            for c in coeffs[1:, :]:
                self.blue_fpa_temperatures *= blue_fpa_volts
                self.blue_fpa_temperatures += c
                self.red_fpa_temperatures *= red_fpa_volts
                self.red_fpa_temperatures += c
            # print("Converting blue_fpa_cnts and red_fpa_cnts")
            # print(self.red_fpa_temperatures)

            self.dau_temperatures = np.zeros_like(dau_volts, dtype=np.float64)
            coeffs = np.asarray(therm.DAU_temps["convert"], dtype=np.float64).T
            self.dau_temperatures[:] = coeffs[0, :].copy()
            for c in coeffs[1:, :]:
                self.dau_temperatures *= dau_volts
                self.dau_temperatures += c
            # print("Converting dau_volts")
            # print(self.dau_temperatures)

            # self.icdu_temperatures = np.zeros_like(icdu_volts, dtype=np.float64)
            # coeffs = np.asarray(therm.DAU_temps["convert"], dtype=np.float64).T
            # self.icdu_temperatures[:] = coeffs[0, :].copy()
            # for c in coeffs[1:, :]:
            #     self.icdu_temperatures *= icdu_volts
            #     self.icdu_temperatures += c
            # print("Converting icdu_volts")
            # print(self.icdu_temperatures)

            self.sds_det_temperatures = np.zeros_like(sds_det_volts, dtype=np.float64)
            coeffs = np.asarray(therm.DAU_SDS_det_temps["convert"], dtype=np.float64).T
            self.sds_det_temperatures[:] = coeffs[0, :].copy()
            for c in coeffs[1:, :]:
                self.sds_det_temperatures *= sds_det_volts
                self.sds_det_temperatures += c
            # print("Converting sds_det_volts")
            # print(self.sds_det_temperatures)


    def append(self, other):
        """Append the CCD data from "other" to this object.

        :param other: A read_oci_1a_class object
        :return: None

        :raises: ValueError if any of the spatial/spectral aggregation are different.

        Typical use::

            oci_therm = READ_OCI_THERMAL_L1A()
            for file in files:
                oci_therm.append(READ_OCI_THERMAL_L1A(file))
            oci_therm.append_conclude()
        """

        if self.path is None:
            self.path = []
            self.aob_temperatures_lst, self.blue_fpa_temperatures_lst = [], []
            self.dau_temperatures_lst, self.dau_tlm_time_lst = [], []
            self.icdu_temperatures_lst, self.red_fpa_temperatures_lst = [], []
            self.sds_det_temperatures_lst, self.tc_tlm_time_lst = [], []

        self.path.append(other.path)
        self.aob_temperatures_lst.append(other.aob_temperatures)
        self.blue_fpa_temperatures_lst.append(other.blue_fpa_temperatures)
        self.dau_temperatures_lst.append(other.dau_temperatures)
        self.dau_tlm_time_lst.append(other.dau_tlm_time)
        # self.icdu_temperatures_lst.append(other.icdu_temperatures)
        self.red_fpa_temperatures_lst.append(other.red_fpa_temperatures)
        self.sds_det_temperatures_lst.append(other.sds_det_temperatures)
        self.tc_tlm_time_lst.append(other.tc_tlm_time)


    def append_conclude(self):
        """Convert lists accumulated by **append** into masked arrays in group 2.

        :return: None

        Until this method is called, the thermal numpy arrays are not defined.
        """
        # Concatenate the lists into arrays.
        self.aob_temperatures = ma.concatenate(self.aob_temperatures_lst, axis=0)
        self.blue_fpa_temperatures = ma.concatenate(self.blue_fpa_temperatures_lst, axis=0)
        self.dau_temperatures = ma.concatenate(self.dau_temperatures_lst, axis=0)
        self.dau_tlm_time = ma.concatenate(self.dau_tlm_time_lst, axis=0)
        # self.icdu_temperatures = ma.concatenate(self.icdu_temperatures_lst, axis=0)
        self.red_fpa_temperatures = ma.concatenate(self.red_fpa_temperatures_lst, axis=0)
        self.sds_det_temperatures = ma.concatenate(self.sds_det_temperatures_lst, axis=0)
        self.tc_tlm_time = ma.concatenate(self.tc_tlm_time_lst, axis=0)

        del self.aob_temperatures_lst,  self.blue_fpa_temperatures_lst
        del self.dau_temperatures_lst
        del self.dau_tlm_time_lst,  self.red_fpa_temperatures_lst
        del self.sds_det_temperatures_lst, self.tc_tlm_time_lst


    def plot(self, variable, ax=None, plt_type='3d'):
        """Method doesn't work. Just return"""
        # return

        """Create a 3d plot of the specified thermal data.

        :param variable: string name of the variable to plot
        :param ax: If not None, 3d axis on which to place the plot. If None (default) make a new figure.
        :param plt_type: Either '3d' for 'imag' to specify the type of display. (Default '3d') 
        :return:

        Valid string values for **variable** are:
           - "aob_temperatures"
           - "blue_fpa_temperatures"
           - "dau_temperatures"
           - "red_fpa_temperatures"
           - "sds_det_temperatures"

        If an axis is provided, it must have been created with projection='3d'
        for a 3d plot and must be created without a projection for an image.
        """
        if plt_type not in ['3d', 'imag']:
            raise ValueError("Plot type must be either '3d' or 'imag'. {:} is neither.".format(plt_type))

        if variable == "aob_temperatures":
            title = "AOB_temperatures"
            time = self.dau_tlm_time
            labels = therm.DAU_AOB_temps["name"]
            temps = self.aob_temperatures.filled(np.nan)
        elif variable == "blue_fpa_temperatures":
            title = "blue_fpa_temperatures"
            time = self.dau_tlm_time
            labels = therm.DAU_FPA_temps["name"]
            temps = self.blue_fpa_temperatures.filled(np.nan)
        elif variable == "dau_temperatures":
            title = "dau_temperatures"
            time = self.dau_tlm_time
            labels = therm.DAU_temps["name"]
            temps = self.dau_temperatures.filled(np.nan)
        elif variable == "red_fpa_temperatures":
            title = "red_fpa_temperatures"
            time = self.dau_tlm_time
            labels = therm.DAU_FPA_temps["name"]
            temps = self.red_fpa_temperatures.filled(np.nan)
        elif variable == "sds_det_temperatures":
            title = "sds_det_temperatures"
            time = self.dau_tlm_time
            labels = therm.DAU_SDS_det_temps["name"]
            temps = self.sds_det_temperatures.filled(np.nan)
        else:
            raise ValueError("Variable given {:} not recognized.".format(variable))

        t_min, t_max = np.nanmin(temps), np.nanmax(temps)
        # If needed, create the figure and axis in the proper form.
        subplot_kw = {"projection": '3d'} if plt_type == '3d' else {}
        if ax is None:
            fig, ax = plt.subplots(1, 1 , subplot_kw=subplot_kw)
        ax.set_title(title)
        x_ticks = range(len(labels))        # Tick and label each channel
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(labels)
        fig.autofmt_xdate(rotation=90)
        if plt_type == 'imag':              # Make the plot
            surf = ax.imshow(temps)
            ax.set_aspect("auto")
        else:

            l, t = np.meshgrid(range(len(time)), x_ticks)
            surf = ax.plot_surface(t, l, temps.T,
                                   vmin=t_min, vmax=t_max,
                                   cmap=mpl.rcParams['image.cmap'])
            ax.set_aspect("auto")
        fig.colorbar(surf, ax=ax)

def main(args):
    """Read the PACE OCI L1A netCDF files in args.paths and plot termperatures and spectrum.

    :param args: Name space from argparse.
    :return: None
    """

    # Read the temperature data from the input files. For each
    # temperature variable, generate 2 plots. A standard image plot (in the
    # code below) and one using the thermal.plot method. They can be compared.

    # A single file can be read using this next line or a single file/multiple
    # files can be loaded as done in this code.
    # thermal = READ_OCI_THERMAL_L1A(paths[0], mask_dark=True)

    thermal = READ_OCI_THERMAL_L1A()  # Create empty instance for appending.
    for file in paths:
        thermal.append(READ_OCI_THERMAL_L1A(file, mask_dark=True))
    thermal.append_conclude()  # Convert the internal lists into arrays.

    plt_type = args.plt_type

    fig, ax = plt.subplots()            # Image generated here.
    fig.suptitle("thermal.aob_temperatures")
    ax.imshow(thermal.aob_temperatures)
    thermal.plot("aob_temperatures", plt_type=plt_type)  # Use built in method.

    fig, ax = plt.subplots()            # Image generated here.
    fig.suptitle("thermal.blue_fpa_temperatures")
    ax.imshow(thermal.blue_fpa_temperatures)
    thermal.plot("blue_fpa_temperatures", plt_type=plt_type)  # Use built in method.

    fig, ax = plt.subplots()            # Image generated here.
    fig.suptitle("thermal.dau_temperatures")
    ax.imshow(thermal.dau_temperatures)
    thermal.plot("dau_temperatures", plt_type=plt_type)  # Use built in method.

    # fig, ax = plt.subplots()            # Image generated here.
    # fig.suptitle("thermal.icdu_temperatures")
    # ax.imshow(thermal.icdu_temperatures)    # Use built in method.

    fig, ax = plt.subplots()            # Image generated here.
    fig.suptitle("thermal.red_fpa_temperatures")
    ax.imshow(thermal.red_fpa_temperatures)
    thermal.plot("red_fpa_temperatures", plt_type=plt_type)  # Use built in method.

    fig, ax = plt.subplots()            # Image generated here.
    fig.suptitle("thermal.sds_det_temperatures")
    ax.imshow(thermal.sds_det_temperatures)
    thermal.plot("sds_det_temperatures", plt_type=plt_type)  # Use built in method.

    # Read the spatial/spectral data from the input files.
    oci = READ_OCI_SPA_SPE_L1A(None, mask_dark=True)
    for file in paths:
        print(file)
        oci.append(READ_OCI_SPA_SPE_L1A(file, mask_dark=True))
    oci.append_conclude()

    oci2 = deepcopy(oci)                # Keep a copy of the original data
    oci2.select__data_type(9)           # Clip the spatial pixels to zone 9
    oci2.subtract_dark(start=50)
    # Mask scan start times less than equal to zero.
    scan_start_time = ma.masked_where(oci2.scan_start_time <= 0, oci2.scan_start_time)

    # Plot results from the original data
    fig, ax = plt.subplots(2, 2, sharex=True)
    fig.suptitle("Red CCD counts vs time from all files. Original data.")
    ax[0, 0].plot(scan_start_time, oci.Sci_red_data[:, :, 0], 'x')
    ax[0, 0].set_title("Sci_red[:, :, 0] vs time")
    ax[1, 0].plot(scan_start_time, oci.Sci_red_data[:, :, 320], 'x')
    ax[1, 0].set_title("Sci_red[:, :, 400] vs time")
    ax[0, 1].plot(scan_start_time, oci.Sci_red_data[:, 0, :], 'x')
    ax[0, 1].set_title("Sci_red[:, 0, :] vs time")
    ax[1, 1].plot(scan_start_time, oci.Sci_red_data[:, 320, :], 'x')
    ax[1, 1].set_title("Sci_red[:, 400, :] vs time")

    Sci_red_data_mean = oci.Sci_red_data.mean(axis=0)  # Compute mean across scan lines
    Sci_red_data_std = oci.Sci_red_data.std(axis=0)
    print("oci.Sci_red_data.shape, Sci_red_data_mean.shape")
    print(oci.Sci_red_data.shape, Sci_red_data_mean.shape)

    fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
    fig.suptitle("Red CCD counts vs time from all files. Original data.")
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Scan Angle')
    ax.set_zlabel('CCD Counts')
    oci.pixel_angle[oci.pixel_angle > 200] = np.nan  # Restrict angles to FOV.
    angle, lam = np.meshgrid(oci.pixel_angle, oci.red_wavelengths)
    ax.plot_surface(lam, angle, ma.filled(Sci_red_data_mean, np.nan))

    # Plot results from the background subtracted data
    Sci_red_data_mean = oci2.Sci_red_data.mean(axis=0)  # Compute mean across scan lines
    fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
    fig.suptitle("Red CCD counts vs time from all files. Background suppressed data.")
    ax.set_xlabel('Wavelength')
    ax.set_ylabel('Scan Angle')
    ax.set_zlabel('CCD Counts')
    # oci2.pixel_angle[oci2.pixel_angle>200] = np.nan   # Not needed because clipped to zone 9
    angle, lam = np.meshgrid(oci2.pixel_angle, oci2.red_wavelengths)
    ax.plot_surface(lam, angle, ma.filled(Sci_red_data_mean, np.nan))
    plt.show()


if __name__ == "__main__":

    paths = ['/Users/cfield/Downloads/PACE_OCI_SPEC_2020076013.20200316T132700.L1A.nc',
            '/Users/cfield/Downloads/PACE_OCI_SPEC_2020076013.20200316T132800.L1A.nc',
            '/Users/cfield/Downloads/PACE_OCI_SPEC_2020076013.20200316T133100.L1A.nc']

    # paths = get_pace_oci_l1a_files("/Users/cfield/Downloads")
    pprint(paths)

    parser = argparse.ArgumentParser(
        description='Test of class A class to read the spatial/spectral data from an OCI l1a file and return it with spectral and pixel lables.',
        epilog=" "
    )

    parser.add_argument('paths', nargs="+",
                        help='Path to files to read. Multiple files can be used '
                             'with command line wildcards.')
    parser.add_argument('--plt-type', '-p', choices=['3d', 'imag'], default='3d',
                        help='Type of temperature plots desired. Default is '
                             '3d. But can provide either.')
    args = parser.parse_args()
    print(args)
    main(args)
