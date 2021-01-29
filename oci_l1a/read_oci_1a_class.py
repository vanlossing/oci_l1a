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
import datetime
from pathlib import Path
from pprint import pprint
import sys

from matplotlib import pyplot as plt
import matplotlib as mpl            # Used to get color map and formatter
# from mpl_toolkits.mplot3d import axes3d
from netCDF4 import Dataset         # https://unidata.github.io/netcdf4-python/netCDF4/index.html
import numpy as np
import numpy.ma as ma

# pprint(sys.path)
# print(__file__)
try:                                # Address differences between PyCharm and command line
    from . import etu_thermal as therm
except ImportError:
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

PPR2NADIR = 110  # Angle (deg) between pulse per revolution index and nadir
# """Nominal degrees between the PPR mark and nadir"""

TIA2UTC = 37        # Subtract this from TIA to yield UTC times (sec).
# """Number of seconds between TIA and UTC since 31 December 2016"""


def get_pace_oci_l1a_files(path, pattern=""):
    """Return from path a sorted list of PACE OCI L1A files that match the pattern.

    :param path: Path to a folder with OCI L1A files.
    :param pattern: Center of pattern to match
    :return: Sorted list of full paths that match pattern.

    If pattern == "" (default), the file match is made against "PACE_OCI_*.L1A.nc".
    If pattern starts with "/", the file match is made against pattern[1:].
    Otherwise, the file match is made against "PACE_OCI_*" + pattern + "*.L1A.nc"

    For example,
    line 2 of the following will return a sorted list of the full
    path of all files that match "PACE_OCI_*.L1A.nc".
    Line 3 returns a sorted list of the full path of all files that match
    "PACE_OCI_20200102.nc" exactly; in this case, a single file.
    Line 4 will return a list of all files generated on 2020 January 05.

    .. code-block:: python
       :linenos:

        import oci_l1a.read_oci_1a_class as oci_ss
        path = "/path/to/my/data/"
        files = oci_ss.get_pace_oci_l1a_files(path)
        files = oci_ss.get_pace_oci_l1a_files(path, "PACE_OCI_20200102.nc")
        files = oci_ss.get_pace_oci_l1a_files(path, "20200105")

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
    elif band == "swir":
        return np.asarray([940., 1038., 1250., 1250.,  1378., 1615.,  1615.,  2130., 2260]), None
        __ =            """940,  1038,  1250SG,1250HG, 1378,  1614SG, 1615HG, 2130,  2260 """
    else:
        raise ValueError("{:} band value invalid. Must be 'blue' or 'red'.".format(band))

    tap_starts = origin + CCDperTAP * DELTA_LAMBDA * np.arange(NUM_TAPS)
    wavelengths = [s + 0.5 * DELTA_LAMBDA * mode + DELTA_LAMBDA * np.arange(0, CCDperTAP, mode)
                   for s, mode in zip(tap_starts, spectral_mode) if mode != 0]
    wavelengths = np.concatenate(wavelengths, axis=0)
    # TODO. This expression should be a warning, not a print; it doesn't indicate what needs to be indicated.
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

        The goal of this class is to create data attributes of the OCI SWIR, RED CCD, and BLUE CCD output.
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

           #. **swir_wavelengths**: The center wavelength of each swir CCD spectral band (nm) (Length 9).
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
           #. **ham_side**: 0 or 1 depending on which side of the HAM is in use.
           #. **scan_start_time**: The values read from scan_start_time (sec of day) (**N_s**).
           #. **scan_start_CCSDS_sec**: The values read from scan_start_CCSDS_sec (**N_s**).
           #. **scan_start_CCSDS_usec**: The values read from scan_start_CCSDS_usec (**N_s**).
           #. **scan_start_TIA**: Scan line start time in TAI (int64 in us) (**N_s**).
           #. **scan_start_UTC**: Scan line start time in UCT computed from scan_start_TIA (int64 in us) (**N_s**).
           #. **pixel_data_type_swir**: Data type of each SWIR pixel. Same shape as **Sci_swir_data**.
           #. **pixel_data_type_red**: Data type of each SWIR pixel. Same shape as **Sci_red_data**.
           #. **pixel_data_type_blue**: Data type of each SWIR pixel. Same shape as **Sci_blue_data**.

        In addition, **path** is defined, the path to a single input file or a list of paths to
        each input file when **append()** is used.

        When data files are accumulated by the append method, the first 6 data attributes in
        group 1 are tested for a match between the previously read file and the current one.
        Any mismatch will raise a ValueError indicating which was the first one to not match.
        Data attributes in Group 1 are all taken from the current file and are always valid.
        However, when appending, data attributes in Group 2 are accumulated in lists
        and don't take the form of numpy masked arrays until after the **append_conclude()**
        method has been called.

        numpy datatime64 objects representing the scan times can be generated by::

            TIA_datetime = np.datetime64('1958-01-01T00:00:00.000000') + self.scan_start_TIA
            midnight = np.datetime64('2020-10-30T00:00:00.000000')
            day_datetime = midnight + (self.oci_1.scan_start_time * 1e6).astype(np.int64)

        .. note:
            TIA can be converted to UST times using the files available by ftp from
            **https://hpiers.obspm.fr/iers/bul/bulc/ntp/**. In particular,
            **https://hpiers.obspm.fr/iers/bul/bulc/ntp/leap-seconds.list**. Also see
            **ftp://ftp.nist.gov/pub/time/leap-seconds.list**

        Typical use::

            import oci_l1a.read_oci_1a_class as oci
            path = "/path/to/my/data/"                  # Source folder
            files = oci.get_pace_oci_l1a_files(path)    # Get file paths
            oci_ccd = oci.READ_OCI_SPA_SPE_L1A(files[0])
            oci_ccd.subtract_dark()                     # Subtract the dark values
            oci_ccd.select_data_type(data_type=1)       # Trim to Earth View pixels only.
            print("Red Wavelengths")
            print(oci_ccd.red_wavelengths)
            print("Red Intensities")
            print(oci_ccd.Sci_red_data)

        .. note::
            The following methods permanently modify the data attributes in place.:
                #. **subtract_dark()**
                #. **select_data_type()**

            To preserve a version of the original use **copy.deepcopy()**. For example::

                from copy import deepcopy
                oci.append(read_oci_1a_class(file))
                oci2 = deepcopy(oci)
                oci2.subtract_dark()
                oci2 = deepcopy(oci)    # Reset oci2 to the original values.
    """
    def __init__(self, path=None, mask_dark=False):
        """Constructor for class READ_OCI_SPA_SPE_L1A"""

        self.path = path
        if path is None:                # Create empty instance for append.
            return
        with Dataset(path, 'r') as fid_cdf:
            dark_group = fid_cdf.groups['onboard_calibration_data']             # Read dark data
            self.DC_swir_data = dark_group.variables['DC_SWIR'][:]
            self.DC_red_data = dark_group.variables['DC_red'][:]
            self.DC_blue_data = dark_group.variables['DC_blue'][:]

            science_group = fid_cdf.groups['science_data']                      # Read the CCD raw data
            self.Sci_swir_data = science_group.variables['sci_SWIR'][()]
            self.Sci_red_data = science_group.variables['sci_red'][()]
            self.Sci_blue_data = science_group.variables['sci_blue'][()]

            scan_line_group = fid_cdf.groups['scan_line_attributes']            # Get scan line attributes
            self.ham_side = scan_line_group.variables['HAM_side'][()]
            self.scan_start_CCSDS_sec = scan_line_group.variables['scan_start_CCSDS_sec'][()]
            self.scan_start_CCSDS_usec = scan_line_group.variables['scan_start_CCSDS_usec'][()]
            self.scan_start_time = scan_line_group.variables['scan_start_time'][()]  # Seconds of day float64.

            # Need the spatial/spectral aggregation first so can calculate dark shift.
            spatial_spectral_group = fid_cdf.groups['spatial_spectral_modes']  # Get aggregation attributes
            self.aux_param_table = spatial_spectral_group.variables['aux_param_table'][()]
            self.blue_spectral_mode = spatial_spectral_group.variables['blue_spectral_mode'][()]
            self.red_spectral_mode = spatial_spectral_group.variables['red_spectral_mode'][()]
            self.spatial_aggregation = spatial_spectral_group.variables['spatial_aggregation'][()]
            self.spatial_zone_data_type = spatial_spectral_group.variables['spatial_zone_data_type'][()]
            self.spatial_zone_lines = spatial_spectral_group.variables['spatial_zone_lines'][()]

        # Compute the wavelengths from the spectral modes.
        self.swir_wavelengths, red_agg = getCCDbands3("", 'swir')
        self.red_wavelengths, red_agg = getCCDbands3(self.red_spectral_mode, 'red')
        self.blue_wavelengths, blue_agg = getCCDbands3(self.blue_spectral_mode, 'blue')
        # Map the Sci_?_data spatial (third index) to scan angle
        self.pixel_angle, self.pixel_data_type, spat_agg = get_CCD_pixel_type(self.spatial_zone_lines,
                                                                              self.spatial_aggregation,
                                                                              self.spatial_zone_data_type,
                                                                              self.Sci_red_data.shape[2])
        self.pixel_data_type_swir = np.broadcast_to(self.pixel_data_type, self.Sci_swir_data.shape)
        self.pixel_data_type_red = np.broadcast_to(self.pixel_data_type, self.Sci_red_data.shape)
        self.pixel_data_type_blue = np.broadcast_to(self.pixel_data_type, self.Sci_blue_data.shape)
        if mask_dark:
            self.Sci_swir_data = ma.masked_where(self.pixel_data_type_red != 9, self.Sci_swir_data)
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
        self.scan_start_UTC = self.scan_start_TIA - TIA2UTC * 1000000

    def append(self, other, droplast=False):
        """Append the CCD data from "other" to this object.

        :param other: A READ_OCI_SPA_SPE_L1A object
        :param droplast: If True, drop the last scan line before appending.
        :return: None

        :raises: ValueError if any of the spatial/spectral aggregation are different.

        It seems that in many files, the last scan line is either masked,
        contains zeros, or otherwise invalid values. Setting **droplast** to
        True, not the default, will exclude that last scan before appending arrays.

        Typical use::

            import oci_l1a.read_oci_1a_class as oci
            path = "/path/to/my/data/"                  # Source folder
            files = oci.get_pace_oci_l1a_files(path)
            oci_ccd = oci.READ_OCI_SPA_SPE_L1A()        # Create empty instance
            for file in files:                          # Loop through the files
                oci_ccd.append(oci.READ_OCI_SPA_SPE_L1A(file))  # and append
            oci_ccd.append_conclude()                   # Wrap up
            oci_ccd.subtract_dark()                     # Subtract the dark values
            oci_ccd.select_data_type(data_type=1)       # Trim to Earth View pixels only.
            print("Red Wavelengths")
            print(oci_ccd.red_wavelengths)
            print("Red Intensities")
            print(oci_ccd.Sci_red_data)
        """
        if self.path is None:
            self.path = []                              # Create all of the lists.
            self._DC_swir_data_lst, self._Sci_swir_data_lst = [], []
            self._DC_red_data_lst, self._Sci_red_data_lst = [], []
            self._DC_blue_data_lst, self._Sci_blue_data_lst = [], []
            self._ham_side_lst = []
            self._scan_start_time_lst, self._scan_start_CCSDS_sec_lst, self._scan_start_CCSDS_usec_lst = [], [], []
            # self.aux_param_table, self.blue_spectral_mode, self.red_spectral_mode = [], [], []
            # self.spatial_aggregation, self.spatial_zone_data_type, self.spatial_zone_lines = [], [], []
            self._scan_start_TIA_lst, self._scan_start_UTC_lst = [], []
        else:                                           # Test match between old and new parameters.
            if not np.array_equal(self.aux_param_table, other.aux_param_table):  #, equal_nan=True):
                raise ValueError("aux_param_table values are different.")
            if not np.array_equal(self.blue_spectral_mode, other.blue_spectral_mode):  #, equal_nan=True):
                raise ValueError("blue_spectral_mode values are different.")
            if not np.array_equal(self.red_spectral_mode, other.red_spectral_mode):  #, equal_nan=True):
                raise ValueError("red_spectral_mode values are different.")
            if not np.array_equal(self.spatial_aggregation, other.spatial_aggregation):  #, equal_nan=True):
                raise ValueError("spatial_aggregation values are different.")
            if not np.array_equal(self.spatial_zone_data_type, other.spatial_zone_data_type):  #, equal_nan=True):
                raise ValueError("spatial_zone_data_type values are different.")
            if not np.array_equal(self.spatial_zone_lines, other.spatial_zone_lines):  #, equal_nan=True):
                raise ValueError("spatial_zone_lines values are different.")

        # Assign those values of "other" that (should) be the same as "self" to self.
        self.path.append(other.path)
        self.aux_param_table = other.aux_param_table

        self.red_spectral_mode = other.red_spectral_mode
        self.blue_spectral_mode = other.blue_spectral_mode
        self.swir_wavelengths = other.swir_wavelengths
        self.red_wavelengths = other.red_wavelengths
        self.blue_wavelengths = other.blue_wavelengths

        self.spatial_aggregation = other.spatial_aggregation
        self.spatial_zone_data_type = other.spatial_zone_data_type
        self.spatial_zone_lines = other.spatial_zone_lines
        self.pixel_angle = other.pixel_angle
        self.pixel_data_type = other.pixel_data_type

        # Append each data element to the corresponding lists.
        if not droplast:
            self._DC_swir_data_lst.append(other.DC_swir_data)
            self._DC_red_data_lst.append(other.DC_red_data)
            self._DC_blue_data_lst.append(other.DC_blue_data)
            self._Sci_swir_data_lst.append(other.Sci_swir_data)
            self._Sci_red_data_lst.append(other.Sci_red_data)
            self._Sci_blue_data_lst.append(other.Sci_blue_data)
            self._ham_side_lst.append(other.ham_side)
            self._scan_start_time_lst.append(other.scan_start_time)
            self._scan_start_CCSDS_sec_lst.append(other.scan_start_CCSDS_sec)
            self._scan_start_CCSDS_usec_lst.append(other.scan_start_CCSDS_usec)
            self._scan_start_TIA_lst.append(other.scan_start_TIA)
            self._scan_start_UTC_lst.append(other.scan_start_UTC)
        else:
            self._DC_swir_data_lst.append(other.DC_swir_data[:-1, :, :])
            self._DC_red_data_lst.append(other.DC_red_data[:-1, :, :])
            self._DC_blue_data_lst.append(other.DC_blue_data[:-1, :, :])
            self._Sci_swir_data_lst.append(other.Sci_swir_data[:-1, :, :])
            self._Sci_red_data_lst.append(other.Sci_red_data[:-1, :, :])
            self._Sci_blue_data_lst.append(other.Sci_blue_data[:-1, :, :])
            self._ham_side_lst.append(other.ham_side[:-1])
            self._scan_start_time_lst.append(other.scan_start_time[:-1])
            self._scan_start_CCSDS_sec_lst.append(other.scan_start_CCSDS_sec[:-1])
            self._scan_start_CCSDS_usec_lst.append(other.scan_start_CCSDS_usec[:-1])
            self._scan_start_TIA_lst.append(other.scan_start_TIA[:-1])
            self._scan_start_UTC_lst.append(other.scan_start_UTC[:-1])

    def append_conclude(self):
        """Convert lists accumulated by **append** into masked arrays in group 2.

        :return: None

        Until this method is called, the data attributes in Group 2 are not valid.
        """
        # Concatenate the lists into arrays, map pixel_data_type to the CCD arrays,
        # and delete the lists that were used to accumulate to de-clutter the namespace.
        self.DC_swir_data = ma.concatenate(self._DC_swir_data_lst, axis=0)
        self.DC_red_data = ma.concatenate(self._DC_red_data_lst, axis=0)
        self.DC_blue_data = ma.concatenate(self._DC_blue_data_lst, axis=0)
        self.Sci_swir_data = ma.concatenate(self._Sci_swir_data_lst, axis=0)
        self.Sci_red_data = ma.concatenate(self._Sci_red_data_lst, axis=0)
        self.Sci_blue_data = ma.concatenate(self._Sci_blue_data_lst, axis=0)
        self.ham_side = ma.concatenate(self._ham_side_lst, axis=0)
        self.scan_start_time = ma.concatenate(self._scan_start_time_lst, axis=0)
        self.scan_start_CCSDS_sec = ma.concatenate(self._scan_start_CCSDS_sec_lst, axis=0)
        self.scan_start_CCSDS_usec = ma.concatenate(self._scan_start_CCSDS_usec_lst, axis=0)
        self.scan_start_TIA = ma.concatenate(self._scan_start_TIA_lst, axis=0)
        self.scan_start_UTC = ma.concatenate(self._scan_start_UTC_lst, axis=0)

        self.pixel_data_type_swir = np.broadcast_to(self.pixel_data_type,
                                                    self.Sci_red_data.shape)
        self.pixel_data_type_red = np.broadcast_to(self.pixel_data_type,
                                                   self.Sci_red_data.shape)
        self.pixel_data_type_blue = np.broadcast_to(self.pixel_data_type,
                                                    self.Sci_blue_data.shape)
        del self._DC_swir_data_lst,  self._DC_red_data_lst,  self._DC_blue_data_lst
        del self._Sci_swir_data_lst,  self._Sci_red_data_lst,  self._Sci_blue_data_lst
        del self._ham_side_lst, self._scan_start_time_lst
        del self._scan_start_TIA_lst, self._scan_start_UTC_lst
        del self._scan_start_CCSDS_sec_lst, self._scan_start_CCSDS_usec_lst

    def subtract_dark(self, start=40, dev=3):
        """Average the dark counts for each scan and subtract from each scan

        :param start: First index to use when computing the dark average (Default=40)
        :param dev: Multiple of sigma to mask when computing dark level (Default=3)
        :return:

        .. note::
            This method permanently modifies **Sci_swir_data**, **Sci_red_data**,
             and **Sci_blue_data** in memory.
        """
        # # For each dark array, compute the mean and std in the scan direction. Use these to mask outliers.
        # # Recompute the along scan mean. Subtract that value from all of the science data.
        # DC_swir_data = self.DC_swir_data[:, :, start:].astype(np.float64)
        # DC_swir_mean = DC_swir_data.mean(axis=2, keepdims=True)
        # DC_swir_std_ = DC_swir_data.std(axis=2, keepdims=True)
        # DC_swir_mean = ma.masked_where((DC_swir_data < DC_swir_mean - dev * DC_swir_std_)
        #                              | (DC_swir_data > DC_swir_mean + dev * DC_swir_std_),
        #                                 DC_swir_data).mean(axis=2, keepdims=True)
        # self.Sci_swir_data = self.Sci_swir_data.astype(np.float64) - DC_swir_mean
        #
        # DC_red_data = self.DC_red_data[:, :, start:].astype(np.float64)
        # DC_red_mean = DC_red_data.mean(axis=2, keepdims=True)
        # DC_red_std_ = DC_red_data.std(axis=2, keepdims=True)
        # DC_red_mean = ma.masked_where((DC_red_data < DC_red_mean - dev * DC_red_std_)
        #                             | (DC_red_data > DC_red_mean + dev * DC_red_std_),
        #                                DC_red_data).mean(axis=2, keepdims=True)
        # self.Sci_red_data = self.Sci_red_data.astype(np.float64) - DC_red_mean
        #
        # DC_blue_data = self.DC_blue_data[:, :, start:].astype(np.float64)
        # DC_blue_mean = DC_blue_data.mean(axis=2, keepdims=True)
        # DC_blue_std_ = DC_blue_data.std(axis=2, keepdims=True)
        # DC_blue_mean = ma.masked_where((DC_blue_data < DC_blue_mean - dev * DC_blue_std_)
        #                             | (DC_blue_data > DC_blue_mean + dev * DC_blue_std_),
        #                               DC_blue_data).mean(axis=2, keepdims=True)
        # self.Sci_blue_data = self.Sci_blue_data.astype(np.float64) - DC_blue_mean

        # For each dark array, compute the mean and std in the scan direction. Use these to mask outliers.
        # Recompute the along scan mean. Subtract that value from all of the science data.
        if 1:               # Old calculation
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
            # DC_red_mean = DC_red_mean.mean(axis=0, keepdims=True)   # ToDo. Temporary line. Remove
            # print(DC_red_mean.shape)
            self.Sci_red_data = self.Sci_red_data.astype(np.float64) - DC_red_mean

            DC_blue_data = self.DC_blue_data[:, :, start:].astype(np.float64)
            DC_blue_mean = DC_blue_data.mean(axis=2, keepdims=True)
            DC_blue_std_ = DC_blue_data.std(axis=2, keepdims=True)
            DC_blue_mean = ma.masked_where((DC_blue_data < DC_blue_mean - dev * DC_blue_std_)
                                        | (DC_blue_data > DC_blue_mean + dev * DC_blue_std_),
                                          DC_blue_data).mean(axis=2, keepdims=True)
            DC_blue_mean = DC_blue_mean.mean(axis=0, keepdims=True)
            # print(DC_blue_mean.shape)
            self.Sci_blue_data = self.Sci_blue_data.astype(np.float64) - DC_blue_mean

        else:
            DC_swir_mean = mask_outliers(self.DC_swir_data[:, :, start:].astype(np.float64)
                                         ).mean(axis=2, keepdims=True)
            self.Sci_swir_data = self.Sci_swir_data.astype(np.float64) - DC_swir_mean

            DC_red_mean = mask_outliers(self.DC_red_data[:, :, start:].astype(np.float64)
                                        ).mean(axis=2, keepdims=True)
            self.Sci_red_data = self.Sci_red_data.astype(np.float64) - DC_red_mean

            DC_blue_mean = mask_outliers(self.DC_blue_data[:, :, start:].astype(np.float64)
                                         ).mean(axis=2, keepdims=True)
            self.Sci_blue_data = self.Sci_blue_data.astype(np.float64) - DC_blue_mean

    def select_data_type(self, data_type=1):
        """Trim data attributes to those spatial pixels with the type **data_type**.

        :param data_type: Scalar data type to keep in the trimmed data. Default=1 Earth view.
        :return: None

        .. note::
            This method permanently modifies the data attributes in place.
        """
        index = self.pixel_data_type[0, 0, :] == data_type
        self.Sci_swir_data = self.Sci_swir_data[:, :, index]
        self.Sci_red_data = self.Sci_red_data[:, :, index]
        self.Sci_blue_data = self.Sci_blue_data[:, :, index]
        self.pixel_data_type = self.pixel_data_type[:, :, index]    # shape = (1, 1, n)
        self.pixel_angle = self.pixel_angle[index]
        self.pixel_data_type_swir = np.broadcast_to(self.pixel_data_type,
                                                    self.Sci_swir_data.shape)
        self.pixel_data_type_red = np.broadcast_to(self.pixel_data_type,
                                                   self.Sci_red_data.shape)
        self.pixel_data_type_blue = np.broadcast_to(self.pixel_data_type,
                                                    self.Sci_blue_data.shape)

    def plot_data(self, data_source, ax, vlim=None, func=np.mean):
        """Display images of spectral data averaged cross track and along track.

        :param data_source: Which image to plot: "blue", "dc_blue", "red", "dc_red", "swir", or "dc_swir"
        :param ax: list like of twp matplotlib.axes._subplots.AxesSubplot. First gets cross track average. Second gets along track average.
        :param vlim: If None (default) auto scale the color map. Otherwise vlim[0], vlim[1]
        :param func: The numpy function to plot. Default numpy.mean. But could be std, min, max, etc.
        :return: tuple of (image on axis 0, image on axis 1)

        The axes must be created by the calling routine. They can be part of a 2D grid.
        When called the first time, this method will place images on the axes, create
        the color bars, and place a text box for the maximum pixel value.
        If the given axes already have an image (or more), the first image on each
        will be updated with the new image data.

        Typical use::

            from matplotlib import pyplot as plt
            import oci_l1a.read_oci_1a_class as oci
            files = oci.get_pace_oci_l1a_files(path)    # Get file paths
            file = files[2]                             # Select one file
            figsize = (11.5, 5.75)                      # Define figure layout
            gridspec_kw = {"left" : 0.1,
                           "right" : 0.99,
                           "top" : 0.90,
                           "bottom" : 0.2,
                           "wspace" : 0.10,
                           "hspace" : 0.20,
                           "height_ratios" : [1, 2]
                           }
            fig, ax = plt.subplots(2, 2, sharex="col", sharey="row",
                figsize=figsize, gridspec_kw=gridspec_kw)
            fig.suptitle(file)
            ax[0, 0].set_title("Mean Computed Cross Track/Along Scan")
            ax[0, 1].set_title("Mean Computed Along Track/Cross Scan")
            ax[0, 0].set_ylabel("Spectral Band (nm)")
            ax[1, 0].set_ylabel("Spectral Band (nm)")
            ax[-1, 0].set_xlabel("Scan Start Time")
            ax[-1, 1].set_xlabel("Pixel Angle (deg)")
            ax[-1, 0].tick_params(axis='x', labelrotation=90)
            ax[-1, 1].tick_params(axis='x', labelrotation=90)

            oci_1 = oci.READ_OCI_SPA_SPE_L1A(file)      # Read a file data
            oci_1.select_data_type(1)                   # Restrict to Earth View
            oci_1.subtract_dark(start=0)                # Dark subtraction.
            try:                # This try is more important when looping
                image0 = oci_1.plot_data("swir", ax=ax[0, :], )
                image1 = oci_1.plot_data("red", ax=ax[1, :])
            except Exception as e:
                print("Failed to plot item {:}, file {:}. "
                    "Skipping to next.".format(i, file))
                raise
                continue        # Use when in a loop instead of raise
            plt.show()
        """

        if data_source == "blue":
            ev = self.Sci_blue_data
            y_ticks = np.linspace(0, ev.shape[1], 16, endpoint=False, dtype=np.int32)
            wavelengths = self.blue_wavelengths[y_ticks]
        elif data_source == "dc_blue":
            ev = self.DC_blue_data
            y_ticks = np.linspace(0, ev.shape[1], 16, endpoint=False, dtype=np.int32)
            wavelengths = self.blue_wavelengths[y_ticks]
        elif data_source == "red":
            ev = self.Sci_red_data
            y_ticks = np.linspace(0, ev.shape[1], 16, endpoint=False, dtype=np.int32)
            wavelengths = self.red_wavelengths[y_ticks]
        elif data_source == "dc_red":
            ev = self.DC_red_data
            y_ticks = np.linspace(0, ev.shape[1], 16, endpoint=False, dtype=np.int32)
            wavelengths = self.red_wavelengths[y_ticks]
        elif data_source == "swir":
            ev = self.Sci_swir_data
            y_ticks = range(9)
            wavelengths = self.swir_wavelengths
        elif data_source == "dc_swir":
            ev = self.DC_swir_data
            y_ticks = range(9)
            wavelengths = self.swir_wavelengths
        else:
            raise ValueError("data_source name {:} not supported.".format(data_source))

        ev_ct_mean = func(ev, axis=2).T
        ev_at_mean = func(ev, axis=0)

        if (vlim == 'auto') or vlim is None:  # Set color range to data range
            vlim = ((ev_ct_mean.min(), ev_ct_mean.max()),
                    (ev_at_mean.min(), ev_at_mean.max()))
        print("{:6.3f}  {:6.3f}    {:6.3f}  {:6.3f}"
              .format(ev_ct_mean.min(), ev_ct_mean.max(),
                      ev_at_mean.min(), ev_at_mean.max()))
        fig = ax[0].get_figure()

        images = ax[0].get_images(), ax[1].get_images()
        print(images)
        if (len(images[0]) == 0) or (len(images[1]) == 0):
            print("Creating images")
            # Create the two images. Add colorbars and the text box (Max Pixel) to the axes.
            # gid='max' is to find it again.
            images = (ax[0].imshow(ev_ct_mean, aspect='auto', origin='lower', cmap="gray"),
                      ax[1].imshow(ev_at_mean, aspect='auto', origin='lower', cmap="gray"))
            fig.colorbar(images[0], ax=ax[0])
            fig.colorbar(images[1], ax=ax[1])
            ax[0].text(1, 1, "", transform=ax[0].transAxes, color='r', ha='left', gid="max0")
            ax[1].text(1, 1, "", transform=ax[1].transAxes, color='r', ha='left', gid='max1')

            # By default, the rightmost tick is outside the range of scan times. So
            # limit the x range so that the right most tick label is valid.
            ax[0].set_xlim(0, ev_ct_mean.shape[1] - 2)
            ax[1].set_xlim(0, ev_at_mean.shape[1] - 2)
        else:
            # Update the images that are there; assumed the first image on each axis
            images = (images[0][0], images[1][0])
            im0_size, im1_size = images[0].get_size(), images[1].get_size()

            if im0_size != ev_ct_mean.shape:        # Pad the new image if smaller than current
                data = np.ma.zeros(im0_size, dtype=ev_ct_mean.dtype)
                data[:] = np.nan
                try:
                    data[:ev_ct_mean.shape[0], :ev_ct_mean.shape[1]] = ev_ct_mean
                except ValueError:
                    raise ValueError("New array is larger than original image.") from None
                images[0].set_data(data)            # Update the image
            else:
                images[0].set_data(ev_ct_mean)      # Otherwise just update the image

            if im1_size != ev_at_mean.shape:    # Pad the new image if smaller than current
                data = np.ma.zeros(im1_size, dtype=ev_at_mean.dtype)
                data[:] = np.nan
                try:
                    data[:ev_at_mean.shape[0], :ev_at_mean.shape[1]] = ev_at_mean
                except ValueError:
                    raise ValueError("New array is larger than original image.") from None
                images[1].set_data(data)            # Update the image
            else:
                images[1].set_data(ev_at_mean)      # Otherwise just update the image

        # Now that the axes have an (updated) image, update the supporting information.
        print(images, vlim)
        images[0].set_clim(vlim[0])
        images[1].set_clim(vlim[1])
        ax[0].findobj(lambda x: x.get_gid() == 'max0')[0].set_text("{:.1f}".format(ev_ct_mean.max()))
        ax[1].findobj(lambda x: x.get_gid() == 'max1')[0].set_text("{:.1f}".format(ev_at_mean.max()))

        # Convert all of the line start times to strings.
        # Todo. Make these times "correct"
        scan_start_times = [
            (datetime.datetime(2020, 12, 8, 0, 0, 0) + datetime.timedelta(seconds=x)).strftime("%H:%M:%S")
            for x in self.scan_start_time.filled(fill_value=0)]
        ax[0].xaxis.set_major_formatter(mpl.ticker.IndexFormatter(scan_start_times))
        ax[0].xaxis.set_major_locator(mpl.ticker.LinearLocator())  # 11 ticks between view limits

        pixel_angle = ["{:.1f}".format(x) for x in (self.pixel_angle - PPR2NADIR)]
        ax[1].xaxis.set_major_formatter(mpl.ticker.IndexFormatter(pixel_angle))

        ax[0].set_yticks(y_ticks)
        ax[0].set_yticklabels(["{:.0f}".format(x) for x in wavelengths])

        return images

    @staticmethod
    def mask_outliers(data, frac=3, axis=2):
        """From a 3D OCI image array, return an array masking values outside frac stds.

        :param data: A 3D OCI image array (DC_swir_data, Sci_red_data, etc)
        :param frac: Number standard deviations to not mask (Default=3)
        :param axis: Axis along which to do mean and std calculation. (Default=2)
        :return: Masked array result

        By default, means and standard deviations are computed along the third axis, the cross track
        scan direction. Any pixels outside frac * std from the mean are masked.
        Setting axis=0 will compute mean and std in the along track/flight direction.

        Typical use::

            import oci_l1a.read_oci_1a_class as oci
            path = "/path/to/my/data/"                  # Source folder
            files = oci.get_pace_oci_l1a_files(path)
            oci_ccd = oci.READ_OCI_SPA_SPE_L1A(file[0]))
            Sci_swir_data = oci_ccd.mask_outliers(oci_ccd.Sci_swir_data)
            Sci_red_data = oci_ccd.mask_outliers(oci_ccd.Sci_red_data)
        """
        mean = data.mean(axis=axis, keepdims=True)
        std = data.std(axis=axis, keepdims=True)
        upper = mean + frac * std
        lower = mean - frac * std
        return np.ma.masked_where((data < lower) | (data > upper), data)


class READ_OCI_THERMAL_L1A(object):
    """Class to read the OCI l1a thermal data.

        :param path: File path, either a Path object or string. Default (None)

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
        each input file when **append()** is used. Note that temperatures are in deg C (I believe.)

        Typical use::

            import oci_l1a.read_oci_1a_class as oci
            path = "/path/to/my/data/"                  # Source folder
            files = oci.get_pace_oci_l1a_files(path)    # Get file paths
            oci_therm = oci.READ_OCI_THERMAL_L1A(files[0])
            print("aob_temperatures")
            print(oci_therm.aob_temperatures)
            print("red_fpa_temperatures")
            print(oci_therm.red_fpa_temperatures)
    """

    def __init__(self, path=None):
        """Constructor for class READ_OCI_THERMAL_L1A"""

        self.path = path
        if path is None:
            return
        with Dataset(path, 'r') as fid_cdf:
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
            # icdu_cnts = engineering_data.variables['ICDU_thermistors'][()]
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

            import oci_l1a.read_oci_1a_class as oci
            path = "/path/to/my/data/"                  # Source folder
            files = oci.get_pace_oci_l1a_files(path)    # Get file paths

            oci_therm = oci.READ_OCI_THERMAL_L1A()
            for file in files:
                oci_therm.append(READ_OCI_THERMAL_L1A(file))
            oci_therm.append_conclude()
            print("aob_temperatures")
            print(oci_therm.aob_temperatures)
            print("red_fpa_temperatures")
            print(oci_therm.red_fpa_temperatures)
       """
        if self.path is None:
            self.path = []
            self._aob_temperatures_lst, self._blue_fpa_temperatures_lst = [], []
            self._dau_temperatures_lst, self._dau_tlm_time_lst = [], []
            self._icdu_temperatures_lst, self._red_fpa_temperatures_lst = [], []
            self._sds_det_temperatures_lst, self._tc_tlm_time_lst = [], []

        self.path.append(other.path)
        self._aob_temperatures_lst.append(other.aob_temperatures)
        self._blue_fpa_temperatures_lst.append(other.blue_fpa_temperatures)
        self._dau_temperatures_lst.append(other.dau_temperatures)
        self._dau_tlm_time_lst.append(other.dau_tlm_time)
        # self._icdu_temperatures_lst.append(other.icdu_temperatures)
        self._red_fpa_temperatures_lst.append(other.red_fpa_temperatures)
        self._sds_det_temperatures_lst.append(other.sds_det_temperatures)
        self._tc_tlm_time_lst.append(other.tc_tlm_time)

    def append_conclude(self):
        """Convert lists accumulated by **append** into masked arrays in group 2.

        :return: None

        Until this method is called, the thermal numpy arrays are not defined.
        """
        # Concatenate the lists into arrays.
        self.aob_temperatures = ma.concatenate(self._aob_temperatures_lst, axis=0)
        self.blue_fpa_temperatures = ma.concatenate(self._blue_fpa_temperatures_lst, axis=0)
        self.dau_temperatures = ma.concatenate(self._dau_temperatures_lst, axis=0)
        self.dau_tlm_time = ma.concatenate(self._dau_tlm_time_lst, axis=0)
        # self.icdu_temperatures = ma.concatenate(self._icdu_temperatures_lst, axis=0)
        self.red_fpa_temperatures = ma.concatenate(self._red_fpa_temperatures_lst, axis=0)
        self.sds_det_temperatures = ma.concatenate(self._sds_det_temperatures_lst, axis=0)
        self.tc_tlm_time = ma.concatenate(self._tc_tlm_time_lst, axis=0)

        del self._aob_temperatures_lst,  self._blue_fpa_temperatures_lst
        del self._dau_temperatures_lst
        del self._dau_tlm_time_lst,  self._red_fpa_temperatures_lst
        del self._sds_det_temperatures_lst, self._tc_tlm_time_lst

    # def plot(self, variable, ax=None, plt_type='3d'):
    #     """Method doesn't work. Just return"""
    #     # return
    #
    #     """Create a 3d plot of the specified thermal data.
    #
    #     :param variable: string name of the variable to plot
    #     :param ax: If not None, 3d axis on which to place the plot. If None (default) make a new figure.
    #     :param plt_type: Either '3d' for 'imag' to specify the type of display. (Default '3d')
    #     :return:
    #
    #     Valid string values for **variable** are:
    #        - "aob_temperatures"
    #        - "blue_fpa_temperatures"
    #        - "dau_temperatures"
    #        - "red_fpa_temperatures"
    #        - "sds_det_temperatures"
    #
    #     If an axis is provided, it must have been created with projection='3d'
    #     for a 3d plot and must be created without a projection for an image.
    #     """
    #     if plt_type not in ['3d', 'imag']:
    #         raise ValueError("Plot type must be either '3d' or 'imag'. {:} is neither.".format(plt_type))
    #
    #     if variable == "aob_temperatures":
    #         title = "AOB_temperatures"
    #         time = self.dau_tlm_time
    #         labels = therm.DAU_AOB_temps["name"]
    #         temps = self.aob_temperatures.fill(np.nan)
    #     elif variable == "blue_fpa_temperatures":
    #         title = "blue_fpa_temperatures"
    #         time = self.dau_tlm_time
    #         labels = therm.DAU_FPA_temps["name"]
    #         temps = self.blue_fpa_temperatures.fill(np.nan)
    #     elif variable == "dau_temperatures":
    #         title = "dau_temperatures"
    #         time = self.dau_tlm_time
    #         labels = therm.DAU_temps["name"]
    #         temps = self.dau_temperatures.fill(np.nan)
    #     elif variable == "red_fpa_temperatures":
    #         title = "red_fpa_temperatures"
    #         time = self.dau_tlm_time
    #         labels = therm.DAU_FPA_temps["name"]
    #         temps = self.red_fpa_temperatures.fill(np.nan)
    #     elif variable == "sds_det_temperatures":
    #         title = "sds_det_temperatures"
    #         time = self.dau_tlm_time
    #         labels = therm.DAU_SDS_det_temps["name"]
    #         temps = self.sds_det_temperatures.fill(np.nan)
    #     else:
    #         raise ValueError("Variable given {:} not recognized.".format(variable))
    #
    #     t_min, t_max = np.nanmin(temps), np.nanmax(temps)
    #     # If needed, create the figure and axis in the proper form.
    #     subplot_kw = {"projection": '3d'} if plt_type == '3d' else {}
    #     if ax is None:
    #         fig, ax = plt.subplots(1, 1, subplot_kw=subplot_kw)
    #     else:
    #         fig = ax.figure
    #     ax.set_title(title)
    #     x_ticks = range(len(labels))        # Tick and label each channel
    #     ax.set_xticks(x_ticks)
    #     ax.set_xticklabels(labels)
    #     fig.autofmt_xdate(rotation=90)
    #     if plt_type == 'imag':              # Make the plot
    #         surf = ax.imshow(temps)
    #         ax.set_aspect("auto")
    #     else:
    #         l, t = np.meshgrid(range(len(time)), x_ticks)
    #         surf = ax.plot_surface(t, l, temps.T,
    #                                vmin=t_min, vmax=t_max,
    #                                cmap=mpl.rcParams['image.cmap'])
    #         ax.set_aspect("auto")
    #     fig.colorbar(surf, ax=ax)


def main(args):
    """Read the PACE OCI L1A netCDF files in args.paths and plot temperatures and spectrum.

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
        thermal.append(READ_OCI_THERMAL_L1A(file))
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
    oci2.select_data_type(9)           # Clip the spatial pixels to zone 9
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
        description='Test of class A class to read the spatial/spectral data from an '
                    'OCI l1a file and return it with spectral and pixel lables.',
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