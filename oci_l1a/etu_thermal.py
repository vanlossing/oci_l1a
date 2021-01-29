#!/usr/bin/env python
"""This file contains (some of) the conversations of OCI L1A thermal data to deg C.

Attempts to describe conversations for
L1A engineering_data: AOB_temperatures, blue_FPA_temperatures, ICDU_thermisters, red_FPA_temperatures,
and SDS_det_temperatures.

The temperature conversions are taken from
**OCI-SYS-TN-0270_OCI Temperature Sensor List for K3 Temperature Calibration Revision G**,
source file, **pace_oci_thermistor_table.txt** from Gene Eplee (21 October 2020,
and code samples **read_pace_oci_icdu_temps.pro** and **read_pace_full_oci_temps.pro**
by Gene Eplee. Note that these files don't currently include the **Vmeas, gnd** correction shown
for some channels in the Excel file **OCI-SYS-TN-0270_OCI**.

Note that the **engineering_data** group contains three time variables: **DAU_tlm_time**,
**DDC_tlm_time**, and **TC_tlm_time**. It is not clear to me which ones go with which data sets.

   - uint16 **blue_FPA_temperatures(tlm_packets, DAU_FPA_temps)**
     'DAU_FPA_temps': <class 'netCDF4._netCDF4.Dimension'>: name = 'DAU_FPA_temps', size = 14,
   - uint16 **red_FPA_temperatures(tlm_packets, DAU_FPA_temps)**
   - uint16 **SDS_det_temperatures(tlm_packets, DAU_SDS_det_temps)**
     'DAU_SDS_det_temps': <class 'netCDF4._netCDF4.Dimension'>: name = 'DAU_SDS_det_temps', size = 16,
   - uint16 **AOB_temperatures(tlm_packets, DAU_AOB_temps)**
     'DAU_AOB_temps': <class 'netCDF4._netCDF4.Dimension'>: name = 'DAU_AOB_temps', size = 9,
   - uint16 **DAU_temperatures(tlm_packets, DAU_temps)**
     'DAU_temps': <class 'netCDF4._netCDF4.Dimension'>: name = 'DAU_temps', size = 14,

Some references
https://www.earthinversion.com/utilities/reading-NetCDF4-data-in-python/

https://www.unidata.ucar.edu/software/netcdf/documentation/NUG/_best_practices.html#bp_Conventions

ToDo. These conversions have not been checked that they generate correct temperatures.
"""

__author__ = "Christopher Field"
__email__ = "Christopher.T.Field@NASA.gov"
__date_assembled__ = "26 October 2020"

DAU_FPA_temps = {"name" : [0,
                           1,
                           2,
                           3,
                           4,
                           5,
                           6,
                           "FPA CCD Right",
                           "FPA CCD Left",
                           "FPA PreAmp Right",
                           "FPA PreAmp Left",
                           "FPA ADC Low",
                           "FPA ADC High",
                           "FPA Clock Driver"],
                 "convert" : 14 * [[255.86, -255.74]]
                 }
"""Dictionary of names (not yet defined) and conversions for dimension **DAU_FPA_temps**.

Used by the variables **blue_FPA_temperatures** and **red_FPA_temperatures**. 
This is a dictionary with 
key **"name"** a list of sensor names and 
key **"convert"** a list of lists of polynomial coefficients from highest to lowest power
to convert volts to deg C (after conversion from raw counts to volts.."""

DAU_SDS_det_temps = {"name" : ["SDS Det. Temp. 1",
                               "SDS Det. Temp. 2",
                               "SDS Det. Temp. 3",
                               "SDS Det. Temp. 4",
                               "SDS Det. Temp. 5",
                               "SDS Det. Temp. 6",
                               "SDS Det. Temp. 7",
                               "SDS Det. Temp. 8",
                               "SDS Det. Temp. 9",
                               "SDS Det. Temp. 10",
                               "SDS Det. Temp. 11",
                               "SDS Det. Temp. 12",
                               "SDS Det. Temp. 13",
                               "SDS Det. Temp. 14",
                               "SDS Det. Temp. 15",
                               "SDS Det. Temp. 16"],
                 "convert" : 8 * [[0.245522553, 3.0735486, 35.69365050, -13.0848956]] +
                                 [[0.302496421, 3.35758291, 35.60259635, -11.61394172]] +
                                 [[0.301284786, 3.34901865, 35.51833920, -12.04654415]] +
                                 [[0.302300868, 3.35297900, 35.56135379, -11.78251807]] +
                                 [[0.301643744, 3.34948270, 35.52123620, -12.10587313]] +
                                 [[0.303496706, 3.36037806, 35.52500658, -12.11985400]] +
                                 [[0.301202923, 3.35492981, 35.57651802, -11.73430196]] +
                                 [[0.245522553, 3.07354860, 35.69365050, -13.0848956]] +
                                 [[0.300985803, 3.35269595, 35.55958500, -11.82720046]]
                 }
"""Dictionary of names and conversions for dimension **DAU_SDS_det_temps**.

Used by the variable **SDS_det_temperatures**. 
This is a dictionary with 
key **"name"** a list of sensor names and 
key **"convert"** a list of list of polynomial coefficients from highest to lowest power
to convert volts to deg C (after conversion from raw counts to volts.."""

DAU_AOB_temps = {"name" : ["AOB Temp. 1",
                           "AOB Temp. 2",
                           "AOB Temp. 3",
                           "AOB Temp. 4",
                           "AOB Temp. 5",
                           "AOB Temp. 6",
                           "AOB Temp. 7",
                           "AOB Temp. 8",
                           "AOB Temp. 9"],
                 "convert" : 4 * [[0.245522553, 3.073548600, 35.69365050, -13.0848956]] +
                                 [[0.301499448, 3.343719770, 35.47645744, -12.39969832]] +
                                 [[0.302564649, 3.347980073, 35.48184146, -12.36104956]] +
                             3 * [[0.245522553, 3.073548600, 35.69365050, -13.0848956]]
                 }
"""Dictionary of names and conversions for dimension **DAU_AOB_temps**.

Used by the variable **AOB_temperatures**. 
This is a dictionary with 
key **"name"** a list of sensor names and 
key **"convert"** a list of list of polynomial coefficients from highest to lowest power
to convert volts to deg C (after conversion from raw counts to volts.."""


DAU_temps = {"name" : ["DDC Temp 1",
                       "DDC Temp 2",
                       "DDC FPGA",
                       "DAUC Board Temp",
                       "Red Board Temp 1",
                       "Red Board Temp 2",
                       "Red Board Temp 3",
                       "Red Board Temp 4",
                       "DDSA P3P3 Temp",
                       "DDSA FW Temp.",
                       "DDSA PN12P0 Temp.",
                       "Red BDSA P3P3 Temp",
                       "Red BDSA FW Temp",
                       "Red BDSA PN12P0 Temp"],
                 "convert" : [[ 0.0, -54.531611,   290.333968,  -606.11493,  622.68509, -365.375416, 141.114874]] +
                             [[ 0.0, -54.531611,   290.333968,  -606.11493,  622.68509, -365.375416, 141.114874]] +
                             [[ 0.0, 541.099116, 26278.28278, -84007.615787, 97243.367511, -49956.383069, 9812.997496]] +
                             [[ 0.0, -19.605129,   130.938995,  -340.447308,  431.156262, -305.048600, 135.874242]] +
                            4 *  [[9.4769, -83.608, 290.03, -517.6, 514.28,  -317.43, 140.66]] +
                            6 * [[0.0, -31.507870, 206.54221, -513.448165, 601.66252, -367.14378, 122.315344]]
                 }
"""Dictionary of names and conversions for dimension **DAU_temps**.

Used by the variable **DAU_temperatures**. 
This is a dictionary with 
key **"name"** a list of sensor names and 
key **"convert"** a list of list of polynomial coefficients from highest to lowest power
to convert volts to deg C (after conversion from raw counts to volts.."""

thermal_dims = {"DAU_FPA_temps" : DAU_FPA_temps,
                "DAU_SDS_det_temps" : DAU_SDS_det_temps,
                "DAU_AOB_temps" : DAU_AOB_temps,
                "DAU_temps" : DAU_temps}
"""Dictionary of thermal channel dimensions. 
Each dimension previously defined in this module can be accessed by
a string key of its name permitting looking it up from the netCDF 
variable data structure."""