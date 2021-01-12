This Python package contains methods to read PACE Ocean Color Instrument L1A thermal and spatial/spectral
data from the netCDF L1A files.


Installation
++++++++++++

1. Activate your Python 3 conda environment.

2. cd to the location you wish to store the package.

3. Type the following

    git clone git@gs490v-gitlab.ndc.nasa.gov:497_OCI/oci_l1a_in_python.git

4. cd to **oci_l1a_in_python**.

5. To install the package as a permanent part of your system,

   pip install .

   If you later wish to remove it

   pip uninstall oci_l1a


6. The document **ocil1a.pdf** contains complete documentation on using the package.

7. Building the HTML documentation requires **Sphinx**. To build the HTML documentation,
   cd to **docs** and type **make html**. The file **docs/_build/html/index.html** is the
   documentation tree starting point. A latex file of the documentation is made by typing
   **make latex**. The latex file needs to be processed to a pdf file.
