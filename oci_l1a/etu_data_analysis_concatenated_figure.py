# To make a combined/concatenated figure....

from matplotlib import pyplot as plt
import oci_l1a.read_oci_1a_class as oci
path = "C:/Users/Vanlossing/Desktop/OBPG_DATA_ANALYSIS/DATA/1/"                  # Source folder
from pprint import pprint
from pathlib import Path

FIG_PATH = Path("C:/Users/Vanlossing/Desktop/OBPG_DATA_ANALYSIS/DATA/Plot_Outputs/")   # Change for your use.
files = oci.get_pace_oci_l1a_files(path)
pprint(files)

#files = oci.get_pace_oci_l1a_files(path)    # Get file paths
file = files[0]              # Select one file
figsize = (11.5, 5.75)       # Define figure layout
gridspec_kw = {"left" :1,
               "right" :1,
               "top" :1,
               "bottom" :1,
               "wspace" :1,
               "hspace" :1,
               "height_ratios" : [1,1]
               }
fig, ax = plt.subplots(2, 2, sharex="col", sharey="row")
ax[0, 0].set_title("Mean Computed Cross Track/Along Scan")
ax[0, 1].set_title("Mean Computed Along Track/Cross Scan")
ax[0, 0].set_ylabel("Spectral Band (nm)")
ax[1, 0].set_ylabel("Spectral Band (nm)")
ax[-1, 0].set_xlabel("Scan Start Time")
ax[-1, 1].set_xlabel("Pixel Angle (deg)")
ax[-1, 0].tick_params(axis='x', labelrotation=90)
ax[-1, 1].tick_params(axis='x', labelrotation=90)
oci_1 = oci.READ_OCI_SPA_SPE_L1A()

oci_1 = oci.READ_OCI_SPA_SPE_L1A()
for i, file in enumerate(files):
    print(i, file)
    oci_1.append(oci.READ_OCI_SPA_SPE_L1A(file), droplast=True)
oci_1.append_conclude()

oci_1.select_data_type(1)
oci_1.subtract_dark(start=0)
image0 = oci_1.plot_data("swir", ax=ax[0, :])
image1 = oci_1.plot_data("red", ax=ax[1, :])

fig.canvas.draw()
path_fig = FIG_PATH / "stray_light_auto_concatenated_2{:03d}.png".format(i)
fig.savefig(path_fig)
plt.show()

