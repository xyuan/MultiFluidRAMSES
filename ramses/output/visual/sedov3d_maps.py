from numpy import array
import pylab
from pymses.analysis.visualization import *
from pymses import RamsesOutput
from pymses.utils import constants as C
# Ramses data
ioutput = 4
ro = RamsesOutput("Data/Sedov3d/output", ioutput)
amr = ro.amr_source(["rho", "P"])
# Map region
## center = [ 0.567811, 0.586055, 0.559156 ]
center = [ 0.5, 0.5, 0.5 ]
# Map processin#g
#cam = Camera()
cam = Camera(center=center, line_of_sight_axis='z', up_vector="y", region_size=(1.0, 1.0), far_cut_depth=0.2, map_max_size=64, distance=0., log_sensitive=True)
from pymses.analysis.visualization import SliceMap, ScalarOperator
op = ScalarOperator(lambda dset: dset["rho"])
map = SliceMap(amr, cam, op, z=0.3) # create a density slice map at z=0.4 depth position
  
#factor = ro.info["unit_density"].express(C.H_cc)
#scale = ro.info["unit_length"].express(C.kpc)
# imshow(map)
pylab.imshow(map)
# Save map into HDF5 file
#mapname = "sedov3d%5.5i"%ioutput
#h5fname = save_map_HDF5(map, cam, map_name=mapname)
# Plot map into Matplotlib figure/PIL Image
#fig = save_HDF5_to_plot(h5fname, map_unit=("H/cc",factor), axis_unit=("kpc", scale), cmap="jet")
# pil_img = save_HDF5_to_img(h5fname, cmap="jet")
# Save into PNG image file
# save_HDF5_to_plot(h5fname, map_unit=("H/cc",factor), axis_unit=("Mpc", scale), img_path="./", # save_HDF5_to_img(h5fname, img_path="./", cmap="jet")


# pylab.show()
