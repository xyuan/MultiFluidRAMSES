### ipython --pylab

import pymses
import pylab
import numpy

ro = pymses.RamsesOutput("Data/Sedov3d/output", 5)
amr = ro.amr_source(["rho"])

#sph_center = [0.5, 0.5, 0.5]
#sph_radius = 0.5
#from pymses.utils.regions import Sphere
#sph = Sphere(sph_center, sph_radius)

#from pymses.analysis import sample_points
#points = sph.random_points(1.0e6) # 1M sampling points

end = [1,1,1]
N = 2 * 2**ro.info["levelmax"]
points = 1. * numpy.array([end]*N)
for i in range(N): points[i] *= (i+1.)/N
point_dset = sample_points(amr, points)
pylab.plot(point_dset.points,point_dset['rho'])

#rho_weight_func = lambda dset: dset["rho"]
#r_bins = numpy.linspace(0.0, sph_radius, 200)
#from pymses.analysis import bin_spherical
#rho_profile = bin_spherical(point_dset, sph_center, rho_weight_func, r_bins, divide_by_counts=True)

#pylab.plot(r_bins[0:len(r_bins)-1],rho_profile)

#from pymses.filters import CellsToPoints
#cell_source = CellsToPoints(amr)
#cells = cell_source.flatten()
#ccenters = cells.points


##import numpy as np
##import matplotlib.pyplot as plt

##plt.plot(rho_profile)
