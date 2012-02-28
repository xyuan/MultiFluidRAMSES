import numpy
import matplotlib
import pylab
import pymses
from pymses.analysis import sample_points
from pymses.analysis.visualization import *
from pymses.sources.ramses.octree import RamsesOctreeReader as ROR
ROR.fields_by_file = {"hydro" : [ROR.Scalar("d" ,  0), \
                                 ROR.Scalar("ux" , 1), \
                                 ROR.Scalar("uy" , 2), \
                                 ROR.Scalar("uz" , 3), \
                                 ROR.Scalar("P" ,  4), \
                                 ROR.Scalar("f" ,  5), \
                                 ROR.Scalar("tS",  6), \
                                 ROR.Scalar("tI",  7), \
                                 ROR.Scalar("tB",  8), \
                                 ROR.Scalar("g" ,  9) ],
}

legend = {"r" : 'radius [normalized]', \
          "d" : 'density [mp/cm3]', \
          "ux": 'x-velocity [km/s]', \
          "uy": 'y-velocity [km/s]', \
          "uz": 'z-velocity [km/s]', \
          "P" : 'pressure [1e12 erg]', \
          "f" : 'ejecta fraction', \
          "tS": 'time since shocked [yr]', \
          "tI": 'ionization time-scale [yr]', \
          "tB": 'radiative losses time-scale [yr]', \
          "g" : 'gamma'}
color_palette = pylab.cm.cool
color_list = {"d": 'r', "ux": 'g', "uy": 'g', "uz": 'g', "P": 'b', "f": 'k', "tS": 'y', "tI": 'c', "tB": 'm', "g": 'k'}

basedir = "/Users/gferrand/data/Tycho_em2/"

def connect(dir="L68_off",out=11,var=["d","ux","uy","uz","P","f","tS","tI","tB","g"]):
    global output, grid
    output = pymses.RamsesOutput(basedir+dir, out)
    print output.info
    grid = output.amr_source(var)
    
def load_profiles(end=[1,1,1],N=None):
    global data
    # define points
    if N==None: N = 2 * 2**output.info["levelmax"]
    points = 1. * numpy.array([end]*N)
    for i in range(N): points[i] *= (i+1.)/N
    # sample data
    data = pymses.analysis.sample_points(grid, points)
    print data.npoints, 'points'
    return data
    
def show_profile(var='d',func=None,log=True,col=None):
    # pick variable
    if func: data_y = func(data)
    else: data_y = data[var]
    # pick colour
    if col == None: color = color_list[var]
    elif col <= 0: color = color_palette(abs(col))
    else: color = col
    # pick axes
    if pylab.gcf().axes != []: 
        plot = pylab.gcf().axes[0]
    else: 
        plot = pylab.axes()
        plot.set_xlabel(legend["r"])
        if not func: plot.set_ylabel(legend[var])
    # plot
    if log: plot.semilogy(data.points,data_y,color=color)
    else:   plot.plot(    data.points,data_y,color=color)
    pylab.gcf().canvas.draw()
    return data_y
    
def loop_profile(dir="L68_off",out=[11],end=[1.,1.,1.],N=None,var='d',log=True,fig=1,over=False):
    # set figure
    if not over:
        pylab.close(fig)
        figure = matplotlib.pyplot.figure(num=fig)
    # loop over outputs
    for i in out:
        print 'output #',i
        connect(dir=dir,out=i,var=[var])
        load_profiles(end=end,N=N)
        if len(out)>1: col = -(i-out[0])/(1.*out[-1]-out[0])
        else: col = 0
        show_profile(var=var,log=log,col=col)
    # colour scale
    if not over:
        bar_plot = figure.add_axes([0.92,0.1,0.015,0.8])
        norm = matplotlib.colors.Normalize(vmin=out[0], vmax=out[-1])
        color_bar = matplotlib.colorbar.ColorbarBase(bar_plot, cmap=color_palette, norm=norm, orientation='vertical')
        color_bar.set_label('output')

def slice(var='d',z=0,log=True):
    cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis='z', region_size=[1., 1.], log_sensitive=log)
    op = ScalarOperator(lambda dset: dset[var])
    map = SliceMap(grid, cam, op, z=z)
    return map

def load_slices(z=0):
    global data
    camera = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis='z', region_size=[1., 1.])
    points = camera.get_slice_points(z)
    data = sample_points(grid, points)
    return data

def show_slice(var="d",log=False):
    N = numpy.sqrt(data.npoints)
    map = data[var].reshape(N,N)
    if log: map = numpy.log(map)
    pylab.imshow(map,origin='lower')
    return map

