import numpy
import matplotlib
import matplotlib.animation as anim
import pylab
import pymses
import os, sys

import re

from matplotlib import rc
from pymses.utils import constants as C
from pymses.analysis import sample_points
from pymses.analysis.visualization import *
from pymses.sources.ramses.octree import RamsesOctreeReader as ROR


#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 20}
font = {'size'   : 20}

matplotlib.rc('font', **font)


ROR.fields_by_file = {"hydro" : [ROR.Scalar("d" ,  0), \
                                 ROR.Scalar("ux" , 1), \
                                 ROR.Scalar("uy" , 2), \
                                 ROR.Scalar("uz" , 3), \
                                 ROR.Scalar("P" ,  4), \
                                 ROR.Scalar("f" ,  5), \
                                 ROR.Scalar("tS",  6), \
                                 ROR.Scalar("tI",  7), \
                                 ROR.Scalar("tB",  8), \
                                 ROR.Scalar("g" ,  9), \
                                 ROR.Scalar("w" , 10) ],
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
          "g" : 'gamma', \
          "w" : 'wcr'}
color_palette = pylab.cm.cool
color_list = {"d": 'r', "ux": 'g', "uy": 'g', "uz": 'g', "P": 'b', "f": 'k', "tS": 'y', "tI": 'c', "tB": 'm', "g": 'k', "w": 'y'}

basedir = "/Users/dkosenko/Work/MultiFluidRAMSES/ramses/output/"

def connect(dir="backreaction_b5_3d3yr",out=9,var=["d","ux","uy","uz","P","f","tS","tI","tB","g","w"]):
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
    print "axes = ", pylab.gcf().axes
    if pylab.gcf().axes != []: 
        plot = pylab.gcf().axes[0]
    else: 
        plot = pylab.axes()
        #plot.set_xlabel(legend["r"])
        plot.set_xlabel("r, pc")
        if not func: plot.set_ylabel(legend[var])
    # plot
    rpc = data.points * output.info['unit_length'].express(C.pc)
    print rpc[-1], data_y[-1], len(rpc)
    print len(rpc[:len(rpc)-1])
    print len(data_y[:len(rpc)-1])
    if log: plot.semilogy(rpc[:len(rpc)-1],data_y[:len(rpc)-1],color=color)
    else:   plot.plot(    rpc,data_y,color=color)
    print max(rpc[-1])
    # pylab.xlim([0, 8])
    pylab.gcf().canvas.draw()
    matplotlib.pyplot.show()
    return data_y
    
def loop_profile(dir="backreaction_b5_3d3yr",out=[4,5,6,7,8,9],end=[1.,1.,1.],N=None,var='d',log=True,fig=1,over=False, eps="dens_profiles.eps"):
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
        pylab.hold(True)
    # colour scale
    if not over:
        bar_plot = figure.add_axes([0.92,0.1,0.015,0.8])
        norm = matplotlib.colors.Normalize(vmin=out[0], vmax=out[-1])
        color_bar = matplotlib.colorbar.ColorbarBase(bar_plot, cmap=color_palette, norm=norm, orientation='vertical')
        color_bar.set_label('output')
    matplotlib.pyplot.savefig(eps)
    matplotlib.pyplot.show()



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
    #print data
    return data

def show_slice(var="d",log=False,end=[1,1,1], rmax=15.0, ax=None, time='1000'):
    print var
    N = numpy.sqrt(data.npoints)
    map = data[var].reshape(N,N)
    map_ej = data["f"].reshape(N,N)

    #print data.npoints, 'points; N=', N
    points2d = 1. * numpy.array([end]*int(N))
    for i in range(int(N)): points2d[i] *= (i+1.)/int(N)
    # sample data
    data2d = pymses.analysis.sample_points(grid, points2d)
    rpc = data2d.points * output.info['unit_length'].express(C.pc) # rpc is 3D[len x 3] now
    
    if log: map = numpy.log(map)
    if rmax > rpc[-1,1]:
        dx = rmax-rpc[-1,1]

        newx = rpc[:,1]
        newmap = numpy.zeros(shape=(len(newx)+1,len(newx)+1))

        newx = numpy.append(newx, rmax)
        for i in range(len(newx)-1):
            m1 = map[i,:]
            m1 = numpy.append(m1, m1[-1])
            newmap[i,:] = m1
        print len(newmap[:,1]),'x',len(newmap[1,:]),'=',numpy.size(newmap)
        
        newmap[:,len(newx)-1] = numpy.ones(len(newx))*newmap[-1,len(newx)-1]
        #newmap = numpy.append(newmap, m2)
        
        img = pylab.pcolormesh(newx,newx,newmap)
        #img = pylab.pcolormesh(newx,newx,newmap, cmap='Spectral')
    else: 
        img = pylab.pcolormesh(rpc[:,1],rpc[:,2],map)
        #img = pylab.pcolormesh(rpc[:,1],rpc[:,2],map, cmap='Spectral')
    #img = pylab.imshow(map,origin='lower')
    pylab.xlim([0, rmax])
    pylab.ylim([0, rmax])
    matplotlib.pyplot.xlabel('r, pc')
    matplotlib.pyplot.ylabel('r, pc')
    jet = pylab.colorbar(cax=ax)
    #matplotlib.pyplot.show()
    return [map,img,jet]


def read_times(dir="test/"):
    f = open(basedir+dir+'params.dat')
    for line in f:  
        if re.search('tout',line): break

    line = line.strip('\n')
    tdata = line.split('=  ')
    tseq = tdata[1].split(', ')
    print tseq
    f.close()
    return tseq

def anim_slices(dir="bckreact_b5_n01_E1_M14_t4d3",out=[2,3,4,5,6,7,8,9],end=[1.,1.,1.],N=None,var='d',log=True,fig=1,over=False, rmax=7.0):
#def anim_slices(dir="backreaction_b5_3d3yr",out=[3,6],end=[1.,1.,1.],N=None,var='d',log=True,fig=1,over=False, rmax=15.0, mp4="dens_maps.mp4"):
    # reading time intervals
    tseq = read_times(dir=dir+'/')

    pylab.close(fig)
    figure = matplotlib.pyplot.figure(num=fig)
    #a = matplotlib.pyplot.gca()
    a = None
    
    ims = []
    for i in out:
        print 'output #',i, '  age:',tseq[i-1]
        connect(dir=dir,out=i,var=[var,'f'])
        load_slices()
        [map,img, jet] = show_slice(var=var,log=False,rmax=rmax,ax=a,time=tseq[i-1])
        #ann = pylab.annotate(tseq[i-1]+' yrs', xy=(.5, .6),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', color="white") 
        ann = pylab.annotate(tseq[i-1]+' yrs', xy=(.8, .9),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', color="white") 
	ims.append((img,ann,))
	a = jet.ax
	#ims.append((img,ann,jet,))
	# artificially expand the longer steps
        #if eval(tseq[i-1])>(5.e2-1.0): 
            #ims.append((img,ann))
        #    if eval(tseq[i-1])>1.5e3-1.0:  ims.append((img,ann))
	fname = '_tmp%03d.png'%i
	print 'Saving frame', fname
	figure.savefig(basedir+dir+'/'+fname)
	pylab.hold(False)

#    im_ani = anim.ArtistAnimation(figure, ims, interval=200, repeat_delay=500, blit=True)

#    print 'Making movie animation.mpg - this make take a while'
#    os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=AAC -oac copy -o animation.mpg")
#    os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=AAC -oac copy -o animation.mpg")
#    matplotlib.pyplot.show()


def eth_evol(dir="test", fig=1, hold=False):
    import evol as e

    rc('text', usetex=True)
    rc('font', family='serif')

    # reading time intervals
    tseq = read_times(dir=dir+'/')

    rslast = 0.0
    fslast = 0.0
    eth = []
    ecr = []
    egs = []
    ekn = []
    esum = []
    print tseq, len(tseq)
    for j in range(2,len(tseq)):
        eth.append(0.0)
        ecr.append(0.0)
        egs.append(0.0)
        ekn.append(0.0)
        esum.append(0.0)

        print 'output #',j, '  age:',tseq[j-1]
        connect(dir=dir,out=j)
        data = load_profiles()
        den = data['d']
        ux = data['ux']
        uy = data['uy']
        uz = data['uz']
        data_p = data['P']
        data_g = data['g']
        data_w = data['w']

        rpc = data.points * output.info['unit_length'].express(C.pc)
        rcm = data.points * output.info['unit_length'].express(C.cm)

        for k in range(len(rpc)-1):
            # if (data_w[k] > 1.e-10): print 'wcr=',data_w[k]
            rad = numpy.sqrt(numpy.sum(rcm[k]*rcm[k]))
            rad_next = numpy.sqrt(numpy.sum(rcm[k+1]*rcm[k+1]))
            ee = data_p[k]/(data_g[k] - 1.0)*rad*rad*(rad_next - rad)*4.0*numpy.pi*1.e-12
            eth[j-2] += ee
            egs[j-2] += ee*(1.0 - data_w[k])
            ecr[j-2] += ee*data_w[k]
            
            vel = numpy.sqrt(ux[k]*ux[k] + uy[k]*uy[k] + uz[k]*uz[k])*1.e5
            ekn[j-2] += den[k]*vel*vel*rad*rad*(rad_next - rad)*4.0*numpy.pi*0.5*1.67e-24
            esum[j-2] = eth[j-2] + ekn[j-2]
    if hold: 
        pylab.close(fig)
        figure = matplotlib.pyplot.figure(num=fig)
    pylab.loglog(tseq[2:],eth, 'r--', label='$\mathrm{E_{th}}$')
    pylab.hold(True)
    pylab.loglog(tseq[2:],egs, 'b:', label='$\mathrm{E_{gas}}$')
    pylab.loglog(tseq[2:],ecr, 'b-.', label='$\mathrm{E_{cr}}$')
    pylab.loglog(tseq[2:],ekn, 'g--', label='$\mathrm{E_{kin}}$')
    pylab.loglog(tseq[2:],esum, 'k-', label='$\mathrm{E_{tot}}$')
    #pylab.xlim([50.0, 2.0e4])
    pylab.xlim([100.0, 2.0e4])
    pylab.ylim([1.e49, 3.e51])
    #pylab.ylabel('E,  erg', fontsize=18)
    pylab.xlabel('t, years')
    pylab.ylabel('E,  erg')
    if hold: 
        pylab.hold(hold)
    else:
        pylab.legend()
    
        matplotlib.pyplot.savefig(basedir+dir+'/eng_t.pdf')
        matplotlib.pyplot.show()


def den_evol(dir="test", fig=1, hold=False):
    import evol as e

    rc('text', usetex=True)
    rc('font', family='serif')

    # reading time intervals
    tseq = read_times(dir=dir+'/')

    rslast = 0.0
    fslast = 0.0
    vard = []
    print tseq, len(tseq)
    for j in range(2,len(tseq)):
        print 'output #',j, '  age:',tseq[j-1]
        connect(dir=dir,out=j)
        data = load_profiles()
        den = data['d']
        ux = data['ux']
        uy = data['uy']
        uz = data['uz']

        rpc = data.points * output.info['unit_length'].express(C.pc)
        rcm = data.points * output.info['unit_length'].express(C.cm)

        dd = max(den) - min(den)
        vtot = 0.0
        mtot = 0.0
        for k in range(len(rpc)-1):
            rad = numpy.sqrt(numpy.sum(rcm[k]*rcm[k]))
            rad_next = numpy.sqrt(numpy.sum(rcm[k+1]*rcm[k+1]))
            vv = rad*rad*(rad_next - rad)*4.0*numpy.pi
            vtot += vv
            mtot += den[k]*vv
        # var(den)/mean(den)
        print dd, mtot, vtot
        vard.append(dd/(mtot/vtot))
    if hold: 
        pylab.close(fig)
        figure = matplotlib.pyplot.figure(num=fig)
    pylab.loglog(tseq[2:],vard, 'k-')
    pylab.xlim([100.0, 4.0e3])
    #pylab.ylim([1.e49, 1.3e51])
    #pylab.ylabel('E,  erg', fontsize=18)
    pylab.xlabel('t, years')
    pylab.ylabel('$\mathrm{\delta n/<n>}$')
    if hold: 
        pylab.hold(hold)
    else:
        matplotlib.pyplot.savefig(basedir+dir+'/varn_t.pdf')
        matplotlib.pyplot.show()

