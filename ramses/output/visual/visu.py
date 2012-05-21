import numpy
import matplotlib
import matplotlib.animation as anim
import pylab
import pymses
import os, sys

import re

from matplotlib import rc
from scipy import fftpack

from pymses.utils import constants as C
from pymses.filters import CellsToPoints
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

def connect(dir="backreaction_b5_3d3yr",out=9,var=["d","ux","uy","uz","P","f","tS","tI","tB","g","w"], verbose=True):
    global output, grid
    output = pymses.RamsesOutput(basedir+dir, out)
    if verbose: print output.info
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
    print 'color=', col
    if log: plot.semilogy(rpc[:len(rpc)-1],data_y[:len(rpc)-1],color=color)
    else:   plot.plot(    rpc,data_y,color=color)
    print max(rpc[-1])
    # pylab.xlim([0, 8])
    pylab.gcf().canvas.draw()
    #matplotlib.pyplot.show()
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



def slice(var='d',z=0.0,log=True):
    #cam = Camera(center=[0.5, 0.5, 0.5], line_of_sight_axis='z', region_size=[1., 1.], log_sensitive=log)
    cam = Camera(center=[0.5, 0.5, z], line_of_sight_axis='z', region_size=[1., 1.], log_sensitive=log)
    op = ScalarOperator(lambda dset: dset[var])
    map = SliceMap(grid, cam, op, z=z)
    return map

def load_slices(z=0.0, size=10):
    global data
    map_size = 2**size
    camera = Camera(center=[0.5, 0.5, z], line_of_sight_axis='z', region_size=[1., 1.], map_max_size=map_size)
    points = camera.get_slice_points(0.0)
    data = sample_points(grid, points)
    #print data
    return data

def show_slice(var="d",log=False,end=[1,1,1], rmax=15.0, vmin = 1.e-3, vmax=2.e0, alpha=0.1, ax=None, ejecta=False, show=False):
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
        newx = rpc[:,1]
        newx = numpy.append(newx, rmax)
        newmap = numpy.zeros(shape=(len(newx),len(newx)))
        newmap_ej = numpy.zeros(shape=(len(newx),len(newx)))
        for i in range(len(newx)-1):
            m1 = map[i,:]
            m1 = numpy.append(m1, m1[-1]) # add an extra cell - double of the last one
            newmap[i,:] = m1
            me1 = map_ej[i,:]
            me1 = numpy.append(me1, me1[-1]) # add an extra cell - double of the last one
            newmap_ej[i,:] = me1
        print len(newmap[:,1]),'x',len(newmap[1,:]),'=',numpy.size(newmap)
        
        newmap[:,len(newx)-1] = numpy.ones(len(newx))*newmap[-1,len(newx)-1]
        newmap_ej[:,len(newx)-1] = numpy.ones(len(newx))*newmap_ej[-1,len(newx)-1]
        
        pylab.hold(True)
        abg = alpha*1.5
        pylab.pcolormesh(newx,newx,newmap, vmin=vmin, vmax=vmax, alpha=abg)
        pylab.hold(True)
    if (ejecta):
        pylab.pcolormesh(rpc[:,1],rpc[:,2],map_ej, cmap='YlOrBr', vmin=1.0e-3, vmax=1.0)
    img = pylab.pcolormesh(rpc[:,1],rpc[:,2],map,  vmin=vmin, vmax=vmax, alpha=alpha)

    pylab.hold(False)
    pylab.xlim([0, rmax])
    pylab.ylim([0, rmax])
    matplotlib.pyplot.xlabel('r, pc')
    matplotlib.pyplot.ylabel('r, pc')
    jet = pylab.colorbar()
    if (show):
        matplotlib.pyplot.show()
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

def anim_slices(dir="bckreact_b5_n01_E1_M14_t4d3",out=[2,3,4,5,6,7,8,9],end=[1.,1.,1.],N=None,var='d',log=True,fig=1,over=False, rmax=7.0, vmax=2.0, alpha=0.6, ejecta=True, dpi=300):
#def anim_slices(dir="backreaction_b5_3d3yr",out=[3,6],end=[1.,1.,1.],N=None,var='d',log=True,fig=1,over=False, rmax=15.0, mp4="dens_maps.mp4"):
    # reading time intervals
    tseq = read_times(dir=dir+'/')

    #a = matplotlib.pyplot.gca()
    a = None
    
    ims = []
    for i in out:
        print 'output #',i, '  age:',tseq[i-1]
        connect(dir=dir,out=i,var=[var,'f'])
        load_slices()
        pylab.close(fig)
        figure = matplotlib.pyplot.figure(num=fig)
        pylab.hold(False)
        [map, img, jet] = show_slice(var=var,log=False,rmax=rmax,ax=a,ejecta=ejecta, vmax=vmax, alpha=alpha)
        #ann = pylab.annotate(tseq[i-1]+' yrs', xy=(.5, .6),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', color="white") 
        ann = pylab.annotate(tseq[i-1]+' yrs', xy=(.8, .9),  xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', color="white") 
	fname = '_tmp%03d.png'%i
	print 'Saving frame', fname
	figure.savefig(basedir+dir+'/'+fname, dpi=dpi, transparent=True)
	pylab.hold(False)

#    im_ani = anim.ArtistAnimation(figure, ims, interval=200, repeat_delay=500, blit=True)

#    print 'Making movie animation.mpg - this make take a while'
#    os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=AAC -oac copy -o animation.mpg")
#    os.system("mencoder 'mf://_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=AAC -oac copy -o animation.mpg")
#    matplotlib.pyplot.show()


def load_3D(dir='test', var='d', rmax=6.0, fig=1, N = 2):
    global data

    pylab.close(fig)
    figure = matplotlib.pyplot.figure(num=fig)
    a=None
    
    for i in range(2*N+1):
        z = 0.5*i/N
        print 'z=', z
	load_slices(z)
        [map, img, jet] = show_slice(var=var,log=False,rmax=rmax, ax=a)
	a = jet.ax
        fname = '_ztmp%3f.png'%z
        print 'Saving frame', fname
        figure.savefig(basedir+dir+'/'+fname)

def sum_3D(dir='test', out=1, size=6, foe=1.e51, slice=False):
    if (slice):
        Nz = 1
    else:
        Nz = 2**(size-1)
    z = []
    ethz  = 0.0
    ecrz  = 0.0
    egsz  = 0.0
    eknz  = 0.0
    volz  = 0.0
    masz  = 0.0
    
    dmax = []
    dmin = []
    
    den3d = []
    rad = []
    # slices along z-axis
    for k in range(2*Nz+1):
        z.append(0.5*k/Nz)
    for k in range(2*Nz):
        dmaxz = []
        dminz = []
	end = [1.,1.,z[k]]
        print 'z=', z[k]
	connect(dir=dir,out=out,verbose=False)
	data = load_slices(z[k], size=size)
	N = numpy.sqrt(data.npoints)
	print 'Nz=', Nz, '  Nxy=', N, ':   k=',k
	xy = 1. * numpy.array([end]*int(N))
	for i in range(int(N)): xy[i][0:2] *= (i+1.)/int(N)
	rad = xy*output.info['unit_length'].express(C.pc)

        den = data['d'].reshape(N,N)
        ux = data['ux'].reshape(N,N)
        uy = data['uy'].reshape(N,N)
        uz = data['uz'].reshape(N,N)
        data_p = data['P'].reshape(N,N)
        data_g = data['g'].reshape(N,N)
        data_w = data['w'].reshape(N,N)
        den3d.append(den)
        

        ethx = 0.0
        ecrx = 0.0
        egsx = 0.0
        eknx = 0.0
        volx = 0.0
        masx = 0.0
        for i in range(len(xy)-1): # x
            dmaxz.append(max(den[i,:]))
            dminz.append(min(den[i,:]))

            ethy = 0.0
            ecry = 0.0
            egsy = 0.0
            ekny = 0.0
            voly = 0.0
            masy = 0.0
            for j in range(len(xy)-1): # y 
		dy = (xy[j+1][1] - xy[j][1])
		pdy = data_p[i][j]*dy  # integrate over y-axis
		ethy += pdy/(data_g[i][j] - 1.0)
		#ecry += ee*data_w[i][j] # this is just not true and totally wrong
		egsy += pdy*(1.0 - data_w[i][j])/(5.0/3.0 - 1.0)
		if (data_w[i][j]>0.9):
		    print '***************************************************'
		    print '********* w_CR: ',k,i,j,data_w[i][j], '***********'
		    print '***************************************************'

		vel2 = ux[i][j]**2 + uy[i][j]**2 + uz[i][j]**2
		ekny += den[i][j]*vel2*dy
		
		voly += dy
		masy += den[i][j]*dy

	    ethx += ethy * (xy[i+1][0] - xy[i][0])  # integrate over x-axis
	    #ecrx += ecry * (xy[i+1][0] - xy[i][0])  # integrate over x-axis
	    egsx += egsy * (xy[i+1][0] - xy[i][0])  # integrate over x-axis
	    eknx += ekny * (xy[i+1][0] - xy[i][0])  # integrate over x-axis
	    volx += voly * (xy[i+1][0] - xy[i][0])  # integrate over x-axis
	    masx += masy * (xy[i+1][0] - xy[i][0])  # integrate over x-axis
	    #print 'x: ',i, 'eknx=', eknx, ',  x=', xy[i][0]
        dmax.append(max(dmaxz))
        dmin.append(min(dminz))

	ethz += ethx * (z[k+1] - z[k]) # integrate over z-axis
	#ecrz += ecrx * (z[k+1] - z[k]) # integrate over z-axis
	egsz += egsx * (z[k+1] - z[k]) # integrate over z-axis
	eknz += eknx * (z[k+1] - z[k]) # integrate over z-axis
	volz += volx * (z[k+1] - z[k]) # integrate over z-axis
	masz += masx * (z[k+1] - z[k]) # integrate over z-axis
	#print 'z: ',k, 'eknz=', eknz, ',  z=', z[k], ',  dz=',(z[k+1] - z[k])

    ethz *= 1.e-12*8.0*output.info['unit_length'].express(C.cm)**3/foe
    #ecrz *= 1.e-12*8.0*output.info['unit_length'].express(C.cm)**3/foe
    egsz *= 1.e-12*8.0*output.info['unit_length'].express(C.cm)**3/foe
    eknz *= 1.e10*0.5*1.67e-24*8.0*output.info['unit_length'].express(C.cm)**3/foe
    esum = ethz + eknz
    
    print 'TOTAL ENERGY: ',esum, foe, 'erg'
    
    ddmax = max(dmax)
    ddmin = min(dmin)
    
    return [ethz, egsz, eknz, (ddmax-ddmin)/(masz/volz), den3d, rad]

def fftden(dir="test", out = 10, fig=1, size=6, col=None):
    rc('text', usetex=True)
    rc('font', family='serif')
    foe = 1.e51

    [a,b,c,r,d3d,rad] = sum_3D(dir=dir,out=out,size=size,foe=foe)
    fd3d = fftpack.fftn(d3d)
    nfd3d = numpy.abs(fd3d)
    wnum = []
    for i in range(len(rad)): wnum.append(1.0/rad[i][0])

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
        plot.set_xlabel("k, pc$^{-1}$")
        plot.set_xlim([0.03, 50.0])
    # plot
    plot.loglog(wnum, nfd3d[:][0][0], color=color)
    pylab.gcf().canvas.draw()

    return [wnum, nfd3d, plot]

def loop_fft(dir="test",out=[4,5,6,7,8,9], size=6,fig=1, over=False, eps="fft_profiles.pdf"):
    # set figure
    if not over:
        pylab.close(fig)
        figure = matplotlib.pyplot.figure(num=fig)

    # loop over outputs
    for i in out:
        print 'output #',i
        if len(out)>1: col = -(i-out[0])/(1.*out[-1]-out[0])
        else: col = 0
        #print len(out), col, i, out[0]
        [wnum, nfd3d, plot] = fftden(dir=dir, out=i, size=size, fig=fig, col=col)
        pylab.hold(True)
    
    # colour scale
    if not over:
        bar_plot = figure.add_axes([0.92,0.1,0.015,0.8])
        norm = matplotlib.colors.Normalize(vmin=out[0], vmax=out[-1])
        color_bar = matplotlib.colorbar.ColorbarBase(bar_plot, cmap=color_palette, norm=norm, orientation='vertical')
        color_bar.set_label('output')

    # initial power spectrum
    dinit = 0.1 * numpy.ones(shape=nfd3d.shape)
    finit = fftpack.fftn(dinit)
    nfinit = numpy.abs(finit)
    #print nfinit[:][0][0]
    #print wnum
    pylab.hold(True)
    plot.loglog(wnum, nfinit[:][0][0], color='k')
    kolm = []
    for w in wnum: kolm.append(1.0e4*w**(11.0/3.0))
    plot.loglog(wnum, kolm, color='k')
    print kolm

    #pylab.xlim([0.3, 30])
    #pylab.xlabel('k, pc$^{-1}$')
    pylab.ylabel('|P|')

    matplotlib.pyplot.savefig(basedir+dir+'/'+eps)
    matplotlib.pyplot.show()


def eth_evol(dir="test", size=6, fig=1, fhold=False, col='k'):
    rc('text', usetex=True)
    rc('font', family='serif')
    foe = 1.e51
    lw = 2
        
    # reading time intervals
    tseq = read_times(dir=dir+'/')
    print tseq, len(tseq)
    eth = []
    egs = []
    ecr = []
    ect = []
    ekn = []
    esum = []
    den = []

    fev = open(basedir+dir+'/eng_t.dat','w')
    fev.write(' '+str(2**size)+' '+str(2**size)+' '+str(2**(size-1))+"\n")
    for m in range(2,len(tseq)):
        print '########'
        print 'output #', m, '  age:',tseq[m-1]
        [a,b,c,r,d3d,rad] = sum_3D(dir=dir,out=m,size=size,foe=foe)
        eth.append(a)
        egs.append(b)
        ecr.append(a-b)
        ect.append(1.0-b/a)
        ekn.append(c)
        esum.append(a + c)
        den.append(r)
        print a,b,c, r
	fev.write(" "+tseq[m-1]+"\t"+str('%9.5e' % a)+"\t"+str('%9.5e' % b)+"\t"+str('%9.5e' % c)+"\t"+str('%9.5e' % r)+"\n")
	
    fev.close()
    figure = matplotlib.pyplot.figure(num=fig)
    if fhold: 
        pylab.loglog(tseq[2:],eth, col+'-.', linewidth=lw)
        pylab.hold(True)
        pylab.loglog(tseq[2:],egs, col+':', linewidth=lw)
        pylab.loglog(tseq[2:],ect, col+':', linewidth=1)
        pylab.loglog(tseq[2:],ekn, col+'--', linewidth=lw)
        pylab.loglog(tseq[2:],esum, col+'-', linewidth=lw)
        #pylab.loglog(tseq[2:],den, col+'-', label='$\mathrm{\delta n/n}$')
        pylab.xlim([100.0, 2.0e4])
        pylab.ylim([0.03, 2.0])
        #pylab.ylabel('E,  erg', fontsize=18)
        pylab.xlabel('t, years')
        pylab.ylabel('$\mathrm{E}, 10^{51}\mathrm{erg}$')

        pylab.hold(fhold)
    else:
        pylab.loglog(tseq[2:],eth, col+'-.', linewidth=lw, label='$\mathrm{E_{th}}$')
        pylab.hold(True)
        pylab.loglog(tseq[2:],egs, col+':', linewidth=lw, label='$\mathrm{E_{gas}}$')
        pylab.loglog(tseq[2:],ect, col+':', linewidth=1, label='$\mathrm{E_{CR}}$')
        pylab.loglog(tseq[2:],ekn, col+'--', linewidth=lw, label='$\mathrm{E_{kin}}$')
        pylab.loglog(tseq[2:],esum, col+'-', linewidth=lw, label='$\mathrm{E_{tot}}$')
        #pylab.loglog(tseq[2:],den, col+'-', label='$\mathrm{\delta n/n}$')
        pylab.xlim([100.0, 2.0e4])
        pylab.ylim([0.03, 2.0])
        #pylab.ylabel('E,  erg', fontsize=18)
        pylab.xlabel('t, years')
        pylab.ylabel('$\mathrm{E}, 10^{51}\mathrm{erg}$')
        pylab.hold(False)

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

