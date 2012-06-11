import numpy
import matplotlib   
import matplotlib.animation as anim
import pylab
import pymses
    
import re   

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

#  fontsize=16
#  pylab.plot([1,2,3],[4,5,6])
#  ax = pylab.gca()
#  for tick in ax.xaxis.get_major_ticks():
#    tick.label1.set_fontsize(fontsize)
#  for tick in ax.yaxis.get_major_ticks():
#    tick.label1.set_fontsize(fontsize)

global lw
lw=2

basedir = "/Users/dkosenko/Work/MultiFluidRAMSES/ramses/output/"
dir = "test/"
    
#### time evolutions
def readdata(dir='test/'):
    f = open(basedir+dir+'shocks.dat')
    header = f.readline()
    header = header.strip('\n')
    header = header.replace('->','')
    header = header.replace('shock:','')
    keys = header.split('  ')
    key = []
    for k in keys:
	if len(k)>1:
	   k=k.strip()
	   key.append(k)

    data = f.read()
    lines = data.split('\n')
    print len(lines)
    print key
    f.close()
    return [key, lines]
    

def dataset(dir='test/', var='radius'):
    [key, lines] = readdata(dir=dir)

    rs = {}
    fs = {}
    for k in key:
	fs[k] = []
	rs[k] = []

    t = []
    j = -1
    for i in range(len(lines)-1):
        #print 'i: ', i, ' lines:', lines[i]
        line = lines[i].split('  ')
        #print 't=',t[j],'; ',line
	if i == (j+1)*4:
	   t.append(eval(lines[i]))
	   j = j+1    
	if i == j*4+1:
	   for k in key:
	       val = line.pop(0)
	       while val == '' or re.search('\-1',val): val = line.pop(0)
	       if re.search('\*',val): val = '0.0'
               #print k,': ',val
	       rs[k].append(eval(val))
	if i == j*4+3:
	   for k in key:
	       val = line.pop(0)
	       while val == '' or re.search('\+1',val): val = line.pop(0)
	       if re.search('\*',val): val = '0.0'
               #print k,': ',val
	       fs[k].append(eval(val))

    print len(t),', ',len(rs[var])
    return [t, rs[var], fs[var]]


def plotevol(dir='test', var='radius', ymax=500.0, col='k', lab='test', fig=1, fhold=False, lloc='lower left'):
    [t, rs, fs] = dataset(dir=dir+'/', var=var)

    #if not fhold: pylab.close(fig)
    #figure = matplotlib.pyplot.figure(num=fig)
    pylab.semilogx(t,fs,col,linewidth=lw,label=lab)
    pylab.hold(fhold)
    #if var != 'Mach': 
    #   pylab.semilogx(t,rs,label='RS')
    #if var == 'velocity': 
    #   pylab.ylim([1.e1, 1.e4])
    pylab.xlim([10.0, 5.0e3])
    pylab.xlabel('t, years')
    ylab = var
    ylab=ylab.replace('_','_{')
    if '_{' in ylab: ylab=ylab+'}'
    pylab.ylim([0.0, ymax])
    pylab.ylabel('$\mathrm{'+ylab+'}$')
    pylab.legend(loc=lloc)
    var=var.replace('/','_')
    print var
    if not fhold: matplotlib.pyplot.savefig(basedir+dir+'/'+var+'_t.eps')
    #if not fhold: matplotlib.pyplot.show()
    
def plotevolall(var='radius', ymax=500.0, lloc='lower left'):
    c = ['g', 'r', 'm', 'b']
    t = ['test6nx003', 'test6nfr02', 'test6nfr1', 'test6nfr5', 'test6nfr0']
    b = ['breact6nx003', 'breact6nfr02', 'breact6nfr1', 'breact6nfr5', 'breact6nfr0']
    
    #fig=1
    #for i in range(len(c)):
    #    plotevol(dir=t[i], col=c[i], fig=fig, fhold=True)
    #plotevol(dir=t[-1], fig=fig)

    fig=2
    figure = matplotlib.pyplot.figure(num=fig)
    for i in range(len(c)):
        plotevol(dir=b[i], var=var, ymax=ymax, col=c[i], lab=b[i], lloc=lloc, fig=fig, fhold=True)
    plotevol(dir=b[-1], var=var, ymax=ymax, lab=b[-1], lloc=lloc, fig=fig)
                