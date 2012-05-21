import numpy
import matplotlib   
import matplotlib.animation as anim
import pylab
import pymses
    
import re   

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 20}

matplotlib.rc('font', **font)

#  fontsize=16
#  pylab.plot([1,2,3],[4,5,6])
#  ax = pylab.gca()
#  for tick in ax.xaxis.get_major_ticks():
#    tick.label1.set_fontsize(fontsize)
#  for tick in ax.yaxis.get_major_ticks():
#    tick.label1.set_fontsize(fontsize)


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


def plotevol(dir='test/', var='radius', fig=1):
    [t, rs, fs] = dataset(dir=dir, var=var)

    pylab.close(fig)
    figure = matplotlib.pyplot.figure(num=fig)
    pylab.semilogx(t,fs,label='FS')
    pylab.hold(True)
    if var != 'Mach': 
       pylab.semilogx(t,rs,label='RS')
    if var == 'velocity': 
       pylab.ylim([1.e1, 1.e4])
    pylab.xlim([200.0, 5.0e3])
    pylab.xlabel('t, years')
    ylab = var
    ylab=ylab.replace('_','_{')
    if '_{' in ylab: ylab=ylab+'}'
    pylab.ylabel('$\mathrm{'+ylab+'}$')
    pylab.legend(loc='lower left')
    var=var.replace('/','_')
    print var
    matplotlib.pyplot.savefig(basedir+dir+var+'_t.pdf')
    #matplotlib.pyplot.show()
    
