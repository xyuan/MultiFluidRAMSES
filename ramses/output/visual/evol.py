import numpy
import matplotlib   
import matplotlib.animation as anim
import pylab
import pymses
    
import re   

basedir = "/Users/dkosenko/Work/MultiFluidRAMSES/ramses/output/"
dir = "test/"
    
#### time evolutions

def evol(dir='test/', var='radius', fig=1):

    f = open(basedir+dir+'shocks.dat')
    header = f.readline()
    header = header.strip('\n')
    header = header.replace('->','')
    header = header.replace('shock:','')
    keys = header.split('  ')
    key = []
    rs = {}
    fs = {}
    for k in keys:
	if len(k)>1:
	   k=k.strip()
	   fs[k] = []
	   rs[k] = []
	   key.append(k)

    data = f.read()
    lines = data.split('\n')
    print len(lines)
    print key
    f.close()
	
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

    pylab.close(fig)
    figure = matplotlib.pyplot.figure(num=fig)
    pylab.semilogx(t,fs[var],label='FS')
    pylab.hold(True)
    if var != 'Mach': 
       pylab.semilogx(t,rs[var],label='RS')
    pylab.xlabel('t, years')
    pylab.ylabel(var)
    pylab.legend()
    var=var.replace('/','_')
    print var
    matplotlib.pyplot.savefig(basedir+dir+var+'_t.pdf')
    #matplotlib.pyplot.savefig(basedir+dir+var+'_t.eps')
    matplotlib.pyplot.show()
    
