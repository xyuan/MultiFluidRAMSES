"""
plotting SNR profiles from Chevalier
Gilles Ferrand, 2012/03/01
"""

import numpy
import pylab

skip = {'old':5, 'new':2}
order = {'old': ["i","r","u","d","P"], 'new': ["r","d","u","P"]}
style = {'old': "*", 'new': "-"}
label = {'old': "integration", 'new': "interpolation",\
          'd': "density", 'u': "velocity", 'P': "pressure"}
colour = {'d': "g", 'u': "b", 'P': "r"}


def plot(dir=".",fname="hydro_cr_V2",types=['old','new'],vars=['d','u','P'],x=[None,None],y=[None,None],xlog=False,ylog=True,norm=True,legend=True,fig=1,size=(20,10),ext='png'):
  """ plots the profiles for each quantity in 'vars' as a function of radius, for different 'types' of data files """
  # DATA
  data = {}
  for type in types:
    filename = dir + "/" + fname + type + ".dat"
    print "loading ",filename
    data[type] = numpy.loadtxt(filename,skiprows=skip[type])
  # FIGURE
  figure = pylab.figure(num=fig,figsize=size)
  if xlog: 
    if ylog: pylab_plot = pylab.loglog
    else:    pylab_plot = pylab.semilogx
  else:
    if ylog: pylab_plot = pylab.semilogy
    else:    pylab_plot = pylab.plot
  plot = {}
  for type in types:
    plot[type] = {}
    i_ref = numpy.where(data[type][:,order[type].index('u')]>0)[0][-1]
    for var in vars:
      j_var = order[type].index(var)
      if norm: 
        print var,"[",i_ref,"] = ", data[type][i_ref,j_var]
        ref = data[type][i_ref,j_var]
      else: 
        ref = 1.
      plot[type][var] = pylab_plot(data[type][:,order[type].index('r')],data[type][:,j_var]/ref,colour[var]+style[type],label=var)
  pylab.axis([x[0],x[1],y[0],y[1]])
  pylab.xlabel("radius [pc]")
  pylab.ylabel("hydrodynamic quantity [a.u.]")
  if legend: 
    legend1 = pylab.legend([plot[types[t ]][vars[-1]] for t in range(len(types))],[label[types[t]] for t in range(len(types))],loc=4)
    legend2 = pylab.legend([plot[types[-1]][vars[ v]] for v in range(len(vars ))],[label[vars [v]] for v in range(len(vars ))],loc=1)
    pylab.gca().add_artist(legend1)
  # FILE
  filename = dir + "/" + fname + "."+ext
  print "Writing ",filename
  pylab.savefig(fname)
  if fig<0: pylab.close(fig)
  #return data
  