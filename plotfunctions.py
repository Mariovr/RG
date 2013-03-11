#!/usr/bin/env python
import numpy as np
import pylab as pl
import os,sys,shutil
import re

import tdasolver as tda

def plot_2col(name,tit = 'title',xname = 'independendvar',yname = 'dependendvar'):
  print('start with the generation of a plot with plot_2col in rgfunctions.py')
  plotdata = open(name,'r')
  aa = []
  bb = []
  for line in plotdata:
    if (line[0] is '#'):
      continue
    data = line.split()
    try:
      bb.append(float(data[1]))
      aa.append(float(data[0])) 
    except:
      pass
  pl.figure()
  pl.plot(aa,bb,'r')
  pl.xlabel(xname)
  pl.ylabel(yname)
  pl.title(tit)
  pl.savefig(name+'.png' )
  pl.close()
  
  
def generate_plot(nlevel,npair,dvar,afhxas,name = 'plotenergy.dat',plotg = True):
  """
  some nice plots to visualize the data with matplotlib, plotg = true if you plot the energylevels of the reduced
  BCS Hamiltonian of the sp levels of a geometry file
  """
  print ('start with the generation of plots')
  plotdata = open(name, 'r')
  name = name.rstrip('.dat')
  aa = []
  bb = []
  cc = []
  dd = []
  ee = []
  for line in plotdata:
    if (line[0] is '#'):
      continue
    data = line.split()
    try:
      bb.append(float(data[1]))
      aa.append(float(data[0])) 
      cc.append(float(data[2]))
      dd.append(float(data[3]))
      ee.append(float(data[4]))
    except:
      pass
    ##cc.append(float(data[2])) 
  pl.figure()
  pl.plot(aa,bb,'r')
  if plotg: 
    pl.xlabel(afhxas)
    pl.ylabel('condensation energy (a.u.)')
    pl.title('interaction constant = %f (a.u.)' %(dvar))
    pl.savefig('c%s.png' %name )
    pl.close()
  else:
    pl.xlabel('g (a.u.)')
    pl.ylabel('condensation energy (a.u.)')
    pl.title('%s = %f (a.u.)' %(afhxas,dvar))
    pl.savefig('c%s.png' %name )
    pl.close()

  pl.figure()
  pl.plot(aa,cc,'r')
  if plotg:
    pl.xlabel(afhxas)
    pl.ylabel('energy spectrum (a.u.)')
    pl.title('interaction constant = %f (a.u.)' %(dvar))
    pl.savefig('g%s.png' %name)
    pl.close()
  else:
    pl.xlabel('g (a.u.)')
    pl.ylabel('groundstate energy (a.u.)')
    pl.title('%s = %f (a.u.)' %(afhxas,dvar))
    pl.savefig('g%s.png' %name ) 
    pl.close()
  if plotg:
    pl.figure()
    pl.plot(aa,dd,'r')
    pl.xlabel(afhxas)
    pl.ylabel('energy of the non-interacting groundstate (a.u.)')
    pl.title('aantal paren = %f' %(npair))
    pl.savefig('nig%s.png' %name)
    pl.close()
    pl.figure()
    try:
      pl.plot(aa,ee)
      pl.xlabel(afhxas)
      pl.ylabel("d (a.u.)")
      pl.title("number of sp levels = %f" %nlevel)
      pl.savefig("%sd.png" % 'meanleveldistance')
      pl.close()
    except:
      print 'the plot of d failed'

def plotrgvars(apair,ref = 'plotenergy.dat',afhvar = 'g (a.u.)',namerg = 'rgvar',stop = None,begin = 0,istart = 3):
  plotf = open(ref, 'r')
  dataf = np.loadtxt(plotf,comments = '#')
  pl.figure()
  for i in xrange(istart,2*apair+istart,2):
    if stop is None:
      pl.plot(dataf[begin:,0],dataf[begin:,i])
    else:
      pl.plot(dataf[begin:stop,0],dataf[begin:stop,i])
  pl.xlabel(afhvar)
  pl.ylabel('real part rgvars (a.u.)')
  pl.title('Richardson-Gaudin variables')
  pl.savefig('re%s.png' %namerg )
  pl.close()
  pl.figure()
  for i in xrange(istart+1,2*apair+istart+1,2):
    if stop is None:
      pl.plot(dataf[begin:,0],dataf[begin:,i])
    else:
      pl.plot(dataf[begin:stop,0],dataf[begin:stop,i])
  pl.xlabel(afhvar)
  pl.ylabel('imaginary part rgvars (a.u.)')
  pl.title('Richardson-Gaudin variables')
  pl.savefig('im%s.png' %namerg)
  pl.close()

def plotrgvarsxi(apair,eendlev,g,ref = 'rergvarxi.dat' , namerg = 'rgvarXI' , redbcs = True,eta = 1.):
  plotre = open(ref, 'r')
  datare = np.loadtxt(plotre,comments = '#')
  pl.figure()
  irange = np.arange(3,2*apair+3,2)
  def calc_sinr(i,red):
    if red:
      return eendlev[i]*2.
    else:
      return eendlev[i]*eendlev[i]*eta
  for i in irange:
    pl.plot(datare[0,i],datare[0,i+1],'g.', markersize = 10)
    pl.plot(datare[len(datare[:,0])-1,i],datare[len(datare[:,0])-1,i+1],'r.',mfc = 'None', markersize = 10)
    pl.plot(datare[:,i],datare[:,i+1],'b-')
  for i in range(len(eendlev)):
    pl.axvline(x = calc_sinr(i,redbcs) , c= 'k' ,linestyle = '--')
  pl.xlabel('real part of rgvars (a.u)')
  pl.ylabel('imaginary part of rgvars (a.u.)')
  pl.title('Richardson-Gaudin variables at g = %f (xi in [0,1])' %(g))
  #pl.xlim((2*eendlev[0]-5*(eendlev[1]-eendlev[0]),2*eendlev[-1]+0.5))
  #pl.ylim((-20,20))
  pl.savefig('%s.png' %namerg )
  pl.close() 
  
def plotrgcloud(name,npair,g,sen):
  namef = open(name,'r')
  namerg = 'rgcould'+name
  pl.figure()
  start = False
  for line in namef:
    if line[0] is '#':
      seniori = re.compile(r'seniority\s\[(.+?)\]')
      d = seniori.search(line)
      if d:
	d = d.group(1)
	lis = re.findall(r'\d',d)
	s = list(lis)
        print 'we found the following seniority: %s' %s
        senlist = map(int,s)      
        if senlist == sen:
	  start = True
        elif start == True and senlist is not sen:
	  break
      continue
    elif start == True:
      index = line.index('}')
      line = line[index+2:]
      larray = line.split()      
      revar = map(float,larray[2:2+npair])
      imvar = map(float,larray[2+npair+1:])
      for i in range(npair):
        pl.plot(revar[i],imvar[i],'g.')       
  pl.xlabel('real part of rgvars (a.u)')
  pl.ylabel('imaginary part of rgvars (a.u.)')
  pl.title('RG variables at g = %f for all ex. states' %(g))
  ax = pl.gca()
  pl.text(0.65,0.85,'sen ='+ str(sen),transform = ax.transAxes)
  pl.savefig('%s%s.png' %(namerg,sen ))
  pl.close()

def generatePlotExited(nlevel,npair,afh,afhxas):
  """
  some nice plots to visualize the data with matplotlib
  """
  gfile = open('plotenergygrondtoestand.dat','r')
  gb = []
  gc = []
  gd = []
  for line in gfile:
    if (line[0] is '#'):
      continue
    data = line.split()
    gb.append(float(data[1]))
    gc.append(float(data[2]))
    gd.append(float(data[3]))
  
  plotfiles =  [x if '.dat' in x and 'plotenergy' in x and 'grondtoestand' not in x and 'rgvar' not in x else 0 for x in os.listdir(os.getcwd())]
  plotfiles = filter (lambda a: a != 0, plotfiles)
  print plotfiles
  pl.figure(1)
  pl.clf()
  pl.figure(2)
  pl.clf()
  pl.figure(3)
  pl.clf()
  aa = []
  bb = []
  cc = []
  dd = []
  ee = []
  for i in range(len(plotfiles)):
    bb.append([])
    cc.append([])
    dd.append([])
    ee.append([])
    aa.append([])
  i = 0
  for name in plotfiles:
    plotdata = open(name, 'r')
    j = 0
    for line in plotdata:
      if (line[0] is '#'):
        continue
      data = line.split()
      try: #super ugly needs to be improved (error prone)
        bb[i].append(float(data[1])-gb[j])
        aa[i].append(float(data[0])) 
        cc[i].append(float(data[2])-gc[j])
        dd[i].append(float(data[3])-gd[j])
      except:
	pass 
      j += 1
    i += 1
  pl.close()
  for i in range(len(plotfiles)):
    fig1= pl.figure(1)
    pl.plot(aa[i],bb[i],'r')  
    fig2= pl.figure(2)
    pl.plot(aa[i],cc[i],'r')
    fig3= pl.figure(3)
    pl.plot(aa[i],dd[i],'r')
    
  pl.figure(1)
  pl.xlabel('g(a.u.)')
  pl.ylabel('the condensation energy spectrum (a.u.)')
  pl.title('dependend value = %f (a.u.)' %(afh))
  pl.savefig('condensation energy spectrum.png' )
  pl.figure(2)
  pl.xlabel('g(a.u.)')
  pl.ylabel('the energyspectrum (a.u.)')
  pl.title('dependend value = %f (a.u.)' %(afh))
  pl.savefig('spectrum.png')
  pl.figure(3)
  pl.xlabel('g (a.u.)')
  pl.ylabel('the non interacting ground state (a.u.)')
  pl.title('number of pairs = %f' %(npair))
  pl.savefig('non-interactinggroundstate.png') 
  
def addlevel(interactionconstant = -0.137):
  
  os.chdir(os.path.join(os.getcwd(), 'dangspectrum5dif'))
  d = [i for i in os.listdir(os.getcwd()) if 'plotenergy' in i and '.dat' in i]
  
  def getal(string):
    regexp = re.compile('\d+')
    match = regexp.search(string)
    g = match.group()
    return int(g)
  
  d = sorted(d,key = getal, reverse = True)
  print d
  grondtoestand = []
  for j in d:
    filep = open(j,'r')
    for line in filep:
      if line[0] == '#':
	continue
      data = line.split()
      if float(data[0]) == interactionconstant:
        grondtoestand.append(float(data[2]))
    filep.close()
      
  print grondtoestand 
  x = range(0,len(grondtoestand))
  print x
  
  A = np.array([ x, np.ones(len(x))])
  # linearly generated sequence
  a,f,g,h = np.linalg.lstsq(A.T,grondtoestand) # obtaining the parameters
  print 'de gevonden rechte = %.10f x + %.10f' %(a[0] , a[1])
  lined = map(lambda g: a[0]*g +a[1],x) # regression line
  print lined
 
  pl.figure()
  pl.plot(x,lined,'r-',label = '%f*x %f' %(a[0],a[1]))
  pl.plot(x,grondtoestand, 'bo',label= 'datapoints')
  pl.legend()
  pl.xlabel('number of added continuum sp levels')
  pl.ylabel('groundstate energy (MeV)')
  pl.title('the groundstate energy of Sn120 with pairingsconstant %f' %interactionconstant)
  pl.savefig('g=%faddlevel.png' %interactionconstant)
  pl.show()
  pl.close()  
  
def main(option, args):
  if option == 'pexcited':
    nlevel = sys.argv[0]
    npair = sys.argv[1]
    afh = sys.argv[2]
    afhxas = sys.argv[3]
    generatePlotExited(nlevel,npair,afh,afhxas)
    
  if option == 'wpairing':
    nlevel,npair,dvar,afhxas = args[0:4]
    try:
      name,nig,plotg = args[4:]
    except:
      pass
    generate_plot(nlevel,npair,dvar,afhxas,name = 'plotenergy.dat',plotg = False)
  
  if option == 'addlevel':
    g = args[0]
    addlevel(interactionconstant = g)
    
  if option == 'rgvar':
    apair, exname = args[0:2]
    ref = 'rergvar.dat'; imf = 'imrgvar.dat' ; afhvar = 'g' ; namerg = 'rgvar' 
    try:
      begin = args[2]
      stop = args[3]
    except:
      begin = 0
      stop = None
    if exname is not None:
      ref = 'rergvar%s.dat' %exname
      imf =  'imrgvar%s.dat' %exname
      namerg = 'rgvar%s' %exname
    plotrgvars(apair,ref = ref,imf = imf,afhvar = 'g' , namerg = namerg,begin = begin,stop = stop)
  
  if option is 'rgcloud':
    name, npair,g,sen = args 
    plotrgcloud(name,npair,g,sen)
  return 0
  
def importmain():
  option = 'rgcloud'
  args = 'DangSn120neutronwindow(5,5)sen=2.dat' , 9, -0.137,[0,0,0,2,0,0]
  main(option,args)
  
if __name__ == '__main__':
  '''
  possible options: 'pexcited' plots all the excited states relative to the groundstate, 'wpairing' plots the 
  results from a writepairing call in writepairing.py(main), 'addlevel' from a set of outputfiles from generating_datak generated
  by adding empty sp levels and get from those files the groundstate energy at a constant g and plot them and perform lin.regression 
  '''
  '''
  option = 'addlevel'
  #args = (6,10,'x','g')
  args = -0.137
  main(option,args)
  '''
  option = 'rgcloud'
  args = 'DangSn120neutronwindow(5,5)sen=2.dat' , 9, -0.137, [1,1,0,0,0,0]
  main(option,args)
  
  
  