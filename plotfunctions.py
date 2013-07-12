#!/usr/bin/env python
import numpy as np
import pylab as pl
import os,sys,shutil
import re

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

def plot_spectrumxichange(dirname,search):
  """
  Plot the entire spectrum at a particular g in function of xi
  args = directory with all the data
  """
  dirlist =	os.listdir(dirname)
  filelist = [i for i in dirlist if search in i]
  pl.figure(0)
  xidata = []
  edata = []
  countgood = 0
  countbad = 0
  for j in filelist:
    file = open(j,'r')
    print 'working in file %s' %(j)
    for line in file:
      try: float(line[0])
      except: continue
      data = line.split()
      xidata.append(float(data[0]))
      edata.append(float(data[1]))
    if xidata[-1] == 1.: 
      pl.plot(xidata,edata,'b') 
      countgood += 1
      print j , countgood , 'good solution'
    else: 
      pl.plot(xidata,edata,'r') 
      print j, countbad, 'bad solution'
      countbad += 1
    file.close()
    xidata = [] ; edata = []
  
  print 'We found %g good solutions and %g tda startdistributions that broke down before xi = 1, we hope that\'s what you expected' %(countgood,countbad)
  pl.xlabel(r'$\xi$')
  pl.ylabel(r'energy spectrum (a.u.)')
  pl.title(r'All tda start distributions $\xi$')
  #Create custom artists
  goodline = pl.Line2D((0,1),(0,0), color='b') 
  badline = pl.Line2D((0,1),(0,0), color='r')
  
  pl.legend([goodline,badline],['solution','breakdown'])
  pl.savefig('spectrum%s.png' %(search.translate(None,'.')))
  pl.close()
  
def generate_plot(nlevel,npair,dvar,name = 'plotenergy.dat',plotg = False):
  """
  some nice plots to visualize the data with matplotlib, plotg = true if you plot the energylevels of the sp levels of a geometry file
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
    #cc.append(float(data[2])) 
  #plot of groundstate energy
  pl.figure()
  pl.plot(aa,bb,'r')
  pl.xlabel(dvar+' (a.u.)')
  pl.ylabel('condensation energy (a.u.)')
  pl.savefig('c%s.png' %name )
  pl.close()
  #plot of condensation energy
  pl.figure()
  pl.plot(aa,cc,'r')
  pl.xlabel(dvar+' (a.u.)')
  pl.ylabel('energy (a.u.)')
  pl.savefig('g%s.png' %name)
  pl.close()
  if plotg:
    pl.figure()
    pl.plot(aa,dd,'r')
    pl.xlabel(dvar)
    pl.ylabel('energy of the non-interacting groundstate (a.u.)')
    pl.title('aantal paren = %f' %(npair))
    pl.savefig('nig%s.png' %name)
    pl.close()
    pl.figure()
    try:
      pl.plot(aa,ee)
      pl.xlabel(dvar)
      pl.ylabel("d (a.u.)")
      pl.title("number of sp levels = %f" %nlevel)
      pl.savefig("%sd.png" % 'meanleveldistance')
      pl.close()
    except:
      print 'the plot of d failed'

def plotrgvarscplane(apair, interval = (-20 , 0) ,infile = 'plotenergy.dat', name = 'rgplane', istart = 3 ):
  dataf , kpunten = readdata(infile, fast = False)
  pl.figure()
  for j in xrange(len(kpunten)):
    for i in xrange(istart,2*apair+istart,2):
      pl.plot(dataf[kpunten[j][1] + interval[0]:kpunten[j][1] + interval[1],i],dataf[kpunten[j][1]+interval[0]:kpunten[j][1] + interval[1],i+1] , 'b')
    pl.xlabel('real part rgvars (a.u.)')
    pl.ylabel('imaginary part rgvgrs (a.u.) ')
    pl.title('Richardson-Gaudin variables')
    pl.savefig('complexplane%d.png' % kpunten[j][0])
    pl.close()

def readdata(ref = 'plotenergy.dat' , fast =True):
  plotf = open(ref, 'r')
  kpunten = []
  try:
    if fast == False:
      raise ValueError
    dataf = np.loadtxt(plotf,comments = '#')
  except:
    dataf = np.array([])
    plotf.seek(0) #go back to beginning of file
    for line in plotf:
      if line[0] == '#':
        analyse = re.search(r'^#\s+((-|\d)\d+\.*\d*)\s+kritisch',line)
        if analyse:
          kpunten.append(float(analyse.group(1)))
        continue
      pline = np.array(map(float,line.split()))
      if len(dataf) <= 1:
        dataf = pline
      else:
        try:
          dataf = np.vstack((dataf,pline))
        except:
          continue
  plotf.close()
  return dataf, kpunten

def plotintofmotion(alevel, apair , ref = 'plotenergy.dat' , afhvar =  'g (a.u.)',namerg = 'integralsofmotion',stop = None,begin = 0,istart = 3):
  istart = istart + 2*apair
  dataf , kpunten = readdata(ref)
  pl.figure()
  for i in xrange(istart,alevel+istart):
    if stop is None:
      pl.plot(dataf[begin:,0],dataf[begin:,i],'b')
    else:
      pl.plot(dataf[begin:stop,0],dataf[begin:stop,i])
  pl.xlabel(afhvar)
  pl.ylabel('integrals of motion (a.u.)')
  pl.title('integrals of motion of the Richardson-Gaudin model')
  pl.savefig('%s.png' %namerg )
  pl.close()

def plotrgvars(apair,ref = 'plotenergy.dat',afhvar = 'g (a.u.)',namerg = 'rgvar',stop = None,begin = 0,istart = 3 , cplane = False):
  dataf , kpunten = readdata(ref)
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
  if cplane:
    pl.figure()
    for i in xrange(istart,2*apair+istart,2):
      if stop is None:
        pl.plot(dataf[begin:,i],dataf[begin:,i+1],'b')
      else:
        pl.plot(dataf[begin:stop,0],dataf[begin:stop,i])
    pl.xlabel('real part rgvars (a.u.)')
    pl.ylabel('imaginary part rgvgrs (a.u.) ')
    pl.title('Richardson-Gaudin variables')
    pl.savefig('complexplane%s.png' %namerg)
    pl.close()

def plotrgvarsxi(apair,eendlev,g,ref = 'rergvarxi.dat' , namerg = 'rgvarXI'):
  plotre = open(ref, 'r')
  datare = np.loadtxt(plotre)
  pl.figure()
  irange = np.arange(1,2*apair+1,2)
  for i in irange:
    pl.plot(datare[0,i],datare[0,i+1],'g.', markersize = 10)
    pl.plot(datare[len(datare[:,0])-1,i],datare[len(datare[:,0])-1,i+1],'r.',mfc = 'None', markersize = 10)
    pl.plot(datare[:,i],datare[:,i+1],'b-')
  for i in range(len(eendlev)):
    pl.axvline(x = eendlev[i]*2. ,c=  'k',linestyle = '--')
  pl.xlabel('real part of rgvars (a.u)')
  pl.ylabel('imaginary part of rgvars (a.u.)')
  pl.title('Richardson-Gaudin variables at g = %f (xi in [0,1])' %(g))
  pl.xlim((2*eendlev[0]-5*(eendlev[1]-eendlev[0]),2*eendlev[-1]+0.5))
  pl.ylim((-20,20))
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

def generatePlotExited(nlevel,npair):
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
  pl.savefig('condensation energy spectrum.png' )
  pl.figure(2)
  pl.xlabel('g(a.u.)')
  pl.ylabel('the energyspectrum (a.u.)')
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
    nlevel,npair,dvar = args[0:3]
    try:
      name,nig,plotg = args[3:]
    except:
      pass
    print nlevel,npair,dvar
    generate_plot(nlevel,npair,dvar,name = 'plotenergy.dat',plotg = False)

  if option == 'addlevel':
    g = args[0]
    addlevel(interactionconstant = g)
    
  if option == 'rgvar':
    apair, exname = args[0:2]
    ref = 'plotenergy.dat';  afhvar = 'g' ; namerg = 'rgvarnauw' 
    try:
      begin = sys.argv[1]
      stop = args[3]
    except:
      begin = 0
      stop = None
    plotrgvars(apair,ref = ref,afhvar = 'g' , namerg = namerg,begin = begin,stop = None,istart=3)
  
  if option is 'rgcloud':
    name, npair,g,sen = args 
    plotrgcloud(name,npair,g,sen)
  
  if option is 'cprgvar':
    npair , name = args
    plotrgvarscplane(npair , (-20,0), infile = name)
  
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
  option = 'cprgvar'
  #args =allstates0.dat, 6,4000, None
  args = 7, 'plotenergy7.dat'
  main(option,args)
  #plot_spectrumxichange(sys.argv[1],sys.argv[2])  
