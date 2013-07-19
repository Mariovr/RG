#!/usr/bin/env python
import numpy as np
import pylab as pl
import os,sys,shutil
import re

def make_newstyle(file , comment = '#', rgindex = 3):
  """
  Function that changes the order of the Richardson-Gaudin variables instead of the old style reval reval reval ... Im imval imval imval ...
  to the new style reval imval reval imval .....
  """
  filehandler = open(file , 'r')
  savenew = open('newstyle' + file , 'w')
  for line in filehandler:
    if line[0] == comment:
      savenew.write(line)
      continue
    data = line.split('\t')
    del(data[rgindex+1]) #delete 'Im'
    rerg = data[rgindex].split()
    imrg =data[rgindex+1].split()
    newrgstr = data[0:rgindex] + [re + '\t' + im for re, im in zip(rerg , imrg)]
    savenew.write('\t'.join(newrgstr) + '\n')
  filehandler.close()
  savenew.close()

    
class Plot_RG_Files(object):
  def __init__(self):
    """
    initialisation of the plotting class of the Richardson-Gaudin solver, this boils
    down to the creation of a figure to which we add one axes
    """
    self.fig = pl.figure(1,figsize=(8,6), dpi=80 , frameon = True , facecolor = '0.75' , edgecolor = 'w')
    self.fig.add_subplot(111 , axisbg = 'w' , projection = 'rectilinear') #if you want to add axes on particular place: fig.add_axes([0.15, 0.1, 0.7, 0.3]) where ->  [begin , bottom to start axes , width , height ]
    self.separated = True #if we have a list and need to plot the plots separated

  def readrgeq(self, filen):
    """
    Some plotting functions need a bit more information this info is extracted from the header of the files
    """
    elev = None 
    with open(filen, 'r') as file:
      for line in file:
        if line.startswith('#'):
          if 'energylevels' in line:
            match = re.search(r':\s*\[(?P<elev>[+\-\d\s.]*)\]\s*' , line) 
            if match:
              elev = list(match.group('elev').split())
          if 'pairs' in line:
            match = re.search(r'pairs[:=\s]*(\d+)',line)
            if match:
              self.apair = int(match.group(1))
          if 'change' in line:
            match = re.search(r'change is\s*:*\s*(\w+)',line)
            if match:
              self.dvar = match.group(1) + ' (a.u.)'
          if 'eta' in line:
            match = re.search(r':*\s*(\d+)', line) #after some time remove first star because then the new __str__ function of the RichardsonEq class is already established so the match is better
            self.eta = int(match.group(1))
        else:
          break
    if elev != None:
      self.alevel = len(elev)
      self.energiel = map(float,elev) 

    print  elev 

  def readdata(self, reflist , comment = '#' , regexp = None , substr = None):
    """
    read data from files in reflist is really important this function has to be called every time because it gathers the plotting data.
    If you put a regular expression the function reads the data in line by line and extracts information it finds in between the data: if the regular expression is definedRemark ,the line by line reading is also used when the data is not compatible with np.loadtxt
    possible regexp are : r'^#\s+((-|\d)\d+\.*\d*)\s+kritisch' to find the critical points 
r'seniority\s\[(.+?)\]' to find the seniority's in an allstates file 
    """
    self.kpunten = []
    datalist = []
    prefixlist = []
    if os.path.isfile(str(reflist)):
      reflist = [reflist] #if we work with only one file this wraps it automatically in right format
    for ref in reflist:
      print('start with the collection of data from file %s' %ref)
      plotf = open(ref, 'r')
      prefixlist.append(re.sub('\.dat$' , '' , ref))
      try:
        if regexp != None:
          raise ValueError
        dataf = np.loadtxt(plotf,comments = comment)
      except:
        dataf = np.array([])
        kpuntenf = []
        plotf.seek(0) #go back to beginning of file
        for line in plotf:
          if regexp is not None:
            analyse = re.search(regexp,line)
            if analyse:
              kpuntenf.append((analyse.group(1), len(dataf) ))
              print 'we found the following matches: %s' % analyse.group(0)
            continue
          if substr != None: 
            line = re.sub(substr, '' , line)
          pline = np.array(map(float,line.split()))
          if len(dataf) <= 1:
            dataf = pline
          else:
            try:
              dataf = np.vstack((dataf,pline))
            except:
              continue

        self.kpunten.append(kpuntenf)
      datalist.append(dataf)

    plotf.close()
    self.readrgeq(reflist[0])
    self.datarg = datalist
    self.prefix = prefixlist

  def readfiles(self, dirname , search , filelist = None , notsearch = 'rgvar' , sortfunction = None , rev = False , regexp = None , substr = None):
    """
    If you want to plot data from a single file use readdata instead, this is a wrapper for readdata if you want to plot data from multiple files
    """
    print('Remark we change self.rgdata because we are reading in new files')
    dirlist =  os.listdir(dirname)
    plotfiles =  [x  for x in dirlist if '.dat' in x and search in x and notsearch not in x]
    if filelist != None:
      plotfiles += filelist

    if sortfunction != None:
      plotfiles = sorted(plotfiles ,key = sortfunction , reverse = rev )

    print plotfiles
    self.readdata(plotfiles , regexp = regexp ,  substr = substr)

  def generate_plot(self):
    """
    some nice plots to visualize the data with matplotlib, plotg = true if you plot the energylevels of the sp levels of a geometry file
    """
    print ('start with the generation of plots')
    #plot of condensation energy
    self.plotwrap(0,1, 'condensation energy (a.u.)' , 'ce',titel = 'the condensation energy (a.u.)')
    self.plotwrap(0,2, 'energy (a.u.)' , 'ge', titel = 'the energy (a.u.)' )

  def plotwrap(self , xindex , yindex , yas , name = None , titel = None ,color = 'r' , sort = '' , label = None ):
    for i in range(len(self.datarg)):
      self.layout(self.dvar , yas , tit = titel)
      self.fig.axes[0].plot(self.datarg[i][:,xindex],self.datarg[i][:,yindex], color+sort , label = label)
      if self.separated == True:
        self.savefig(name, filenum = i)
    if self.separated == False:
      self.savefig(name + 'together')

  def plotrgvarscplane(self, interval = (-20 , 0), label = None):
    for k in xrange(len(self.datarg) ):
      for j in xrange(len(self.kpunten[filenum])):
        self.layout('real part rgvars (a.u.)' , 'imaginary part rgvgrs (a.u.) ', tit = 'Richardson-Gaudin variables')
        for i in xrange(self.rgindex,2*self.apair+self.rgindex,2):
          self.fig.axes[0].plot(self.datarg[k][self.kpunten[filenum][j][1] + interval[0]:self.kpunten[k][j][1] + interval[1],i],self.datarg[k][self.kpunten[k][j][1]+interval[0]:self.kpunten[k][j][1] + interval[1],i+1] , 'b' , label = label)
        self.savefig('%f' % (float(self.kpunten[k][j][0])), filenum = k) # you never want this together

  def plotrgvars(self,cplane = False , begin = 0 , stop = None):
    print('starting to plot the Richardson-Gaudin variables')
    self.plotrgwrap(self.rgindex, 2*self.apair+self.rgindex , self.dvar , 'real part rgvars (a.u.)',tit =  'Richardson-Gaudin variables', name = 're' , begin = begin , stop =stop)
    self.plotrgwrap(self.rgindex+1, 2*self.apair+self.rgindex+1 , self.dvar ,'imaginary part rgvars (a.u.)', tit = 'Richardson-Gaudin variables', name = 'im', begin = begin , stop = stop )
    if cplane:
      self.plotrgwrap(self.rgindex, 2*self.apair+self.rgindex ,'real part rgvars (a.u.)' ,'imaginary part rgvars (a.u.)',tit =  'Richardson-Gaudin variables', name = 'cp', begin= begin , stop = stop)

  def plotrgwrap(self, columnstart ,columnend , xas , yas , tit = None , begin = 0  , stop = None, name = '' , color = 'b' , sort = '.' ,label = None):
    for j in xrange(len(self.datarg)):
      for i in xrange(columnstart,columnend,2):
        if 'cp' in name:
          self.fig.axes[0].plot(self.datarg[j][begin:stop,i],self.datarg[j][begin:stop,i+1], color+sort , label = label)
        else:
          self.fig.axes[0].plot(self.datarg[j][begin:stop,0],self.datarg[j][begin:stop,i] , label = label)
      self.layout(xas , yas , tit = tit)
      if self.separated == True:
        self.savefig(name , filenum = j)
    if self.separated == False:
      self.savefig(name + 'together')
    
  def plotintofmotion(self,namerg = 'iom',stop =None,begin = 0 , xlim = None , ylim = None):
    columns = self.rgindex + 2 * self.apair
    for j in xrange(len(self.datarg)):
      for i in xrange(columns,self.alevel+ columns):
        self.fig.axes[0].plot(self.datarg[j][begin:stop,0],self.datarg[j][begin:stop,i],'b')
      self.layout(self.dvar , 'integrals of motion (a.u.)', tit = 'integrals of motion of the Richardson-Gaudin model' , xlim = xlim , ylim = ylim)
      if self.separated == True:
        self.savefig(namerg , filenum = j)
    if self.separated == False:
      self.savefig(namerg + 'together')


  def normalize_to_groundstate(self):
    print('Warning we normalize all the excited states to the groundstate energy')
    gronddat = self.datarg[0]
    for i in range(1,len(self.datarg)):
      dif = np.shape(gronddat )[0] - np.shape(self.datarg[i])[0]
      print dif
      if dif < 0 :
        self.datarg[i] = self.datarg[i][0:dif ,:] 
      elif dif > 0:
        gronddat = gronddat[: -1.*dif , :]
      print np.shape(gronddat) , np.shape(self.datarg[i])
      self.datarg[i][:,1:] = self.datarg[i][:,1:] - gronddat[:,1:] #we made sure that the data of the groundstateenergy is first in the rgdata list
    del(self.datarg[0], self.prefix[0])

  def layout(self ,  xlab , ylab , xlim = None , ylim = None , tit = None , axnum = 0 , legendhand = None , legendlab = None , legendpos = 'best' , finetuning = False):
    print('We are starting with the layout')
    self.fig.axes[axnum].set_xlabel(xlab)
    self.fig.axes[axnum].set_ylabel(ylab)
    if xlim != None:
      self.fig.axes[axnum].set_xlim(xlim) #good value for xlim in the case of a xi path is : (2*self.rgeq.energiel[0]-5*(self.rgeq.energiel[1]-self.rgeq.energiel[0]),2*self.rgeq.energiel[-1]+0.5) 
    if ylim != None:
      self.fig.axes[axnum].set_ylim(ylim)
    if tit != None:
      self.fig.axes[axnum].set_title(tit)
    if legendlab != None:
      self.fig.axes[axnum].legend(legendhand , legendlab, loc = legendpos)
    """
    if you forgot to add a label to a line with linenumber: lnum you can do: self.fig.axes[axnum].lines[lnum].set_label('this is my new label')
    the underneath is the same as : h , l = self.fig.axes[axnum].get_legend_handles_labels()
                                    self.fig.axes[axnum].legend(h,l)
    """
    leg = self.fig.axes[axnum].legend(loc = legendpos) #draws the legend on axes[axnum] all the plots that you labeled are now depicted in legend
    if finetuning == True:
      # the matplotlib.patches.Rectangle instance surrounding the legend
      frame  = leg.get_frame()  
      frame.set_facecolor('0.80')    # set the frame face color to light gray

      # matplotlib.text.Text instances you can change all properties of labels
      for t in leg.get_texts():
        t.set_fontsize('small')    # the legend text fontsize

      # matplotlib.lines.Line2D instances
      for l in leg.get_lines():
        l.set_linewidth(1.5)  # the legend line width

  def savefig(self , name , filenum = 0):
    """
    After we are satisfied with our figure we save it with this function: dpi = pixels per inch, we delete all . in the prefix because it will give faults when we use
    the savefig function of fig
    """
    self.fig.savefig('%s%s.png' %(self.prefix[filenum].translate(None,'.'), name ), dpi = 80 , facecolor = 'w' , edgecolor = 'w')
    self.fig.clf()
    self.fig.add_subplot(111 , axisbg = 'w' , projection = 'rectilinear')

  def least_sqr_fit(self,x, y):
    """
    Calculates the least square fit of a list of independend variables x and dependend variables y.
    It returns a list of function values of the best fitted straight line, with the given x values as independend variables and also a list with the parameters
    that define the line. It's also possible to fit at the same time multiple datasets with the same xvalues just give y the form [(v1 , v2 , v3) , (v1 , v2 , v3), ... ]
    Where the first tuple consists of the function values of x1 the second of x2 .... , So you get immediately three fitted lines, with the coefficients in a[0][0] , a[0][1]
    , a[0][2] for the first, second and third rico for the three lines same for the bisection point with y axis
    """
    A = np.array([ x, np.ones(len(x))])
    # linearly generated sequence
    a,f,g,h = np.linalg.lstsq(A.T,y) # obtaining the parameters
    print 'de gevonden rechte = %.10f x + %.10f' %(a[0], a[1])
    lined = map(lambda g: a[0]*g +a[1],x) # regression line
    return lined , a

  def standard_plot(self , rgw = True , intm = True):
    self.generate_plot()
    if rgw:
      self.plotrgvars(cplane = False , begin = 0 , stop = None)
    if intm:
      self.plotintofmotion()


class Plot_Geo_File(Plot_RG_Files):
  """
  remark before the Richardson-Gaudin variables start this file has 6 columns, extra: nig , meandistance , number of levels
  """
  def __init__(self , name = 'x', dvar = None):
    self.dvar = dvar 
    self.rgindex = 6
    if os.path.isdir(name):
      self.readfiles(pname , 'plotenergy')
    elif os.path.isfile(name):
      self.readdata([pname])
    super(Plot_Geo_File,self).__init__()

  def generate_plot(self):
    super(Plot_Geo_File,self).generate_plot(dvar)
    print('plot non-interacting groundstate')
    self.plotwrap(0,3, 'energy of the non-interacting groundstate (a.u.)','nig', titel = 'aantal paren = %f' %(self.apair))
    try:
      self.plotwrap(0,4,"d (a.u.)" ,'meandistance', titel = "number of sp levels = %f" %self.alevel)
    except:
      print 'the plot of d failed'

class Plot_Data_File(Plot_RG_Files):
  def __init__(self, name = 'x'):
    self.rgindex = 3
    if os.path.isdir(name):
      self.readfiles(pname , 'plotenergy')
    elif os.path.isfile(name):
      self.readdata([pname])
    super(Plot_Data_File,self).__init__()

  def addlevel(self , g ):
    genergy = [k[0][0] for k in self.kpunten]
    x = range(0,len(genergy ))
    y , coef = self.least_sqr_fit(x,genergy ) 
    self.fig.axes[0].plot(x, y ,'r-',label = '%f*x %f' %(coef[0],coef[1]))
    self.fig.axes[0].plot(x, genergy, 'bo',label= 'datapoints')
    print genergy
    self.layout('number of added continuum sp levels', 'groundstate energy (MeV)', tit =  'the groundstate energy of Sn120 with i.c.:  %.3f' %g )
    self.savefig('g=%fal.png' % g)

class Plot_Xi_File(Plot_RG_Files):
  def __init__(self, chardata):
    self.rgindex = 2
    self.chardata = chardata
    self.dvar = 'xi (a.u.)'
    super(Plot_Xi_File,self).__init__()

  def plot_spectrumxichange(self):
    """
    Plot the entire spectrum at a particular g in function of xi
    args = directory with all the data
    """
    countgood = 0 ; countbad = 0
    for idata in self.datarg:
      if idata[-1, 0] == 1.: 
        self.fig.axes[0].plot(idata[0:,0], idata[0: ,1] ,'b') 
        countgood += 1
        print  countgood , 'good solution'
      else: 
        self.fig.axes[0].plot(idata[0:,0], idata[0: ,1] ,'r') 
        print countbad, 'bad solution'
        countbad += 1
    print 'We found %g good solutions and %g tda startdistributions that broke down before xi = 1, we hope that\'s what you expected' %(countgood,countbad)
    #Create custom artistsr[goodline,badline],['solution','breakdown']
    goodline = pl.Line2D((0,1),(0,0), color='b') 
    badline = pl.Line2D((0,1),(0,0), color='r')
    self.layout(self.dvar , r'energy spectrum (a.u.)' , tit = r'All tda start distributions $\xi$' , legendhand = [goodline , badline] , legendlab = ['solution', 'breakdown'] )
    self.savefig('xispec')

  def plotrgvarsxi(self, namerg = 'rgvxi' ):
    for j in xrange(len(self.datarg)):
      for i in np.arange(self.rgindex,2*self.apair+self.rgindex,2):
        self.fig.axes[0].plot(self.datarg[j][0,i],self.datarg[j][0,i+1],'g.', markersize = 10) #Richardson-Gaudin solutions (xi = 1)
        self.fig.axes[0].plot(self.datarg[j][len(self.datarg[j][:,0])-1,i],self.datarg[j][len(self.datarg[j][:,0])-1,i+1],'r.',mfc = 'None', markersize = 10) # Corresponding tda solutions (xi = 0 )
        self.fig.axes[0].plot(self.datarg[filenum][:,i],self.datarg[filenum][:,i+1],'b-') # intermediate values of xi
      if self.eta == None:
        sing = np.array(self.energiel)* 2
      else:
        sing = self.eta * np.array(self.energiel) * np.array(self.energiel)
      for i in range(self.alevel):
        self.fig.axes[0].axvline(x = sing[i] ,c=  'k',linestyle = '--')
      self.layout('real part of rgvars (a.u)', 'imaginary part of rgvars (a.u.)', tit = 'Richardson-Gaudin variables at g = %f (xi in [0,1])' %(self.chardata))
      if self.separated == True:
        self.savefig(namerg , filenum = j)
    if self.separated == False:
      self.savefig(namerg + 'together')

class Plot_All_File(Plot_RG_Files):
  def __init__(self, g):
    self.chardata = g
    self.rgindex = 2
    super(Plot_All_File,self).__init__()

  def plotrgcloud(self):
    """
    This function needs it own datareader because it's to specific
    """
    print self.kpunten
    for i in range(len(self.kpunten[0])):
      self.fig.axes[0].text(0.65,0.85,'sen ='+ self.kpunten[0][i][0],transform = self.fig.axes[0].transAxes)
      if i == len(self.kpunten[0]) -1 :
        end = None
      else:
        end = self.kpunten[0][i+1][1] + 1
      print end
      self.plotrgwrap( self.rgindex,2*self.apair+self.rgindex,'real part of rgvars (a.u)' , 'imaginary part of rgvars (a.u.)', tit ='RG vars g = %f all states'%(self.chardata) , begin = self.kpunten[0][i][1]  , stop = end , name = 'cpcloud'+ self.kpunten[0][i][0] , filenum = 0)
  
  
def main(option, args):
  plotter = Plot_Data_File()
  plottergeo = Plot_Geo_File()
  if option == 'pexcited':
    plotter = Plot_Data_File()
    plotter.readfiles(os.getcwd() ,'plotenergy', notsearch = 'rgvar' , sortfunction = lambda x : -1. if 'ground' in x else 0) #sortfunction makes sure the groundstate is first this is important for the normalization
    plotter.standard_plot(True, True)
    plotter.normalize_to_groundstate()
    plotter.separated = False
    plotter.generate_plot()
    
  if option == 'wpairing':
    if args[1] == True:
      plottergeo.readdata([args[0]])
      plottergeo.generate_plot()
    else:
      plotter.readdata([args[0]])
      plotter.generate_plot()
      
  if option == 'addlevel':
    plotter.readfiles( '.' , 'plotenergy' ,   sortfunction = lambda s : int(re.search(r'\d+' , s).group()), rev = True , regexp = r'^%f\s+[\-+\.\d]+\s+([\-+\.\d]+)\s' % args[0])
    plotter.addlevel(args[0])
    
  if option == 'rgvar':
    ref = args[0] 
    begin = 0
    stop = None
    cp = args[1]
    plotter.readdata([ref])
    plotter.plotrgvars(cplane = cp , begin = begin , stop = stop)
  
  if option is 'rgcloud':
    name = 'newstyleDang120neutronwin5_5sen2.dat'
    plottera = Plot_All_File(-0.137)
    plottera.readdata([name] , regexp = r'seniority\s\[(.+?)\]' , substr = r'\{.*\}')
    plottera.plotrgcloud()
  
  if option is 'cprgvar':
    ref = args[0]
    plotter.readdata([ref], regexp = r'^#\s+((-|\d)\d+\.*\d*)\s+kritisch', begin = 1)
    plotter.plotrgvarscplane(interval = (-20,0))

  if option is 'intmotion':
    plotter.readdata([args])
    plotter.plotintofmotion(xlim = (-0.5 , 0) , ylim = (-10 , 5))

  if 'xi' in option:
    plotterxi = Plot_Xi_File(args[1])
    if option is 'xipath':
      plotterxi.readfiles('.',args[0])
      plotterxi.plotrgvarsxi()

    if option is 'specxichange':
      plotterxi.readfiles('.', args[0])    
      plotterxi.plot_spectrumxichange()
  
def importmain():
  option = 'addlevel'
  args = -0.137
  main(option,args)
  
if __name__ == '__main__':
  '''
  possible options: 'pexcited' plots all the excited states relative to the groundstate, 'wpairing' plots the 
  results from a writepairing call in writepairing.py(main), 'addlevel' from a set of outputfiles from generating_datak generated
  by adding empty sp levels and get from those files the groundstate energy at a constant g and plot them and perform lin.regression 
  '''
  option = 'pexcited'
  #args =  -0.137 , None
  args = 'plotenergy.dat'
  #make_newstyle( 'Dang120neutronwin5_5sen2.dat', comment = '#', rgindex = 3)
  main(option,args)
