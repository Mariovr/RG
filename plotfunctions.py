# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Mario Van Raemdonck, 2013;
# (c) Ghent University, 2013
#!/usr/bin/env python
import numpy as np
import pylab as pl
import os,sys,shutil
import re
import matplotlib
import math

import datareader as dr

#matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']}) #adjust fonts
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#matplotlib.rc('text', usetex=True)

def makemovie(name = None):
  # makes a movie from all the .png files in the current directory
  print 'Starting to create a movie, with all the .png files in directory: %s ' %str(os.getcwd())
  if name != None:
    dirname = name
  else:
    dirname = str(os.getcwd())
  command = ('mencoder',
           'mf://*.png',
           '-mf',
           'type=png:w=800:h=600:fps=5',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           dirname+'.avi')

  os.spawnvp(os.P_WAIT, 'mencoder', command)

class File_Collector(object):
  def __init__(self, rootdir , search , notsearch = '.png' , notdir = 'xyvwa' , filelist = None , sortfunction = None , rev = False):
    if filelist != None:
      self.plotfiles = filelist
    else:
      self.plotfiles = []
    self.sortfunction = sortfunction
    self.readfiles(rootdir , search , notsearch = notsearch , notdir = notdir)
    self.sortplotfiles(rev)
    print self.plotfiles

  def addfiles(self , *args):
    for i in args:
      self.plotfiles.append(i)
  
  def sortplotfiles(self, rev = False):
    if self.sortfunction != None:
      self.plotfiles = sorted(self.plotfiles ,key = self.sortfunction , reverse = rev )
    else:
      print 'No sort function given so the order of the files doesn\'t matter for the figure'

  def readfiles(self, dirname , search , notsearch = 'rgvar' , notdir = 'xyvwa'):
    """
    If you want to plot data from a single file use readdata instead, this is a wrapper for readdata if you want to plot data from multiple files
    """
    print('We are in the following directory: %s looking for files that contain %s and not %s' %(dirname, search , notsearch))
    dirlist =  os.listdir(dirname)
    for filep in dirlist:
      filep = os.path.join(dirname,filep) 
      if os.path.islink(filep):
        pass
      elif os.path.isdir(filep):
        m = re.search(notdir , filep)
        if m is None:
          self.readfiles(filep , search, notsearch = notsearch, notdir = notdir )
      elif os.path.isfile(filep) and '.dat' in filep: 
        nm = re.search(notsearch, filep)
        m = re.search(search , filep)
        #print m , nm
        if m is not None and nm is None:
          self.plotfiles.append(filep)
      else:
        pass
    
class Plot_RG_Files(object):
  def __init__(self):
    """
    initialisation of the plotting class of the Richardson-Gaudin solver, this boils
    down to the creation of a figure to which we add one axes
    """
    self.fig = pl.figure(1,figsize=(8,6), dpi=80 , frameon = True , facecolor = '0.75' , edgecolor = 'w')
    self.fig.add_subplot(111 , axisbg = 'w' , projection = 'rectilinear') #if you want to add axes on particular place: fig.add_axes([0.15, 0.1, 0.7, 0.3]) where ->  [begin , bottom to start axes , width , height ]
    self.separated = True #if we have a list and need to plot the plots separated

  def add_axes(self,pos = [0.5 , 0.2 , 0.4 , 0.3], axisbg = None , projection = 'rectilinear'):
    self.fig.add_axes(pos , axisbg = axisbg, projection = projection)

  def readdata(self, reflist , comment = '#' , regexp = None , substr = None, filename = True):
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
      if not filename:
        prefixlist.append( os.path.dirname(ref) + '/')
      else:
        prefixlist.append(re.sub('\.dat$' , '' , ref))
      try:
        if regexp != None:
          raise ValueError
        dataf = np.loadtxt(plotf,comments = comment)
        print 'we readed data in with np.loadtxt'
      except:
        print('reading in data with numpy loadtxt failed or use reg exp to extract information')
        dataf = np.array([])
        kpuntenf = []
        plotf.seek(0) #go back to beginning of file
        for line in plotf:
          if regexp is not None:
            analyse = re.search(regexp,line)
            if analyse:
              kpuntenf.append((analyse.group(1),len(dataf) ))
              print 'we found the following matches: %s' % analyse.group(0)
          if substr != None: 
            line = re.sub(substr, '' , line)
          if line[0] != comment:
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
    self.datarg = datalist
    self.prefix = prefixlist
    self.reader = dr.ReaderOutput(reflist[0]) #Some plotting functions need a bit more information this info is extracted from the header of the files
    #self.reader.depvar['depvar'] += ' (a.u.)'

  def procesfiles(self, dirname , search , notsearch = r'\.sw*|\.png', notdir = 'awfwfr', sortfunction = None , rev = False , regexp = None , substr = None , filelist = None , filename = True):
    filecol =File_Collector(dirname , search , notsearch = notsearch ,filelist = filelist , sortfunction = sortfunction , rev =rev )
    self.readdata(filecol.plotfiles, regexp = regexp ,  substr = substr, filename = filename)

  def generate_plot(self, xlimg = None , ylimg =None , exname = '' , prefix = True , save = True):
    """
    some nice plots to visualize the data with matplotlib, plotg = true if you plot the energylevels of the sp levels of a geometry file
    """
    print ('start with the generation of plots')
    #plot of condensation energy
    self.plotwrap(0,2, 'energy' , name = 'ge'+ exname, titel = 'the energy ', xlim = xlimg , ylim = ylimg , prefix = prefix ,save = save )
    self.plotwrap(0,1, 'condensation energy ' , name = 'ce' + exname ,titel = 'the condensation energy ',xlim = xlimg , ylim = ylimg  , prefix = prefix,save = save )

  def plotwrap(self, xindex, yindex, yas, name = None, titel = None ,color = 'r' , sort = '' , label = None , xlim = None , ylim = None , prefix = False,save = True):
    for i in range(len(self.datarg)):
      self.fig.axes[0].plot(self.datarg[i][:,xindex],self.datarg[i][:,yindex], color+sort , label = label)
      if self.separated == True and save:
        self.layout(self.reader.depvar['depvar'] , yas , tit = titel, xlim = xlim , ylim = ylim)
        self.savefig(name, filenum = i , prefix = prefix)
    if self.separated == False and save:
      self.layout(self.reader.depvar['depvar'] , yas , tit = titel, xlim = xlim , ylim = ylim)
      self.savefig(name + 'together' , prefix = prefix)

  def plotrgvarscplane(self, interval = (-20 , 0), label = None):
    for k in xrange(len(self.datarg) ):
      for j in xrange(len(self.kpunten[filenum])):
        self.layout('real part rgvars ' , 'imaginary part rgvgrs ', tit = 'Richardson-Gaudin variables')
        for i in xrange(self.rgindex,2*self.reader.npair+self.rgindex,2):
          self.fig.axes[0].plot(self.datarg[k][self.kpunten[filenum][j][1] + interval[0]:self.kpunten[k][j][1] + interval[1],i],self.datarg[k][self.kpunten[k][j][1]+interval[0]:self.kpunten[k][j][1] + interval[1],i+1] , 'b' , label = label)
        self.savefig('%f' % (float(self.kpunten[k][j][0])), filenum = k) # you never want this together

  def plotrgvars(self,cplane = False , begin = 0 , stop = None, name = '' , save = True , axnum = 0, xlim = None , ylim = None , prefix = True):
    print('starting to plot the Richardson-Gaudin variables')
    self.plotrgwrap(self.rgindex, 2*self.reader.npair+self.rgindex , self.reader.depvar['depvar'] , 'real part rgvars ',axnum = axnum ,tit =  'Richardson-Gaudin variables', name = 're'+ name , begin = begin , stop =stop , save = save, xlim = xlim , ylim = ylim, prefix = prefix)
    self.plotrgwrap(self.rgindex+1, 2*self.reader.npair+self.rgindex+1 , self.reader.depvar['depvar'] ,'imaginary part rgvars ',axnum = axnum , tit = 'Richardson-Gaudin variables', name = 'im'+ name, begin = begin , stop = stop  , save = save, xlim = xlim , ylim = ylim, prefix = prefix)
    if cplane:
      self.plotrgwrap(self.rgindex, 2*self.reader.npair+self.rgindex ,'real part rgvars ' ,'imaginary part rgvars ',axnum = axnum ,tit =  'Richardson-Gaudin variables', name = 'cp' + name, begin= begin , stop = stop , save = save, xlim = xlim , ylim = ylim, prefix = prefix)

  def plotrgwrap(self, columnstart ,columnend  ,xas , yas , axnum = 0 ,tit = None , begin = 0  , stop = None, name = '' , color = 'b' , sort = '-' ,label = None , save = True , xlim = None , ylim = None, prefix = True):
    for j in xrange(len(self.datarg)):
      #self.plotstar( number = 6 , length = 1  , sort = ['dashed', 'dashdot' ] , color = ['b','r']) #used to create the nice star plot in my factorisable interaction paper
      for i in xrange(columnstart,columnend,2):
        if 'cp' in name:
          """
          sort = ':'
          if j % 2 == 1:
            color = 'r'
            if i == columnstart:
              label = r'$g > \frac{-1}{7}$'
            else:
              label = None
          else:
            sort = '-'
            color = 'b'
            if i == columnstart:
              label = r'$g < \frac{-1}{7}$' #to create the legend uncomment the automatical legend line in the layout
            else:
              label = None
            """
          sort = '.'
          self.fig.axes[axnum].plot(self.datarg[j][begin:stop,i],self.datarg[j][begin:stop,i+1], color+sort , label = label , markersize = 3)#, mfc = 'None')
        else:
          self.fig.axes[axnum].plot(self.datarg[j][begin:stop,0],self.datarg[j][begin:stop,i] , color, label = label)
      if self.separated == True and save:
        self.layout(xas , yas , tit = tit, xlim = xlim , ylim = ylim)
        self.savefig(name , filenum = j, prefix = prefix)
    if self.separated == False and save:
      self.layout(xas , yas , tit = tit, xlim = xlim , ylim = ylim)
      self.savefig(name + 'together', prefix = prefix)
    
  def plotstar(self, number = 6 , length = 2  , sort =[ 'dashed', 'dashdot'], color = ['b','r']):
    colorv  = color[0] ; sortv = sort[0]
    for shift in [0, math.pi/number]:
      for angle in [math.pi *2/6. * i + shift for i in range(number)]:
        x = [0, math.cos(angle) * length]
        y = [0, math.sin(angle) * length]
        self.fig.axes[0].plot(x,y, colorv , linestyle = sortv)
      colorv = color[1] ; sortv = sort[1]
    print 'plotted star'

  def plotintofmotion(self,name = 'iom',stop =None,begin = 0 , xlim = None , ylim = None , samedir = False , colormap = None, axbg = None):
    columns = self.rgindex + 2 * self.reader.npair
    if colormap != None:
      cm = pl.cm.get_cmap(colormap)
      normf = self.normalizefunction([ dat[begin,2] for dat in self.datarg ])   
    for j in xrange(len(self.datarg)):
      for i in xrange(columns,self.reader.nlevel+ columns):
        lines = self.fig.axes[0].plot(self.datarg[j][begin:stop,0],self.datarg[j][begin:stop,i] , c = 'b')
        if colormap != None:
          pl.setp(lines, color = cm(normf(self.datarg[j][begin,2])))
      if self.separated == True:
        self.layout(self.reader.depvar['depvar'] , 'integrals of motion ', tit = 'integrals of motion of the Richardson-Gaudin model' , xlim = xlim , ylim = ylim)
        self.savefig(name , filenum = j)
    if self.separated == False:
      self.layout(self.reader.depvar['depvar'] , 'integrals of motion ', tit = 'integrals of motion of the Richardson-Gaudin model' , xlim = xlim , ylim = ylim , axbg = axbg)
      if colormap != None:
        sm = pl.cm.ScalarMappable(cmap= 'hot', norm=pl.normalize(vmin=0, vmax=1))
        # fake up the array of the scalar mappable. Urgh...
        sm._A = []
        pl.colorbar(sm)
      self.savefig(name + 'together' , samedir = samedir)

  def perezlattice(self, xlim = None , ylim = None , name = 'perezl'):
    """
    If the datafiles of all the states are read in with the right regexp, kpunten contains all the list indices of self.datarg  of the interaction constants of interest
    use the following if you want the perezlattice at particular g (if used regexp to extract info in self.kpunten): to determine the index number in self.datarg: self.kpunten[k][j][1] remark: if this is done the problem exists that the particular g value searched for doens't exist because of change in stepwidht because of circumvention of critical points
    REMARK: at the moment the integrals are colorcoded: full green is integral of motion corresponding to lowest sp level, full red is integral of motion corresponding to the highest sp level, in between is the transition
    """
    colstart = self.rgindex + 2 * self.reader.npair ; nplots = 300 
    for j in range(nplots):
      rowindex = int(len(self.datarg[0])/float(nplots))*j + 7 ; intc = self.datarg[0][rowindex,0]
      for i in xrange(colstart,self.reader.nlevel+ colstart):
        for k in range(len(self.datarg)):
          try:
            lines = self.fig.axes[0].plot(self.datarg[k][rowindex,2],self.datarg[k][rowindex,i] , c = ((i-colstart)/float(self.reader.nlevel),1- (i-colstart)/float(self.reader.nlevel),0), marker = '.')
          except IndexError:
            pass
      self.layout( 'Energy ', 'integrals of motion ', tit = 'integrals of motion of the Richardson-Gaudin model at g = %s' % intc , xlim = xlim , ylim = ylim)
      namek = str(intc).translate(None,'.-') + name 
      self.savefig(namek , filenum = j , prefix = False)
      #print 'we plotted a perez lattice around %s' %self.kpunten[0][j][0]

  def normalizefunction(self , values):
    """
    normalizes the values between 0 and 1
    """
    maxv = np.max(values)
    minv = np.min(values)
    def f(x):
      return (x - minv)/(maxv-minv)
    return f

  def scatterplot(self ,  xvars , yvars , colorvars , colormap = 'hot' ):
      cm = pl.cm.get_cmap(colormap)
      sc = self.fig.axes[0].scatter(xvars ,yvars, c=colorvars, cmap = cm )
      pl.colorbar(sc)

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
      self.datarg[i][:,1:3] = self.datarg[i][:,1:3] - gronddat[:,1:3] #we made sure that the data of the groundstateenergy is first in the rgdata list
    del(self.datarg[0], self.prefix[0])

  def slow_butsure_normalization(self):
    print('Warning we normalize all the excited states to the groundstate energy')
    gronddat =  self.datarg[0][:,2] 
    depvals = list(self.datarg[0][:,0] )
    for i in range(1,len(self.datarg)):
      j = 0 ; gj = 0 ; end = np.shape(self.datarg[i])[0]
      while j < end :
        if depvals[gj] != self.datarg[i][j,0]:
          try:
            gj = depvals.index(self.datarg[i][j,0])
          except ValueError:
            self.datarg[i] = np.delete(self.datarg[i],j, axis=0)
            end -= 1
            print 'skipped some non-matching values for the normalization with the ground-state'
            continue
        self.datarg[i][j,2] = self.datarg[i][j,2] - gronddat[gj] #we made sure that the data of the groundstate-energy is first in the rgdata list
        j += 1 ; gj += 1
    del(self.datarg[0], self.prefix[0])

  def layout(self ,  xlab , ylab , xlim = None , ylim = None , tit = None , axnum = 0 , legendhand = None , legendlab = None , legendpos = 'best' , finetuning = False , axbg = None , fs = 22, ticksize = 10):
    """
    In this function we finetune some aspects of the axes for all the tuning possibilitys see: http://matplotlib.org/api/axes_api.html
    especially the set functions ;)
    """
    print('We are starting with the layout')
    self.fig.axes[axnum].set_xlabel(xlab, fontsize = fs)
    self.fig.axes[axnum].set_ylabel(ylab , fontsize = fs)
    if xlim != None:
      self.fig.axes[axnum].set_xlim(xlim) #good value for xlim in the case of a xi path is : (2*self.rgeq.energiel[0]-5*(self.rgeq.energiel[1]-self.rgeq.energiel[0]),2*self.rgeq.energiel[-1]+0.5) 
    if ylim != None:
      self.fig.axes[axnum].set_ylim(ylim)
    if tit != None:
      self.fig.axes[axnum].set_title(tit , fontsize = fs)
    if legendlab != None:
      self.fig.axes[axnum].legend(legendhand , legendlab, loc = legendpos)  #if you want to add extra info
    
    #self.fig.axes[axnum].ticklabel_format(style='sci', axis='y') #force scientifique notation for y axis
    #self.fig.axes[axnum].yaxis.major.formatter.set_powerlimits((0,0))
    for tick in self.fig.axes[axnum].xaxis.get_major_ticks():
      tick.label.set_fontsize(ticksize) 
    for tick in self.fig.axes[axnum].yaxis.get_major_ticks():
      tick.label.set_fontsize(ticksize) 

    leg = self.fig.axes[axnum].legend(loc = legendpos) #draws the legend on axes[axnum] all the plots that you labeled are now depicted in legend

    if axbg != None:
      self.fig.axes[axnum].set_axis_bgcolor(axbg)
    """
    if you forgot to add a label to a line with linenumber: lnum you can do: self.fig.axes[axnum].lines[lnum].set_label('this is my new label')
    the underneath is the same as : h , l = self.fig.axes[axnum].get_legend_handles_labels()
                                    self.fig.axes[axnum].legend(h,l)
    """
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

  def savefig(self , name , filenum = 0 , samedir = False , prefix = True):
    """
    After we are satisfied with our figure we save it with this function: dpi = pixels per inch, under a name determined by the savestring function().
    """
    #REMARK watch out with the translation of the dot to nothing when you gave as arguments the current working directory '.' because
    #if you do this it is not possible to save the file in the appropriate place because the folder doesn't exist anymore 
    #because the first . dissapeared you can only remove . from floats or extensions not from current dir (maybe build in check that if the first letter of the filename is a dot then that dot is not removed)
    figname = self.savestring(name , filenum , samedir = samedir , prefix = prefix )
    self.fig.savefig(figname , dpi = 80 , facecolor = 'w' , edgecolor = 'w')
    self.fig.clf()
    self.fig.add_subplot(111 , axisbg = 'w' , projection = 'rectilinear')

  def savestring(self , name , filenum , samedir = False , prefix = True):
    """
    This function generates the name whereunder the figure is going to be saved
    """
    if prefix == True:
      if samedir:
        """
        Making use of some implementation detail of savefig, if we read in files from all different directory's, the prefixes contain the path of those files relative to the rootdirectory. So if you save the file we save it with first the prefix and then the name , so the figures end up in the same directory as the files. If you don't want this behaviour we need to remove the / in the prefixs so fig.savefig will not recognize it as a path so all the figures end up in the current working directory. Remark we only remove the / because if all the figures end up in same dir we need the path information to distinguish them.
        """
        self.prefix = [pre.translate(None , '/.')  for pre  in self.prefix]
      return '%s%s%d.png' %(self.prefix[filenum], name, filenum )
      #return '%s%s.png' %(self.prefix[filenum], name)
    else:
      return '%s%d.png' %(name, filenum )

  def writetext(self ,text , pos , axnum = 0, hor = None  ,ver = None , rot = None ,fs =14 , transform = None):
    self.fig.axes[axnum].text(pos[0] ,pos[1] ,text , rotation = rot ,horizontalalignment = hor, verticalalignment = ver , fontsize = fs, transform = transform) #, color = 'black', style = 'italic')

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
      #cplane eventjes op True gezet
      self.plotrgvars(cplane = True , begin = 0 , stop = None)
    if intm:
      self.plotintofmotion()

class Plot_Geo_File(Plot_RG_Files):
  """
  remark before the Richardson-Gaudin variables start this file has 6 columns, extra: nig , meandistance , number of levels
  """
  def __init__(self , name = 'x', searchstr = 'plotenergy'):
    self.rgindex = 6
    if os.path.isdir(name):
      self.procesfiles(name , searchstr)
    elif os.path.isfile(name):
      self.readdata([name])
    super(Plot_Geo_File,self).__init__()

  def generate_plot(self):
    super(Plot_Geo_File,self).generate_plot()
    print('plot non-interacting groundstate')
    self.plotwrap(0,3, 'energy of the non-interacting groundstate ','nig', titel = 'aantal paren = %f' %(self.reader.npair))
    try:
      self.plotwrap(0,4,"d " ,'meandistance', titel = "number of sp levels = %f" %self.reader.nlevel)
    except:
      print 'the plot of d failed'

class Plot_Data_File(Plot_RG_Files ):
  def __init__(self, name = 'x', searchstr = 'plotenergy' , notsearch =r'\.swp|\.png' , regexp = None, sortfunction = None):
    self.rgindex = 3
    if os.path.isdir(name):
      self.procesfiles(name , searchstr, notsearch = notsearch , regexp = regexp , sortfunction = sortfunction)
    elif os.path.isfile(name):
      self.readdata([name], regexp = regexp)
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

  def plotrgcloud(self ,begin = 0, step = 1 ,  colormap = 'hot' , xlim = None , ylim = None):
    while begin <= np.shape(self.datarg[0])[0]:
      revars = [rerg for dat in self.datarg for rerg in dat[begin,self.rgindex:self.rgindex+2*self.reader.npair:2]]
      imvars = [imrg for dat in self.datarg for imrg in dat[begin,self.rgindex+1:self.rgindex+2*self.reader.npair + 1:2]]
      energy = [[dat[begin,2]] *self.reader.npair for dat in self.datarg ]
      self.scatterplot(revars , imvars , energy , colormap = colormap)
      self.layout( 'real part of rgvars ' ,  'imaginary part of rgvars ', xlim = xlim , ylim = ylim, tit = 'RG vars g = %f all states'%(self.datarg[0][begin , 0]) , axnum = 0 , legendhand = None , legendlab = None , legendpos = 'best' , finetuning = False)
      self.savefig('allstates%f' % (self.datarg[0][begin,0]) , samedir = True)
      begin += step
    makemovie(name = 'allstatesrgcloud')    
  
  def plot_spectrum(self,xlim = None , ylim = None, search = 'plotenergy', rgw = True, intm = True, name = 'spectrum', readgreen = False, standard = True, save = True):
    self.procesfiles(os.getcwd(), search, notsearch = r'\.swp|\.png|readgreen', sortfunction = lambda x : -1. if '/0/' in x or 'ground' in x or 'grond' in x else 0.) #sortfunction makes sure the groundstate is first this is important for the normalization
    if standard: self.standard_plot(rgw , intm)
    #self.normalize_to_groundstate()
    self.slow_butsure_normalization()
    self.separated = False
    for data in self.datarg:
      data[:,2] /= (self.reader.npair*4.)**2.
      data[:,0] += 1/(self.reader.npair*2.+2)
    if readgreen: mine , pre = self.readgreen()                                                       
    self.generate_plot(xlimg = xlim , ylimg = ylim, prefix = False , exname = name, save = save)
    if readgreen: print mine , ' prefix is :', pre

  def plot_gapsurvey(self,xlim = None , ylim = None, search = 'plotenergy', rgw = True, intm = True, name = 'spectrum', readgreen = True, standard = True, save = False, dir = '.'):
    for i in os.listdir(dir):
      if os.path.isdir(i):
        os.chdir(i)
        self.plot_spectrum(xlim = xlim, ylim = ylim, search = search, rgw = rgw, intm = intm, name = name , readgreen = readgreen, standard = standard , save = save )
        os.chdir('..')
    for pos,text in [((-0.039,0.08),'6p' ) ,((-0.03,0.08),'10p'),((-0.021,0.08),'15p'),((-0.015,0.08),'25p'),((-0.0065,0.08),'40p')]:
      self.writetext(text , pos ,axnum = 0, hor = 'left', ver = 'bottom', rot = 0,fs =14,transform = self.fig.axes[0].transData)

    self.layout( r'$g+ \frac{1}{2N+2}$' ,r'$\frac{E_{ex}}{(4N)^2}$'  , tit = 'exploring the gap', xlim = xlim , ylim = ylim)
    self.savefig( 'master' , filenum = 0 , samedir = False , prefix = False)

  def readgreen(self):
    readgreenpoint = self.reader.eta/(2.*(self.reader.npair-1) -2.*np.sum(np.array(self.reader.degeneracies)/4. -np.array(self.reader.seniorities)))
    step = abs(self.datarg[0][0,0] -self.datarg[0][1,0])  
    self.fig.axes[0].axvline(x = readgreenpoint,ymin = 0 , ymax = 1, c = 'b', linewidth = 1)
    datareadgreen = [row[2]  for data in self.datarg for row in data if ((row[0] - step/2. <= readgreenpoint) and (row[0] + step/2.  >= readgreenpoint ) )] 
    #datareadgreen = []
    #for data in self.datarg:
    #  f = False
    #  for row in data:
    #    if ((row[0] - step/2.< readgreenpoint) and (row[0] + step/2. > readgreenpoint ) ):
    #      datareadgreen.append(row[2])
    #      f = True
    #      break
    #  if f == False:
    #    datareadgreen.append(1e9)
    #assert(len(datareadgreen)  == len(self.prefix))
    lowest = min(datareadgreen) 
    file = open('readgreenfile.dat' , 'w')
    file.write('#the read-green point is at: %f \n#the value of the lowest excited state is: %f\n#the filename of this state is: %s \n' %(readgreenpoint ,lowest, self.prefix[datareadgreen.index(lowest)] ))
    file.close()
    return lowest,  self.prefix[datareadgreen.index(lowest)]     

  
class Plot_Xi_File(Plot_RG_Files):
  def __init__(self, name , search , regexp =r'constant:\s*([\-0-9.]+)'):
    self.rgindex = 2
    if os.path.isdir(name):
      self.procesfiles(name ,search, notsearch =r'\.swp|\.png', regexp = regexp)
    elif os.path.isfile(name):
      self.readdata([name], regexp = regexp)
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
    self.layout(self.reader.depvar['depvar'] , r'energy spectrum ' , tit = r'All tda start distributions $\xi$' , legendhand = [goodline , badline] , legendlab = ['solution', 'breakdown'] )
    self.savefig('xispec')

  def plotrgvarsxi(self, name = 'rgvxi' ,xlim = None , ylim = None):
    for j in xrange(len(self.datarg)):
      for i in np.arange(self.rgindex,2*self.reader.npair+self.rgindex,2):
        self.fig.axes[0].plot(self.datarg[j][0,i],self.datarg[j][0,i+1],'b.', markersize = 23) #Richardson-Gaudin solutions (xi = 1)
        self.fig.axes[0].plot(self.datarg[j][len(self.datarg[j][:,0])-1,i],self.datarg[j][len(self.datarg[j][:,0])-1,i+1],'b.',mfc = 'None', markersize = 23) # Corresponding tda solutions (xi = 0 )
        self.fig.axes[0].plot(self.datarg[j][:,i],self.datarg[j][:,i+1],'b-' , lw =2) # intermediate values of xi
      if self.reader.eta == None:
        sing = np.array(self.reader.elevels)* 2
      else:
        sing = self.reader.eta * np.array(self.reader.elevels) * np.array(self.reader.elevels)
      for i in range(3):#self.reader.nlevel):
        self.fig.axes[0].axvline(x = sing[i] ,c=  'k',linestyle = '--')
      if self.separated == True:
        self.layout('real part of rgvars ', 'imaginary part of rgvars ', xlim =xlim, ylim = ylim, tit = 'g = %s ' %(self.kpunten[j][0][0]) , fs = 20)
        self.savefig(name , filenum = j, prefix = False)
    if self.separated == False:
      self.layout('real part of rgvars ', 'imaginary part of rgvars ', xlim =xlim, ylim = ylim, tit = 'g = %s ' %(self.kpunten[j][0][0]) , fs = 20)
      self.savefig(name + 'together' , prefix = False )

class Plot_All_File(Plot_RG_Files):
  def __init__(self,name, g , regexp = r'seniority\s\[(.+?)\]',substr = r'\{.*\}'):
    self.chardata = g
    self.rgindex = 2
    super(Plot_All_File,self).__init__()
    self.readdata(name, regexp = regexp ,substr = substr)

  def plotrgcloud(self):
    """
    This function needs it own datareader because it's to specific
    """
    import itertools
    print self.kpunten
    for i in range(len(self.kpunten[0])):
      self.writetext('sen ='+ self.kpunten[0][i][0], (0.60,0.95), axnum = 0, hor = 'left',ver = 'bottom', rot = None ,fs =14 , transform = self.fig.axes[0].transAxes)
      if i == len(self.kpunten[0]) -1 :
        end = None
      else:
        end = self.kpunten[0][i+1][1] 
      begin = self.kpunten[0][i][1] 
      revars = list(itertools.chain.from_iterable(self.datarg[0][begin:end,self.rgindex:self.rgindex+2*self.reader.npair:2]))
      imvars = list(itertools.chain.from_iterable(self.datarg[0][begin:end,self.rgindex+1:self.rgindex+2*self.reader.npair + 1:2]))
      energy = list(itertools.chain.from_iterable([[line[0]] * self.reader.npair for line in self.datarg[0][begin:end]]))
      self.scatterplot(revars , imvars , energy , colormap = 'hot')
      self.layout( 'real part of rgvars ' ,  'imaginary part of rgvars ', xlim = None, ylim =None, tit = 'RG vars g = %f all states'%(self.reader.g) , axnum = 0 , legendhand = None , legendlab = None , legendpos = 'best' , finetuning = False)
      self.savefig('allstates%f' % (self.reader.g) , samedir = True)
  
def plot(name, ylim = None):
  #plots a two column file column 1 is x-axis, column 2 is y-axis
  for file in os.listdir(os.getcwd()):
    print file
    if name in file and '.dat' in file:
      data = np.loadtxt(file)
      if 'FCI' in file:
        print file
        data[:,0] *= 1./0.52917721092
        pl.plot(data[:1500,0] ,data[:1500,1] , label = file.strip('.dat').strip(name))
      else:
        pl.plot(data[:100,0] ,data[:100,1] , label = file.strip('.dat').strip(name))

  if ylim != None: pl.ylim( ylim)
  pl.legend()
  pl.savefig('%s.png' %name.strip('.dat'))

def main(option, args):
  plotter = Plot_Data_File()
  plottergeo = Plot_Geo_File()
  if option == 'pexcited':
    plotter.plot_spectrum(xlim = (-0.3,0), ylim = (0,1000), search = 'plotenergy', rgw = True, intm = True, name = 'spectrum' , readgreen = True, standard = False, save = False)

  if option == 'gapsurvey':
    plotter.plot_gapsurvey(xlim = (-0.04,0.04), ylim = (0,0.1), search = 'plotenergy', rgw = False, intm = True, name = 'spectrum' , readgreen = False, standard = False, save = False,dir = '.')

  if option == 'wpairing':
    if args[1] == True:
      plottergeo.procesfiles(args[0], 'plotenergy')
      plottergeo.generate_plot()
    else:
      plotter.readdata(args[0])
      #plotter.procesfiles(args[0],'plotenergy')
      plotter.generate_plot(xlimg = None, ylimg = None)

  if option == 'inset':
    """
    Example of how you need to draw a small inset in a larger plot of a particular
    area of interest with matplotlib
    """
    plotter.readdata([args[0]])
    plotter.reader.depvar['depvar'] = r'$\eta$ '  #change the future x-axis label to latex 
    begin =0
    stop = None
    plotter.plotrgvars(begin = begin , stop = stop , name = 'etanul2', save = False)
    begin =  9880
    stop = None
    plotter.rgindex = 5
    plotter.reader.npair = 9
    plotter.add_axes([0.5,0.2,0.3,0.3])
    plotter.fig.axes[1].xaxis.set_major_locator(matplotlib.ticker.LinearLocator(5))
    #see also: http://scipy-lectures.github.io/intro/matplotlib/matplotlib.html#ticks
    plotter.plotrgvars(begin = begin , stop =stop, axnum = 1)

  if option == 'rgclouddata':
    plotter.procesfiles(args[0] , 'plotenergy' , notdir = 'movie')
    plotter.plotrgcloud(step = 100, xlim = (-14.,-7.) , ylim = (-2.5, 2.5))
      
  if option == 'addlevel':
    plotter.procesfiles( '.' , 'plotenergy' ,   sortfunction = lambda s : int(re.search(r'\d+' , s).group()), rev = True , regexp = r'^%f\s+[\-+\.\d]+\s+([\-+\.\d]+)\s' % args[0])
    plotter.addlevel(args[0])
    
  if option == 'rgvar':
    ref = args[0] 
    begin =0
    stop = None
    cp = args[1]
    plotter.procesfiles(args[0],'plotenergy',filename = False)
    #plotter.readdata(args[0])
    plotter.reader.depvar['depvar'] = 'g'  #change the future x-axis label to latex 
    plotter.separated = True
    plotter.plotrgvars(cplane = cp , begin = begin , stop = stop , name = '', xlim = None, ylim = None, prefix = True)
  
  if option is 'rgcloud':
    name = 'pairingtinsen=0.dat'
    plottera = Plot_All_File(name, -0.217 , regexp = r'seniority\s\[(.+?)\]',substr = r'\{.*\}')
    plottera.plotrgcloud()
  
  if option is 'cprgvar':
    ref = args[0]
    plotter.readdata([ref], regexp = r'^#\s+((-|\d)\d+\.*\d*)\s+kritisch', begin = 1)
    plotter.plotrgvarscplane(interval = (-20,0))

  if option is 'intmotion':
    #plotter.readdata([args])
    plotter.separated = True
    if args[3] != 'perez':
      plotter = Plot_Data_File(args[0] , args[1] , args[2])
      plotter.plotintofmotion(name = 'intmotion',xlim = (-1.5,0.), ylim = (-2 , 2) , samedir =True , colormap ='hot' , axbg = 'g')
    else:
      plotter = Plot_Data_File(args[0] , args[1] , args[2])# , regexp = r'^(-0.003100|-0.063100|-0.123100|-0.213100|-0.363100|-0.993100)')
      plotter.perezlattice()

  if 'xi' in option:
    plotterxi = Plot_Xi_File(args[0], args[1], regexp = r'constant:\s*([\-0-9.]+)')
    if option is 'xipath':
      #plotterxi.procesfiles(args[0],args[1] , regexp = r'constant:\s*([\-0-9.]+)')
      plotterxi.separated = True
      plotterxi.plotrgvarsxi(ylim = None , xlim = None)

    if option is 'specxichange':
      #to plot entire spectra with broken down in red and states who went from xi = 0 to xi =1 in blue
      plotterxi.plot_spectrumxichange()
  
def defineoptions():
  '''
  possible options: 'pexcited' plots all the excited states relative to the groundstate, 'wpairing' plots the 
  results from a writepairing call in writepairing.py(main), 'addlevel' from a set of outputfiles from generating_datak generated
  by adding empty sp levels and get from those files the groundstate energy at a constant g and plot them and perform lin.regression 
  '''
  #common options are: wpairing, rgvar, intmotion
  option = 'rgcloud'
  #args =  -0.137 , None
  args =  '.','xipath' ,True , 'xipath',True, 'xipath' ,  False,'.',r'constant:\s*([\-0-9.]+)', r'xi[0-9\.a-zA-Z\-]+.dat$','g'  ,False 
  main(option,args)

if __name__ == '__main__':
  defineoptions()
  #makemovie()
