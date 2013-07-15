import os,re
import numpy as np
import pylab as pl

#REMARK don't forget to install mencoder on the computer that exectutes this script otherwise no movie
def makemovie():
  # makes a movie from all the .png files in the current directory
  print 'Starting to create a movie, with all the .png files in directory: %s ' %str(os.getcwd())
  dirname = 'moviexivar'
  command = ('mencoder',
           'mf://*.png',
           '-mf',
           'type=png:w=800:h=600:fps=15',
           '-ovc',
           'lavc',
           '-lavcopts',
           'vcodec=mpeg4',
           '-oac',
           'copy',
           '-o',
           dirname+'.avi')

  os.spawnvp(os.P_WAIT, 'mencoder', command)
  
  
def plotrgvarsxi(filepath, argp):
  cwd = os.getcwd()
  plotre = open(filepath, 'r')
  print 'we are at : ' ,filepath
  dirname, namerg = os.path.split(filepath)
  os.chdir(dirname)
  apair = argp
  eendlev = np.arange(0,12)
  print namerg
  match = re.search('-\d+[.]\d+',namerg)
  g = match.group()
  print ' we plot rgvar with g = ' , g
  datare = np.loadtxt(plotre)
  pl.figure()
  irange = np.arange(1,2*apair+1,2)
  for i in irange:
    try:
      pl.plot(datare[0,i],datare[0,i+1],'g.', markersize = 10)
      pl.plot(datare[len(datare[:,0])-1,i],datare[len(datare[:,0])-1,i+1],'r.',mfc = 'None', markersize = 10)
      pl.plot(datare[:,i],datare[:,i+1],'b-')
    except IndexError:
      print 'We encountered an IndexError at %s' %filepath
      pass
  for i in range(len(eendlev)):
    pl.axvline(x = eendlev[i]*2. ,c=  'k',linestyle = '--')
  pl.xlabel('real part of rgvars (a.u)')
  pl.ylabel('imaginary part of rgvars (a.u.)')
  pl.title('Richardson-Gaudin variables at g = %f (xi in [0,1])' %(float(g)))
  pl.xlim((-4,22))
  pl.ylim((-10,10))
  pl.savefig('%s.png' %namerg )
  pl.close()
  os.chdir(cwd)
  
def find(func,rootdir,arg = None):
  #call function for all rgvar files in rootdir
  files = os.listdir(rootdir)
  for file in files:
    filepath = os.path.join(rootdir,file)
    mov = True
    if os.path.islink(filepath):
      pass
    elif os.path.isdir(filepath):
      find(func,filepath,arg)
    elif '.png' in file and 'moviedir' in rootdir:
      mov = False
      break
    elif os.path.isfile(filepath) and 'moviedir' in rootdir and '.png' not in filepath and 'tdafile' not in filepath:
      func(filepath,arg)
    else:
      print 'find cann\'t tread', filepath
  if 'moviedir' in rootdir and mov == True:
    cwd = os.getcwd()
    os.chdir(rootdir)
    makemovie()
    os.chdir(cwd)


if __name__ == '__main__':
  npair = 6
  rootdir = '/home/beheerder/sambataro2'
  find(plotrgvarsxi,rootdir,arg = npair)
  
  
