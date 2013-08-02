#! /bin/env python
from itertools import combinations,combinations_with_replacement
import numpy as np
import sys
import os,re
import pylab as pl
import plotfunctions as pf
import rgfunctions as rg

"""
REMARK this file is totally DEPRECATED from 30/07/2013 on because everything it's able to do except the creation of new files with rearranged data (which is useless because the data is already there) can now be done by the generalized plotting classes of plotfunctions.py (Normally almost all the plotting should be done from that file because it's utterly powerfull). See -> in class Plot_Data_File() the method plotrgcloud()

"""

def extract(string , regexp, group = none):
  regex = re.compile(regexp)
  match = regex.match(string)
  if group == none:
    return match.group(0)
  else:
    return match.group(group)

def collectrgcloud(filepath, g , npair , collectlist):
  print 'we are collecting the information of the rg solutions at one particular g : %f in the list collectlist' %g
  data = rg.readlevels(filepath,g,nauw = 0.003)
  #state = extract(filepath,r'.+(\{.+\})',group = 1)
  try:
    file = open('correspondentie.dat','r')
  except ioerror:
    printtda()
    file = open('correspondentie.dat','r')
  idnum = re.search( r'/(\d+)/.+dat$',filepath).group(1)
  idnum = int(idnum)+1
  for line in file:
    if line[0] == '#':
      continue
    line = line.split('\t')
    if idnum == int(line[0]):
      state = line[1].strip('\n')  
  if len(data) ==  2*npair + 2:
    collectlist.append(state+'  ' +'  '.join(map(str,data))+'\n') 
  file.close()
  
#def printtda(nlev = 12, apar = 10 , ontaardingen = [4,4,4,8,4,4,8,8,4,8,4,8,12] , senioriteit = none, startw = -1,  endw = 310000):
def printtda(nlev = 12, apar = 6 , ontaardingen = [2,2,2,2,2,2,2,2,2,2,2,2] , senioriteit = none, startw = -1,  endw = 310000):
  '''
  generates a file that contains the connection between the wind parameter (the number of the state as determined by combinations with replacement)
  and the distribution of the pairs over the sp levels at low enough g so the sp levels that contain a pair behave as the corresponding sp levels.
  '''
  tdacombinations = combinations_with_replacement(np.arange(0,nlev),apar)
  if senioriteit == none:
    senioriteit = np.zeros(nlev)
  file = open('correspondentie.dat', 'w')
  wind = 0
  for tdadict in tdacombinations:
    tdastartd = {}
    goodsol = true
    for i in tdadict:
      a = tdadict.count(i)      
      tdastartd[i] = a
      if a*2 + senioriteit[i]*2 > ontaardingen[i]:
        goodsol = false
    if goodsol == true:
      wind += 1
      if wind >= startw and wind < endw:
        print 'the %g th state corresponds to %s ' %( wind,str(tdastartd))
        file.write('%g\t%s\n' %(wind,str(tdastartd)))
  file.close()

def plotrgcloud(name,npair,g):
  print 'starting to plot rgcloud at g = %f and from file %s' %(g , name)
  namef = open(name,'r')
  namerg = 'rgcould'+name
  pl.figure()
  start = false
  revar = [] ; imvar = [] ; energie = []
  for line in namef:
    if line[0] is '#':
      continue
    index = line.index('}')
    line = line[index+2:]
    larray = line.split()      
    karray = map(float,larray[2:])
    energie += npair * [float(larray[1])]
    revar  += [ i for i in karray if karray.index(i) %2 == 0]
    imvar += [i for i in karray if karray.index(i) % 2 == 1]
  cm = pl.cm.get_cmap('hot')
  sc = pl.scatter(revar,imvar, c= energie,cmap =cm)       
  pl.colorbar(sc)
  pl.xlabel('real part of rgvars (a.u)')
  pl.ylabel('imaginary part of rgvars (a.u.)')
  pl.title('rg variables at g = %f for all ex. states' %(g))
  #ax = pl.gca()
  #pl.text(0.65,0.85,'sen ='+ str(sen),transform = ax.transaxes)
  pl.savefig('%s.png' %(namerg))
  pl.close()

def find(func,rootdir,arg = none):
  #call function for all rgvar files in rootdir
  print 'entering find function with %s as rootdir ,  %s as the function that is going to be executed' %(rootdir, func)
  files = os.listdir(rootdir)
  for file in files:
    filepath = os.path.join(rootdir,file)
    if os.path.islink(filepath):
      pass
    elif os.path.isdir(filepath) and not 'moviedir' in filepath:
      find(func,filepath,arg)
    elif os.path.isfile(filepath) and 'plotenergy' in filepath and '.png' not in filepath: 
      func(filepath,arg[0] , arg[1],arg[2])
    else:
      pass

def sortfunction(obj):
  d = obj.split('}') 
  e = d[1]
  f = e.split()
  return float(f[1])

def main():
  g = -0.14410 ; npair = 6
  while g >= -1.:
    collectlist = []
    rootdir = os.path.join(os.environ['home'], 'data/eta1picketfenceneg/')
    find(collectrgcloud,rootdir,arg = (g,npair,collectlist ))
    collectlist = sorted(collectlist , key = sortfunction ) #sort collectlist according to sortfunction
    #print collected information in collectlist to file
    collectfile = open('allstates%f.dat' %g ,'w')
    collectfile.write('#this file contains the entire spectrum at interaction constant: %f\n#cenergy\tenergy\trgvarreal\trgvarim\trgvarreal\trgvarim... \n' %g)
    for i in collectlist:
      collectfile.write(i)
    collectfile.close()
    plotrgcloud('allstates%f.dat' %g , npair,g )
    g -= 0.003
  rg.makemovie()


if __name__ == '__main__':
  main() 
  #printtda()





