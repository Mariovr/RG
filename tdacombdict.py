#! /bin/env python
from itertools import combinations,combinations_with_replacement
import numpy as np
import sys
import os,re
import pylab as pl
import plotfunctions as pf
import rgfunctions as rg

def collectrgcloud(filepath, g ,collectfile):
  data = rg.readlines(filepath,g)
  collectfile.write('  '.join(map(float,data))) 

  
def printtda(nlev = 12, apar = 10 , ontaardingen = [4,4,4,8,4,4,8,8,4,8,4,8,12] , senioriteit = None, startw = -1,  endw = 310000):
  '''
  generates a file that contains the connection between the wind parameter (the number of the state as determined by combinations with replacement)
  and the distribution of the pairs over the sp levels at low enough g so the sp levels that contain a pair behave as the corresponding sp levels.
  '''
  tdacombinations = combinations_with_replacement(np.arange(0,nlev+1),apar)
  if senioriteit == None:
    senioriteit = np.zeros(nlev)
  file = open('correspondentie.dat', 'w')
  wind = 0
  for tdadict in tdacombinations:
    tdastartd = {}
    goodsol = True
    for i in tdadict:
      a = tdadict.count(i)      
      tdastartd[i] = a
      if a*2 + senioriteit[i]*2 > ontaardingen[i]:
        goodsol = False
    if goodsol == True:
      wind += 1
      if wind >= startw and wind < endw:
        print 'the %g th state corresponds to %s ' %( wind,str(tdastartd))
        file.write('%g\t%s\n' %(wind,str(tdastartd)))
    file.close()

def find(func,rootdir,arg = None):
  #call function for all rgvar files in rootdir
  files = os.listdir(rootdir)
  for file in files:
    filepath = os.path.join(rootdir,file)
    if os.path.islink(filepath):
      pass
    elif os.path.isdir(filepath):
      find(func,filepath,arg)
    elif os.path.isfile(filepath) and 'plotenergy' in filepath: 
      func(filepath,arg[0] , arg[1])
    else:
      print 'find cann\'t read', filepath

def main():
  g = -0.01
  collectfile = open('allstates%d.dat' %g ,'w')
  rootdir = '$HOME/data/eta=1picketfenceneg/'
  collectfile.write('#This file contains the entire spectrum at interaction constant: %f\n#condensatione\tenergy\trgvarreal\trgvarim\trgvarreal\rgvarim... ' %g)
  find(collectrgcloud,rootdir,arg = (g,collectfile ))
  collectfile.close()


if __name__ == '__main__':
  main() 





