# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Mario Van Raemdonck, 2013;
# (c) Ghent University, 2013
# -*- coding: utf-8 -*-
#! /bin/env python
import sys,os
import numpy as np
import pylab as pl
from matplotlib.font_manager import FontProperties
from itertools import permutations , combinations

import tdasolver as td
import rgfunctions as rf
import richardsongaudin as rg
import datareader as dr

"""
Remark this file calculates overlaps in a general fashion suitable for tdastates, rich. red. bcs. states, and
factorisable interaction states. This method is very slow at the moment, because we just calculate permanents withnested for loops which is incredibly slow. Faster method is to make use of the Muir-Brochard theorem or the formulas of Caux et. al. The file correlationfaribault contains an implementation of those formulas for the reduced BCS Hamiltonian. TODO extend that files with the formulas for the factorisable interaction Hamiltonian so we are also able to calculate fast overlaps and correlationcoefficients for the factorisable interaction hamiltonian.
"""

class State_Calculator(object):
  """
  class that calculates the state by making use of the calc_coef(self, tdasol) that every instance of equation need to have
  """
  def __init__(self,equation,npairs,tdadict = None, normalize = True):
    self.eqvars = equation
    self.npair = npairs
    if tdadict == None:
      if isinstance(equation,td.TdaSolver):
        self.tdadict = self.eqvars.tdadict
      else:
        self.tdadict =  dict(zip(range(self.npair) , [1] * self.npair))
    else:
      self.tdadict = tdadict
    self.state = self.construct_state(normalize = normalize)

  def set_tdadict(self,tdadict):
    self.eqvars.set_tdadict(tdadict)

  def calc_coeflist(self):
    oplossingen = self.eqvars.get_solutions()
    calclist = []
    for i in self.tdadict.keys():
      if self.tdadict[i] != 0:
        for k in range(self.tdadict[i]): #only important for tdastates because then more then one pair can fit in a level, in case of RG variables self.tdadict[i] =1 always
          calclist.append( self.eqvars.calc_coef(oplossingen[i])) #calc_coef returns an array with [1/(2e0 - soli) , 1/(2e1-soli) , ...]
    return calclist

  def calc_statecoef(self,state , coeflist ):
    """
    calculates the coefficient of a sp state (state) in the tdastate with apair pairs making use of coeflist which contains  arrays of calcxicoefficients with constant  sol but changing sp level as in the summation for the creation of the state
    """
    perms = permutations(state)
    totcoef = 0.
    for perm in perms:
      i = 0
      coef = 1.
      for coefs in coeflist:
        coef *= coefs[perm[i]]
        #print coef , perm , coefs[perm[i]]
        i += 1
      totcoef += coef
    #print totcoef
    return totcoef

  def construct_state(self, normalize = True):
    """
    Constructs the vector with the number of pairs that self.dict defines.
    TO DO generalize for general degeneracies
    """
    nlevel = len(self.eqvars.energiel)
    coeflist = self.calc_coeflist()
    assert self.npair < sum(self.eqvars.ontaardingen), 'number of pairs is bigger then the state can contain'
    assert sum(self.eqvars.ontaardingen) == 2*len(self.eqvars.energiel), 'the degeneracy is larger then 2 please update the construction of the state to general degeneracys before moving along'
    #picklist = []
    #for i in range(nlevel):
    #  picklist += [i] * self.eqvars.ontaardingen[i]/2
    spstates = combinations( np.arange(nlevel),self.npair)
    dim = rf.binoml(nlevel,self.npair)
    state = np.zeros(dim,complex)
    i = 0
    for state_id in spstates:
      state[i] = self.calc_statecoef(state_id, coeflist)  
      i += 1 
    assert i == dim , 'i = %g , dim-1 = %g' %(i , dim-1) 
    if normalize:
      state = state / np.linalg.norm(state)
    return state

  def calc_overlap(self, state):
    print 'calculating overlap '
    return np.dot(self.state , state.state).real

  def get_norm(self):
    return np.linalg.norm(self.state)


def create_tdafile(rgeq , dvar , begin = -0.0001 , end = -2. , step = -0.001 , filen = 'alltdas.dat'):
  rsolve = rg.RichardsonSolver(rgeq) # quick hack to let the richardson solver class generate the appropriate tdasolver 
  tda = rsolve.tda
  file = open(filen , 'w')
  file.write(str(tda))
  file.write('#%s\ttda1 .... tdan\n' %dvar)
  while(begin > end):
    tda.setvar(dvar, begin)
    tda.bisect_tda()
    file.write('%f %s\n' %(tda.getvar(dvar), '  '.join(map(str,tda.get_solutions()))))
    begin += step
  file.close()

def list_to_dict(tdalist, nlevel):
  return dict(zip(np.arange(nlevel), tdalist))

def plotoverlaps(tdadictlist , cvar , filen = 'overlap.dat', legend = True , xlim = None , ylim = None, fs = 17):
  plotf = open(filen, 'r')
  cvar = '|'+cvar + '|'
  try:
    dataf = np.loadtxt(plotf,comments = '#')
  except:
    dataf = np.array([])
    plotf.seek(0) #go back to beginning of file
    for line in plotf:
      if line[0] == '#':
        continue
      pline = np.array(map(float,line.split()))
      if len(dataf) <= 1:
        dataf = pline
      else:
        try:
          dataf = np.vstack((dataf,pline))
        except:
          continue
    print np.shape(dataf)
  pl.figure()
  ax = pl.subplot(111)
  for i in xrange(len(dataf[0,:])-1):
    if legend == True:
      ax.plot(abs(dataf[:,0]),abs(dataf[:,i+1]), label = str(tdadictlist[i]))
    else:
      ax.plot(np.log(abs(dataf[:,0])),abs(dataf[:,i+1]))
  if legend == True:
    writelegend(ax , cvar)
  else:
    writetext(ax)
  pl.ylabel('Overlaps RG ground-state with tdastates', fontsize = fs)
  pl.xlabel( 'log|g|' , fontsize = fs)
  #pl.title('Overlaps tdastates with Richardson-Gaudin ground-state', fontsize = fs)
  pl.xlim(xlim)
  pl.ylim(ylim)
  pl.savefig('%s2.png' %filen.strip('.dat'))
  pl.close()

def writelegend(ax):
  box = ax.get_position()
  ax.set_position([box.x0 - box.width * 0.05 , box.y0 , box.width * 0.80 , box.height])
  fontL = FontProperties()
  fontL.set_size('x-small')
  ax.legend(loc = 2 , prop = fontL , bbox_to_anchor = (1.00,1.00) , fancybox = True , shadow = True)

def writetext(ax , fs= 12):
  ax.text(-4.501,0.95,r'(1111110$\ldots$0)' , horizontalalignment = 'left', fontsize = fs)
  ax.text(-4.501,0.75,r'(111120$\ldots$0)' , rotation = 10, verticalalignment = 'bottom' , horizontalalignment = 'left', fontsize = fs)
  ax.text(-4.5011,0.54,r'(11130$\ldots$0)' , verticalalignment = 'bottom' , horizontalalignment = 'left', rotation =20 ,fontsize = fs)
  ax.text(-4.5011,0.40,r'(211110$\ldots$0)' , rotation = 20,verticalalignment = 'bottom' , horizontalalignment = 'left', fontsize = fs)
  ax.text(0.91,0.65,r'(60$\ldots$0)' , horizontalalignment = 'right', fontsize = fs)
  ax.text(0.8,0.42,r'(510$\ldots$0)' , horizontalalignment = 'left', fontsize = fs)
  ax.text(-1.30,0.90,r'(150$\ldots$0)' , rotation =-60 , fontsize = fs , horizontalalignment = 'left')
  ax.text(1.43,0.25,r'(1140$\ldots$0)' , rotation =-15 ,horizontalalignment = 'right', fontsize = fs)
  ax.text(0.910, 0.20,r'(4110$\ldots$0)' , rotation = -20 ,horizontalalignment = 'right', fontsize = fs)
  ax.text(-6.9011,0.39,r'(31110$\ldots$0)' ,horizontalalignment = 'left', fontsize = fs)
  ax.text(1.43 ,0.123 ,r'(21120$\ldots$0)' , rotation = -8 ,horizontalalignment = 'right', fontsize = fs) #, color = 'black', style = 'italic')

def main():
  #create_tdafile(rgeq , 'g' ,filen = 'alltdasredbcspicketfence.dat')
  linen = 5
  linestep = 5
  filen = 'plotenergy.dat'
  outputf = 'overlapfac5ti.dat' 
  outfile = open(outputf, 'w')
  ##tdadictlist = [{0:1 , 1:1 , 2:1 , 3:1 ,4:1 , 5:1}, {0:6} , {0:1 , 1:5} , {0:5 , 1:1},{0:1 , 1:1 , 2:4} ,{0:1 , 1:1 , 2:1 , 3:1 ,4:2},{0:4 , 1:1 , 2:1}, {0:2 , 1:1 , 2:1 , 3: 2}, {0:2 , 1:1 , 2:1 , 3:1 ,4:1},{0:3 , 1:1 , 2:1 , 3:1 },{0:1 , 1:1 , 2:1 , 3:3 } ]
  #tdadictlist = [{0:1 , 1:1 } , {0:2} ]#, {0:1 , 1 :2} , {0:2 , 1:1} ]
  tdadictlist = [{0:1}]
  outputreader = dr.ReaderOutput(filen, inputline = '#')
  outputreader.readrgvars(afhvar = None, linenr = linen, startrg = 3)
  rgeq = outputreader.make_rgeq()
  tda = rg.RichardsonSolver(rgeq).tda
  tda.bisect_tda()
  outfile.write(str(rgeq))
  outfile.write('#%s  %s\n' %(outputreader.depvar['depvar'] ,'  '.join(map(str,tdadictlist))))
  while outputreader.depvar['depval'] >= -1.:
    overlaplist = []
    for tdad in tdadictlist:
      tdastate = State_Calculator(tda , rgeq.apair , tdadict =tdad)
      rgeqstate = State_Calculator(rgeq, rgeq.apair)
      ovlap = rgeqstate.calc_overlap(tdastate)
      overlaplist.append(ovlap)
    outfile.write('%f  %s\n' %(outputreader.depvar['depval'], '  '.join(map(str,overlaplist))) )
    linen += linestep
    try:
      outputreader.readrgvars(afhvar = None, linenr = linen, startrg = 3)
    except IndexError:
      break #indexError indicates that we are at the end of the outputfile filen
    rgeq = outputreader.make_rgeq()
    tda.setvar(outputreader.depvar['depvar'], outputreader.depvar['depval'])
    tda.bisect_tda()

  outfile.close()
  plotoverlaps(tdadictlist ,outputreader.depvar['depvar'], filen = outputf)

if __name__ == "__main__":
  #main()
  tdadictlist = [{0:1 , 1:1 , 2:1 , 3:1 ,4:1 , 5:1}, {0:6} , {0:1 , 1:5} , {0:5 , 1:1},{0:1 , 1:1 , 2:4} ,{0:1 , 1:1 , 2:1 , 3:1 ,4:2},{0:4 , 1:1 , 2:1}, {0:2 , 1:1 , 2:1 , 3: 2}, {0:2 , 1:1 , 2:1 , 3:1 ,4:1},{0:3 , 1:1 , 2:1 , 3:1 },{0:1 , 1:1 , 2:1 , 3:3 }]
  plotoverlaps(tdadictlist , 'g' , filen = 'overlapfacintpicketfence.dat', legend = False, xlim = (-4.5 , 1.5))

