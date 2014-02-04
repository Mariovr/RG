# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Mario Van Raemdonck, 2013;
# (c) Ghent University, 2013
#! /bin/env python 
# -*- coding: utf-8 -*-

import numpy as np
import pylab as pl

import newtonraphson as nr
import richardsongaudin as rg

class Alpha:
  '''
  class that defines the system of equations, that need to be solved to determine
  the correct alpha values for the calculation of the correlation coefficients
  '''
  def __init__(self,eenlevi,rgvari):
    '''
    initialisation of the alpha class
    '''
    self.eenlev = eenlevi
    self.rgvar = rgvari
    
  def __call__(self,alpha):
    '''
    creates the system to define the alpha variables
    and makes the class callabel
    '''
    npair = len(self.rgvar)
    nlevel = len(self.eenlev)
    stelsel = np.zeros(npair,complex)
    for u in xrange(npair):
      cu = 0
      for i in xrange(nlevel):
        cu += (1./(2.*self.eenlev[i]-self.rgvar[u]))**2.
      su = 0
      term2 = 0
      for eta in xrange(npair):
        if eta != u:
          su += (1./(self.rgvar[eta]-self.rgvar[u]))**2
          term2 +=  alpha[eta]*(1./(self.rgvar[eta]-self.rgvar[u]))**2
      stelsel[u] = (cu - 2*su)*alpha[u]+ term2 -1.
    return stelsel
  
  def plotVgl(self,alpharange,nalpha,nvgl):
    '''
    Plots one equation of the stelsel to determine alpha 
    in function of one alpha
    '''
    x = np.arange(-100,100,0.1)
    y = np.zeros(len(x),float)
    for i in xrange(len(x)):
      alpharange[nalpha] = x[i]
      y[i] = self(alpharange)[nvgl].real
    pl.figure()
    pl.plot(x,y)
    pl.savefig('Alphatestplot')
  
  def solveAlpha(self):
    '''
    Function that solves for the alpha variables that are needed to determine the the correlation coefficients
    '''
    print '-------------------------------------------'
    print 'start solving for alpha'
    print '-------------------------------------------'
    alphaguess = self.generateAlphaguess(self.eenlev,self.rgvar)
    alpha = nr.solve(self,alphaguess)
    return alpha
    
  def generateAlphaguess(self,eenlev,rgvar):
    '''
    Function that generates the first guess for the apair alpha variables.
    REMARK: we need to think over a good guess because now it's very rudimentary
    '''
    alphag = np.ones(len(rgvar))  *10
    return alphag


  
def cijApprox(eenlev,rgvar,alpha , i , j):
  '''
  Calculates the correlation coefficient cij by an approximation method
  for the details see the appendix of the review: spectroscopy of discrete energylevels in ultrasmall metalic grains
  '''
  som = 0.
  som2 = 0.
  for b in xrange(len(rgvar)):
    som += 1./(2.*eenlev[j]-rgvar[b])
  for a in xrange(len(rgvar)):
    som2 += alpha[a] * 1./(2.*eenlev[i]-rgvar[a])
  cij = som2 * som
  return cij

def corApproxFault(g, eenlev, alpha, rgvar):
  '''
  Calculates the fault of the used approximation method of the correlation coefficients
  if the variable it returns is big -> bad approximation
  '''
  som = 0
  for i in xrange(len(rgvar)):
    somel = 0
    for j in xrange(len(eenlev)):
      somel += 1./(2.*eenlev[j]-rgvar[i])
    somel = somel**2.
    som += (somel *g**2 -1)*alpha[i]
  return som/float(g)

def writeCor(cij, fault,afh, filename = 'correlaties'):
  '''
  write correlationcoefficients to a file
  '''
  cfile = open(filename,'a')
  filename.write('%f\t%f\t%f\n' %(afh,cij,fault))
  cfile.close()

def userinput():
  """
  function that asks the user for some input
  """
  aantalel = input("Give the number of sp levels")
  #g is negatief voor aantrekkende interactie en positief voor afstotende interactie
  g = input("Give the interaction constant")
  apair = input('Give the number of pairs')
  eendlev = np.zeros(aantalel,np.float)
  ontaardingen = np.zeros(aantalel,np.float)
  for i in xrange(aantalel):
    a = input("Give the sp state?")
    b = input("Give the degeneration of that sp state")
    eendlev[i] = a
    ontaardingen[i] = b
  return g,eendlev,ontaardingen,apair
  
  
def main():
  """
  Some test's of the implementation of the correlation-coefficients from the Richardson-Gaudin equation's
  """
  
  #g,eendlev,ontaardingen,apair = userinput()
  g = -0.1
  apair = 2
  eendlev = [0.,1.,2.,3.]
  eendlev = np.array(eendlev)
  kk = True
  tdastartd = {}
  tdastartd = {0:apair}
  if kk is True:
    for i in range(apair):
      tdastartd[i] = 1.
  aantalel = len(eendlev)
  senioriteit = np.zeros(aantalel,float)
  ontaardingen = np.ones(aantalel,np.float)*2.
  rgenergie, rgvars = rg.main_rgsolver(eendlev,ontaardingen,g,apair,tdastartd,senioriteit,firstguess = None,rgwrite = False)
  checkstelsel = Alpha(eendlev,rgvars)
  #alphaguess = checkstelsel.generateAlphaguess(eendlev,rgvars)
  #checkstelsel.plotVgl(alphaguess,0,1)
  alphas = checkstelsel.solveAlpha()
  print 'the alpha variables are: '
  print alphas
  cor = np.zeros((aantalel,aantalel),complex)
  faulcor = np.zeros((aantalel,aantalel),complex)
  for i in xrange(aantalel):
    for j in xrange(aantalel):
      cor[i,j] =  cijApprox(eendlev,rgvars,alphas , i , j)
      faulcor[i,j] = corApproxFault(g, eendlev, alphas, rgvars)
  print '---------------the correlation coefficients --------------------------------'
  print cor
  print '-----------------------------------------------------------------------------'
  print '---------------the corresponding faults are ---------------------------------'
  print faulcor
  print '-----------------------------------------------------------------------------'

  
if __name__ == "__main__":
  main()

