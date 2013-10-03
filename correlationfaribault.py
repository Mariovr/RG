# -*- coding: utf-8 -*-
#! /usr/bin/env python
import numpy as np
from numpy import array, zeros, ones, arange,copy
from numpy import append, delete
from math import sqrt
import pylab as pl

import richardsongaudin as rg
import rgfunctions as rgf

class CorrelationFunction(object):
  """
  Implementation of correlation functions as described in:
  Exact mesoscopic correlation functions of the Richardson pairing model, Faribault, Calabrese, Caux, PRB 77 (2008)
  """
  def __init__(self, rgeq):
    self.rgeq = rgeq
    self.gaudinmatrix = self.get_gaudinmatrix()

  def gaudinfunction(self,a,b):
    if a != b:
      return 2./(self.rgeq.rgsolutions[a] -self.rgeq.rgsolutions[b])**2
    else:
      esom = np.sum(1./(self.rgeq.rgsolutions[a] *self.rgeq.energiel)**2) - 2. *np.sum(1./ (self.rgeq.rgsolutions[a] - np.delete(self.rgeq.rgsolutions,a))**2 ) 
      return esom

  def get_gaudinmatrix(self): 
    gmat = zeros((self.rgeq.apair, self.rgeq.apair), np.complex)
    for i in range(self.rgeq.apair):
      for j in range(self.rgeq.apair):
        gmat[i,j] = self.gaudinfunction(i,j)
    return  gmat 

  def K(self, l,q , alpha):
    return (self.rgeq.rgsolutions[l] - self.rgeq.energiel[alpha])/(self.rgeq.rgsolutions[l] - self.rgeq.rgsolutions[q] )

  def get_dmatrix(self,alpha,beta,q):
    dmat = np.zeros((self.rgeq.apair , self.rgeq.apair), np.complex)
    for i in range(0,self.rgeq.apair):
      print dmat
      if i > q:
        dmat[:,i] = self.gaudinmatrix[:,i]
      elif i ==q:
        dmat[:,i] = 1./(self.rgeq.rgsolutions - self.rgeq.energiel[alpha])**2 
      elif i < q-1:
        dmat[:,i] = self.gaudinmatrix[:,i] - self.K(i,q ,alpha)/ self.K(i+1 , q, alpha) * self.gaudinmatrix[:,i+1]
      elif i == q-1:
        dmat[:,i] = self.gaudinmatrix[:,i] + 2. * (self.rgeq.rgsolutions[q] - self.rgeq.energiel[beta])* (self.rgeq.rgsolutions[q-1] - self.rgeq.energiel[alpha])/ (self.rgeq.rgsolutions[q-1] - self.rgeq.rgsolutions[q])*(2.*self.rgeq.rgsolutions - self.rgeq.energiel[alpha] -self.rgeq.energiel[beta])/((self.rgeq.rgsolutions - self.rgeq.energiel[alpha])**2 * (self.rgeq.rgsolutions - self.rgeq.energiel[beta])**2 )
    return dmat

  def get_correlationmp(self , alpha, beta):
    som = 0
    for q in range(self.rgeq.apair):
      dmat = self.get_dmatrix(alpha,beta,q)
      som += (self.rgeq.rgsolutions[q] - self.rgeq.energiel[alpha])/(self.rgeq.rgsolutions[q] - self.rgeq.energiel[beta]) * np.linalg.det(dmat)
    return som 

  def get_correlationzz(self , alpha, beta):
    som = 0
    norm = self.normstate()
    eersteterm = norm/4. 
    for q in range(self.rgeq.apair):
      dmatab = self.get_dmatrix(alpha,beta,q)
      dmatba = self.get_dmatrix(beta,alpha,q)
      som +=  np.linalg.det(dmatab) + np.linalg.det(dmatba)
    som = -1./2.*som
    som += eersteterm
    return som 

  def normstate(self):
    return np.linalg.det(self.gaudinmatrix)

def maintest():
  '''
  test main for main_rgsolver (which is an extraordinary important function) making that function
  more efficient is a huge timewinst. And making it more flexible will cause much better code
  '''
  #picket fence model
  nlev = 6
  eendlev = np.arange(1,nlev+1)
  ontaardingen = np.ones(nlev,float)*2
  senioriteit = np.zeros(nlev,float)
  apair = nlev/2
  #artikel Stefan test
  '''
  eendlev = array([0.04,0.08,0.16,0.20,0.32,0.36,0.40,0.52,0.64,0.68,0.72,0.8,1])
  eendlev = np.sqrt(eendlev)
  senioriteit = zeros(len(eendlev),float)
  ontaardingen = [4,4,4,8,4,4,8,8,4,8,4,8,12]
  apair = 10
  g = -0.075
  '''
  eta = 1.
  g = -0.0100
  tdastartd = {0:apair }
  tdastartd = rgf.tdadict_kleinekoppeling(apair,ontaardingen,senioriteit)
  alevel = len(eendlev)
  #print ontaardingen,senioriteit,eendlev
  alpha = 0
  beta =0
  assert(len(ontaardingen) == len(eendlev) and len(ontaardingen) == len(senioriteit))
  rgeq = rg.RichRedBcs(eendlev,ontaardingen,senioriteit,g,apair)
  a = rg.RichardsonSolver(rgeq)
  rgeq = a.main_solve(tdastartd,xistep = 0.01,xival = 1.,rgwrite = True,plotrgvarpath = True , plotepath = True,xlim = None , ylim = None)
  corcal = CorrelationFunction(rgeq)
  mp = [] ; zz = [] ; xw = []
  for i in range(nlev):
    cmpp = corcal.get_correlationmp(alpha,i)
    czz = corcal.get_correlationzz(alpha,i)
    print 'The particle particle correlation is :' , cmpp
    print 'The z z correlation is: ', czz
    mp.append(cmpp.real)
    zz.append(czz.real)
    xw.append(i/float(nlev))
  pl.plot(xw , mp )
  pl.title('mp')
  pl.plot(xw , zz )
  pl.title('zz')
  pl.show()

if __name__ == "__main__":
  maintest()


