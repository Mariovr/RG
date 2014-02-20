# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Mario Van Raemdonck, 2013;
# (c) Ghent University, 2013
# -*- coding: utf-8 -*-
#! /usr/bin/env python
import numpy as np
from numpy import array, zeros, ones, arange,copy , fromiter , reshape
from numpy import append, delete
from math import sqrt
import pylab as pl

import richardsongaudin as rg
import rgfunctions as rgf
import writepairing as wp
import timergfunctions as tf
import overlaps as ov
import datareader as dr

class CorrelationFunction(object):
  """
  Implementation of correlation functions as described in:
  Exact mesoscopic correlation functions of the Richardson pairing model, Faribault, Calabrese, Caux, PRB 77 (2008)
  REMARK REMARK REMARK: The energylevels faribault et. al. use are defined as the pair energies they are twice the sp energies used in my solver
  so I need to multiply rgeq.energiel with 2 to replace the energy in their formulas to obtain the correct answers for the correlationcoefficients.
  """
  def __init__(self, rgeq , calc_2rdm = True):
    if rgeq != None:
      self.rgeq = rgeq 
      self.re_calc(calc_2rdm)
    else:
      pass
    
  def set_rgeq(self,rgeq ,calc_2rdm = True):
    self.rgeq = rgeq
    self.re_calc(calc_2rdm)

  def re_calc(self, calc_2rdm = True ):
    #self.rgeq.solve() #check that the Richardson-Gaudin variables are already defined and if something else is set recalc the RG vars, not necessary because we command that a CoorrelationFunction calculator has a rgeq with a proper solution.
    self.gaudinmatrix = self.get_matrix(self.rgeq.apair, self.gaudinfunction)
    self.norm = self.get_norm()
    self.rdm1 = self.get_1rdm()
    if calc_2rdm:
      self.calc_2rdm()

  def calc_2rdm(self):
    self.nn = self.get_nn()
    self.pm = self.get_pm()

  def get_element_2rdm(self,i,j,k,l):
    if i== j and k == l:
      return self.pm[i,k]
    elif i == k and j == l : 
      return self.nn[i,j]
    else:
      return 0.

  def get_matrix(self, dim , function ,*args):
    gmat = fromiter((function(i,j, *args) for i in xrange(dim) for j in xrange(dim)), np.complex)         
    gmat = reshape(gmat , (dim,dim))
    return gmat
  
  def gaudinfunction(self,a,b):
    if a != b:
      return 2./(self.rgeq.rgsolutions[a] -self.rgeq.rgsolutions[b])**2.
    else:
      return np.sum(1./(self.rgeq.rgsolutions[a] -self.rgeq.energiel*2.)**2.) - 2. *np.sum(1./ (self.rgeq.rgsolutions[a] - np.delete(self.rgeq.rgsolutions,a))**2. ) 

  #2 functions to calculate 1rdm (which is diagonal)
  def qfunc(self,b,a,i):
    if a == b:
      return np.prod(np.delete(self.rgeq.rgsolutions,a)- self.rgeq.rgsolutions[a] ) * (1./2. *sum(1./(self.rgeq.energiel*2 - self.rgeq.rgsolutions[a])**2)- sum(1./(np.delete(self.rgeq.rgsolutions,a) - self.rgeq.rgsolutions[a])**2) - 1./(self.rgeq.energiel[i]*2 - self.rgeq.rgsolutions[a] )**2)
    else:
      return np.prod(np.delete(self.rgeq.rgsolutions,a)- self.rgeq.rgsolutions[a] ) * (1./(self.rgeq.rgsolutions[b] - self.rgeq.rgsolutions[a])**2 - 1./(self.rgeq.energiel[i]*2 - self.rgeq.rgsolutions[b] )**2)
      
  def get_1rdm(self): #This is tested and working!
    #REMARK is diagonal in the seniority zero space so we return a one dimensional list
    prod = 1. ; rdm = []
    for b in range(self.rgeq.apair):
      prod *= np.prod(self.rgeq.rgsolutions[b] - np.delete(self.rgeq.rgsolutions,b))
    for i in range(self.rgeq.alevel):
      qmat = self.get_matrix(self.rgeq.apair,self.qfunc, i)
      rdm.append(1./2.-2**(self.rgeq.apair-1)/self.get_norm() * np.linalg.det(qmat)/prod )
    return rdm

  def get_sz(self,i): #This S_z is working
    return self.rdm1[i] - self.rgeq.ontaardingen[i]/4.
    
  #2 functions to calculate the overlap with other pairstates (not yet working)
  def jfunction(self,a,b, rgsol2):
    return (rgsol2[b] - self.rgeq.rgsolutions[b]) / (rgsol2[a] - self.rgeq.rgsolutions[b])*(np.sum(1./((rgsol2[a] - self.rgeq.energiel*2)*(self.rgeq.rgsolutions[b]- self.rgeq.energiel*2))) - 2 *np.sum(1./((rgsol2[a]-np.delete(rgsol2,a))*(self.rgeq.rgsolutions[b] -  np.delete(rgsol2,a)))))

  def get_overlap(self, rgsol2):
    jmat = self.get_matrix(self.rgeq.apair , self.jfunction, rgsol2.rgeq.rgsolutions)
    prod = 1. ; prod2 = 1.; prod3 = 1.
    for b in range(self.rgeq.apair):
      prod *= np.prod(rgsol2.rgeq.rgsolutions[b] - np.delete(self.rgeq.rgsolutions,b))
      if b < self.rgeq.apair-1:
        prod2 *= np.prod(self.rgeq.rgsolutions[b] - self.rgeq.rgsolutions[b+1:])
        prod3 *= np.prod(rgsol2.rgeq.rgsolutions[b+1:] - rgsol2.rgeq.rgsolutions[b])
    return prod/(prod2*prod3)* np.linalg.det(jmat) /sqrt(self.get_norm() *rgsol2.get_norm())

  #functions to calculate the 2rdm
  def K(self, l,q , alpha):
    return (self.rgeq.rgsolutions[l] - 2.*self.rgeq.energiel[alpha])/(self.rgeq.rgsolutions[l] - self.rgeq.rgsolutions[q] )

  def get_dmatrix(self,alpha,beta,q): #tested and working :)
    assert(alpha != beta)
    dmat = np.zeros((self.rgeq.apair , self.rgeq.apair), np.complex)
    for i in range(0,self.rgeq.apair):
      if i > q:
        dmat[:,i] = self.gaudinmatrix[:,i]
      elif i ==q:
        dmat[:,i] = 1./(self.rgeq.rgsolutions - 2.*self.rgeq.energiel[alpha])**2 
      elif i < q-1:
        dmat[:,i] = self.gaudinmatrix[:,i] - self.K(i,q ,alpha)/ self.K(i+1 , q, alpha) * self.gaudinmatrix[:,i+1]
      elif i == q-1:
        dmat[:,i] = self.gaudinmatrix[:,i] + 2.*(2.*self.rgeq.energiel[beta]- self.rgeq.rgsolutions[q] )*(2.*self.rgeq.energiel[alpha] - self.rgeq.rgsolutions[q-1])/ (self.rgeq.rgsolutions[q-1] - self.rgeq.rgsolutions[q]) *(2.*self.rgeq.rgsolutions - 2.*self.rgeq.energiel[alpha] -2.*self.rgeq.energiel[beta])/ ((self.rgeq.rgsolutions - 2.*self.rgeq.energiel[alpha])**2 * (self.rgeq.rgsolutions - 2.*self.rgeq.energiel[beta])**2 )
    return dmat

  def rdm2_func(self , alpha, beta):
    #<S+_i S_j>
    som = 0
    if alpha != beta:
      for q in range(self.rgeq.apair):
        dmat = self.get_dmatrix(alpha,beta,q)
        som += (self.rgeq.rgsolutions[q] - 2.*self.rgeq.energiel[alpha])/(self.rgeq.rgsolutions[q] - 2.*self.rgeq.energiel[beta]) * np.linalg.det(dmat)
      som /= self.norm
    else:
      som = self.rdm1[alpha]
    return som 
    
  def nn_func(self , alpha, beta): 
    #<n_in_j>
    som = 0
    if alpha != beta:
      for q in range(self.rgeq.apair):
        dmatab = self.get_dmatrix(alpha,beta,q) 
        dmatba = self.get_dmatrix(beta,alpha,q) 
        som +=  np.linalg.det(dmatab) + np.linalg.det(dmatba)
      som = -1./2.*som/ self.norm
      som += 1/2. * (self.rdm1[alpha] + self.rdm1[beta])
    else:
      som = self.rdm1[alpha]
    return som  

  def get_ss(self): #<S_i S_j> = <n_i n_j > -1/2 (n_i + n_j ) + ||w||^2/4 , remark the n_i = a^+_i a_i and not a^+_i a_i+a^+_\bar{i} a_\bar{i}
    #REMARK trivially : <(S^z_alpha)^2> = 1./4.
    nn = self.get_nn()
    nex = 1./2. * array([[i + j for i in self.rdm1] for j in self.rdm1])
    nn -= nex 
    nn += 1./4.
    return nn

  def get_nn(self): # tested and working
    return self.get_matrix(self.rgeq.alevel, self.nn_func)

  def get_pm(self):
    return self.get_matrix(self.rgeq.alevel, self.rdm2_func)

  #the norm
  def get_norm(self):
    return np.linalg.det(self.gaudinmatrix)

def maintest():
  '''
  test main for the correlationcoef calculator : reproduces the plots on page 11 of :
  exact mesoscopic correlation functions of the Richardson pairing model. Faribault et. al. PRB 2008
  '''
  #picket fence model
  nlev = 20
  eendlev = [i/2. for i in range(1,nlev+1)] #remark to have the same energy levels as Faribault et. al. we should take the half because they use another convention
  ontaardingen = np.ones(nlev,float)*2
  senioriteit = np.zeros(nlev,float)
  apair = nlev/2
  #tdastard = {0:apair }
  tdastartd = rgf.tdadict_kleinekoppeling(apair,ontaardingen,senioriteit)
  alevel = len(eendlev)
  alpha = 0
  beta =0
  assert(len(ontaardingen) == len(eendlev) and len(ontaardingen) == len(senioriteit))
  rgeq = rg.RichRedBcs(eendlev,ontaardingen,senioriteit,-0.001,apair)
  rgf.generate_dir("cortest",None)
  for g in [-0.1 , -0.2 , -0.4, -0.5 , -0.6, -1.]: #The interaction constants choosen by Faribault et. al. in their PRB, 2008 publication
    energierg , rgeq =  wp.generating_datak(rgeq,tdastartd,'g',-0.001,g,rgwrite = True ,tdafilebool = False,exname = '',moviede = False, intofmotion = False, printstep = 30 )       
    corcal = CorrelationFunction(rgeq)
    norm = corcal.norm 
    print corcal.norm
    #tf.timer(corcal.get_gaudinmatrix2 , repetitions = 100)
    #tf.timer(corcal.get_gaudinmatrix, repetitions = 100)
    print '1rdm is :'
    print corcal.get_1rdm()
    mp = [] ; zz = [] ; xw = []
    cmpp = corcal.get_pm()
    czz = corcal.get_ss() 
    print 'The particle particle correlation is :' , cmpp
    print 'The z z correlation is: ', czz
    mp = cmpp[alpha,alpha+1:]
    zz = czz[alpha,alpha+1:]
    xw = [i/float(nlev) for i in range(2,nlev+1)]
    pl.plot(xw , mp )
    pl.title('mp')
    pl.savefig('mp%s.png' %(str(g)))
    pl.close()
    pl.plot(xw , zz )
    pl.title('zz')
    pl.savefig('zz%s.png' %str(g))
    pl.close()

def testoverlap():
  #test of the 1rdm and 2rdm elements with the little program of Paul Johnson (PASSED), and a test for the overlap calculation with a brute force overlap calculation (OVERLAPS ARE NOT YET EQUAL)
  step = -0.001 ; afhv = -0.238 ; end = -0.238
  file1 = 'plotenergy.dat'
  file2 = 'plotenergy2.dat'
  outputf = 'overlap_1rdm.dat' 
  outfile = open('overlap_1rdm.dat' ,'w')
  outfile2 =open('2rdm.dat' , 'w')
  outputreader = dr.ReaderOutput(file1, inputline = '#')
  outputreader.readrgvars(afhvar = afhv,linenr = None, startrg = 3)
  rgeq = outputreader.make_rgeq()
  outputreader2 = dr.ReaderOutput(file2, inputline = '#')
  outputreader2.readrgvars(afhvar = afhv,linenr = None, startrg = 3)
  #outputreader.tdadist contains the start tda distribution of the state
  rgeq2 = outputreader2.make_rgeq()
  outfile.write(str(rgeq))
  outfile.write('#%s  corr_far  bruteforce\n' %(outputreader.depvar['depvar'] ))
  corfunc = CorrelationFunction(rgeq)
  corfunc2 = CorrelationFunction(rgeq2)
  while outputreader.depvar['depval'] >= end:
    #overlap with brute force from the overlaps.py file
    rgeqstate = ov.State_Calculator(corfunc.rgeq, corfunc.rgeq.apair, normalize = False)
    ovlap = rgeqstate.calc_overlap(ov.State_Calculator(corfunc2.rgeq,corfunc2.rgeq.apair, normalize = False))
    #overlap with fast method
    overlap = corfunc.get_overlap(corfunc2)
    rdm =  corfunc.get_1rdm()
    rdmnn = corfunc.get_nn()
    rdmpm = corfunc.get_pm()
    print corfunc.get_dmatrix(7,5,1)
    outfile.write('%f  %f  %f  %f  %f  %s\n' %(outputreader.depvar['depval'], overlap ,ovlap, corfunc.get_norm(), rgeqstate.get_norm(), rdm) )
    outfile2.write('%f  %s  %s\n' %(outputreader.depvar['depval'], str(rdmnn) ,str(rdmpm) ))
    afhv += step
    try:
      outputreader.readrgvars(afhvar = afhv, linenr = None, startrg = 3)
      outputreader2.readrgvars(afhvar = afhv, linenr = None, startrg = 3)
    except IndexError:
      break #indexError indicates that we are at the end of the outputfile filen
    corfunc.set_rgeq(outputreader.make_rgeq())
    corfunc2.set_rgeq(outputreader2.make_rgeq())

  outfile.close()
  outfile2.close()

if __name__ == "__main__":
  maintest()
  #testoverlap()
