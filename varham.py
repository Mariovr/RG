# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Mario Van Raemdonck, 2014;
# (c) Ghent University, 2014
#
# We make use of the horton package for theoretical chemistry developed by Toon Verstraelen et. al. (http://theochem.github.io/horton/tut_getting_started.html)
#
# -*- coding: utf-8 -*-
#! /usr/bin/env python 
from horton import *

import correlationfaribault as cf
import richardsongaudin as rg
import rgfunctions as rgf
import writepairing as wp
import datareader as dr
import varRG as vrg

import pylab as pl
import sys
import os

"""
This file contains some implementations of Hamiltonians that can be approximated by our variational RG method: for example the nonrelativistic quantum chemical Hamiltonian, and the pairing Hamiltonian with a level dependend interaction constant.
"""

class Hamiltoniaan(object):
  """
  Abstract base class of Hamiltonians we are going to investigate variationally with our Richardson-Gaudin ansatz.
  """
  def calc_energy(self,corfunc):
    raise Exception("Implement this function in a subclass from Reader")

  def get_2p_integrals(self):
    raise Exception("Implement this function in a subclass from Reader")
    
  def get_1p_integrals(self): 
    raise Exception("Implement this function in a subclass from Reader")

  def get_rgeq_char(self):
    raise Exception("Implement this function in a subclass from Reader")

class Chem_Ham(Hamiltoniaan):
  """
  The non-relativistic quantum chemical Hamiltonian.
  """
  def __init__(self , filename = None , atoms = [4] , positions =[[0,0,0]], basis = '3-21G'):
    #if filename is not None we read the molecular system data from a .xyz file
    self.horton_ham = self.meanfield(filename , atoms , positions , basis = basis)
    self.eint , self.eri   = self.transforms()
    self.nucnuc = self.horton_ham.system.compute_nucnuc()
    self.nbasis = self.horton_ham.system.get_kinetic().nbasis

  def calc_energy(self , rgcor):
    som = 0
    for i in range(self.nbasis):
      som += 2.*(self.eint[i,i]) * rgcor.rdm1[i] #1rdm is diagonal
    for i in range(self.nbasis):
      for j in range(self.nbasis):
        if i != j:
          som += (2.*self.eri[i,j,i,j] - self.eri[i,j,j,i] )*rgcor.get_element_2rdm(i,j,i,j)
          som += self.eri[i,i,j,j] *rgcor.get_element_2rdm(i,i,j,j)
      som += self.eri[i,i,i,i] *rgcor.get_element_2rdm(i,i,i,i)
    som += self.nucnuc
    print 'The energy = %f' %som
    return som
  
  def get_2p_integrals(self):
    return self.eri
    
  def get_1p_integrals(self): 
    return self.eint

  def get_rgeq_char(self):
    for hfspbasis in self.horton_ham.system.wfn._iter_expansions(): #For restricted calculations this runs only over the alpha electrons because the beta are equal
      elev = list(hfspbasis.energies)
      ap = int(sum(hfspbasis.occupations))
    deg = np.array([2.] *len(elev))
    sen = np.array([0.] *len(elev))
    return elev , deg , sen , ap

  def meanfield(self,filen ,atoms, positions ,  basis ):
    #possible bases:  '3-21G' , STO-3G
    if filen != None:
      # Load the coordinates from file (.xyz file)
      sys = System.from_file(filen, obasis= basis)
    else:
      #Or define the system by hand
      sys = System(positions, atoms ,  obasis= basis)  #standard one berryllium atom
    #print sys.get_kinetic()._array ; #print sys.get_overlap()._array ; #print sys.get_nuclear_attraction()._array ; #print sys.get_electron_repulsion()._array
    # Initialize the closed-shell wfn
    setup_mean_field_wfn(sys, charge=0)
    # Initial WFN guess
    guess_hamiltonian_core(sys)
    # Construct a Hamiltonian
    ham = Hamiltonian(sys, [HartreeFockExchange()])
    # Converge WFN with SCF
    converged = converge_scf(ham)
    # Compute the energy
    log.set_level(log.high)
    ham.compute()
    log.set_level(log.medium)
    #ham.system.to_file('ber.h5')
    return ham

  def transforms(self):
    eint  = self.horton_ham.system.get_kinetic().dot(self.horton_ham.system.wfn.get_exp('alpha').coeffs.T ,self.horton_ham.system.wfn.get_exp('alpha').coeffs )
    eint += self.horton_ham.system.get_nuclear_attraction().dot(self.horton_ham.system.wfn.get_exp('alpha').coeffs.T ,self.horton_ham.system.wfn.get_exp('alpha').coeffs )
    eri = np.tensordot(self.horton_ham.system.get_electron_repulsion()._array , self.horton_ham.system.wfn.get_exp('alpha').coeffs ,axes =  ([0],[0]))
    eri = np.tensordot(eri , self.horton_ham.system.wfn.get_exp('alpha').coeffs ,axes =  ([0],[0]))
    eri = np.tensordot(eri , self.horton_ham.system.wfn.get_exp('alpha').coeffs ,axes =  ([0],[0]))
    eri = np.tensordot(eri , self.horton_ham.system.wfn.get_exp('alpha').coeffs ,axes =  ([0],[0]))
    return eint , eri

  def get_overlap():
    return self.horton_ham.system.get_overlap()._array

  def create_densematrices(self ,corrgeq):
    #Create dense matrices (to work with horton) from 1rdm and 2rdm matrices from a Richardson-Gaudin ground state
    dtwb = DenseTwoBody(corrgeq.rgeq.alevel)
    dob = DenseOneBody(corrgeq.rgeq.alevel)

    for i in range(corrgeq.rgeq.alevel): 
      dob.set_element(i,i,rgeqdm.rdm1[i])
    for i in range(corrgeq.rgeq.alevel): 
      for j in range(i,corrgeq.rgeq.alevel): 
        for k in range(j,corrgeq.rgeq.alevel): 
          for l in range(k,corrgeq.rgeq.alevel): 
            dtwb.set_element(i,j, k , l,rgeqdm.get_element_2rdm(i,j,k,l))
    dtwb.check_symmetry()
    dob.check_symmetry()
    return dob, dtwb

class General_Pairing_Ham(Hamiltonian):
  """
  The pairing Hamiltonian with a level dependend interaction constant.
  """
  def __init__(self, eps ,deg , sen , gij , pair):
    """
    gij determines how level i couples to level j (2d array)
    eps are the sp energies (1d array)
    """
    self.gij = gij
    self.apair = pair
    self.gij = self.set_twofold_deg( gij , deg)
    self.eps, self.deg,self.sen, self.alevel = rgf.uncheckdegeneratie(eps, deg)

  def calc_energy(self, rgcor):
    som = 0
    for i in range(self.alevel):
      som += 2.*(self.eps[i]) * rgcor.rdm1[i] #1rdm is diagonal
    for i in range(self.alevel):
      for j in range(self.alevel):
        if i != j:
          som += self.gij[i,j] *rgcor.get_element_2rdm(i,i,j,j)
      som += self.gij[i,i] *rgcor.get_element_2rdm(i,i,i,i)
    print 'The energy = %f' %som
    return som

  def get_2p_integrals(self):
    return self.gij
    
  def get_1p_integrals(self): 
    return self.eps

  def get_rgeq_char(self):
    return self.eps , self.deg , self.sen , self.apair

  def set_twofold_deg(self, gij , deg):
    totaln = sum(deg)/2.
    newgij = np.zeros((totaln , totaln))
    isom = 0 
    for i in range(len(deg)):
      jsom = 0
      for j in range(len(deg)):
        for d in range(int(deg[j]/2)):
          newgij[isom,jsom +d] = gij[i][j]
        jsom += d+1
      for k in range(int(deg[i]/2)):
        newgij[isom+k] = newgij[isom]
      isom += k+1 
    return newgij

def nucmain():
  name = 'Sn_116'
  nuc2 = [[-0.24625, -0.16486, -0.14600, -0.18328, -0.23380],
  [-0.16486, -0.23543, -0.19953, -0.36971, -0.22500],
  [-0.14600, -0.19953, -0.72440, -0.17412, -0.24855],
  [-0.18328, -0.36971, -0.17412, -0.17665, -0.17615],
  [-0.23380, -0.22500, -0.24855, -0.17615, -0.20315]] #the pairing Ham changes this matrix automatically to twofold degenerate levels
  enl=  [-6.1210, -5.5080, -3.8910,  -3.7780 , -3.749] ; deg= [8, 6, 2, 12, 4] ; sen= [0, 0, 0, 0, 0 ] ; ap = 8
  vrg.GoldVarRG(None).create_output_file(name,enl , deg, sen,ap, -0.001 , -0.0001 , -0.4)
  outputreader = dr.ReaderOutput('plotenergy.dat')
  nucham = General_Pairing_Ham(outputreader.elevels , outputreader.degeneracies, outputreader.seniorities , nuc2 , outputreader.npair)
  outputreader.elevels , outputreader.degeneracies , outputreader.seniorities, outputreader.npair  = nucham.get_rgeq_char()
  elist = [] ; x = []
  end = -0.1 ; linestep = -10 ; linen = -1
  with open('%s.dat' %name , 'w') as outputfile:
    outputreader.readrgvars(afhvar = None ,linenr = linen, startrg = 3) #reading of solutions: remark xifile has startrg = 2, plotgeofile has startrg = 6 , plotenergyfile has startrg = 3 
    while outputreader.depvar['depval'] <= end:
      outputreader.readrgvars(afhvar = None ,linenr = linen, startrg = 3) #reading of solutions: remark xifile has startrg = 2, plotgeofile has startrg = 6 , plotenergyfile has startrg = 3 
      rgeq = outputreader.make_rgeq() 
      linen += linestep
      rgeqdm = cf.CorrelationFunction(rgeq, calc_2rdm = True)
      elist.append(nucham.calc_energy(rgeqdm) )
      x.append(outputreader.depvar['depval'])
      outputfile.write('%f\t%f\n' % (outputreader.depvar['depval'] , elist[-1]))
  pl.plot(x, elist)
  pl.savefig('%s.png' % name)
  

def chemmain(*args , **kwargs):
  #possible bases: ['STO-3G' , '3-21G' , '3-21++G*' , '6-31++G**', '6-31G**','6-31+G*' , 'ANO' ,'aug-cc-pVDZ', 'aug-cc-pVTZ', 'aug-cc-pVQZ', 'cc-pVDZ', 'cc-pVTZ', 'cc-pVQZ'  ]
  #to add other bases download from 
  name = 'Be'
  chemham = Chem_Ham(filename = None ,atoms = [4], positions = [[0,0,0]],  basis =  '3-21G') #3-21G
  ehf = chemham.horton_ham.system.extra['energy']
  #efci = -127.9189885238 #get this from gaussian fci ne3-21G
  efci =-14.5314443939 #Be fci 3-21G
  #efcibeaug = -14.6157221969 
  elev, deg , sen, ap  = chemham.get_rgeq_char()
  g = 0.0001 ; elist = []
  nlevel,deg = rgf.checkdegeneratie(elev ,list(deg))
  sen = nlevel* [0.]
  rgeq = vrg.VarSolver(None).create_rgeq(elev,deg, sen , g , ap, tda = 'weak')
  rgeq.energiel , rgeq.ontaardingen ,rgeq.senioriteit ,rgeq.alevel= rgf.uncheckdegeneratie(rgeq.energiel , rgeq.ontaardingen)
  with open('%s.dat' %name , 'w') as d:
    for g in [ 0.001 * i for i in range(1,201)]:
      rgeqdm = cf.CorrelationFunction(rgeq, calc_2rdm = True)
      elist.append(chemham.calc_energy(rgeqdm))
      corp = (ehf -elist[-1])/(  ehf - efci)
      d.write('%f\t%f\t%f\n' % (g , elist[-1], corp))
      energierg , rgeq =  wp.generating_datak(rgeqdm.rgeq, None,'g',0.001,g,rgwrite = True ,tdafilebool = False,exname = '',moviede = False, intofmotion = False, printstep = 30 )       
  pl.plot([ 0.001 * i for i in range(1,201)], elist)
  pl.savefig('%s.png' % name)

if __name__ == "__main__":
  #chemmain()
  nucmain()
