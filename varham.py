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
from scipy import linalg

import pylab as pl
import sys
import os

import numpy as np

sys.path.append('../../CIFlow/rundir')
import read_psi as rp

"""
This file contains some implementations of Hamiltonians that can be approximated by our variational RG method: for example the nonrelativistic quantum chemical Hamiltonian, and the pairing Hamiltonian with a level dependend interaction constant.
"""
class Hamiltoniaan(object):
  """
  Abstract base class of Hamiltonians we are going to investigate variationally with our Richardson-Gaudin ansatz.
  """
  def calc_energy(self,corfunc):
    raise Exception("Implement this function in a subclass from Hamiltoniaan")

  def get_2p_integrals(self):
    raise Exception("Implement this function in a subclass from Hamiltoniaan")
    
  def get_1p_integrals(self): 
    raise Exception("Implement this function in a subclass from Hamiltoniaan")

  def get_rgeq_char(self):
    raise Exception("Implement this function in a subclass from Hamiltoniaan")

  def calc_energy_rgeq(self , rgeq):
    #rgeq.energiel , rgeq.ontaardingen ,rgeq.senioriteit ,rgeq.alevel= rgf.uncheckdegeneratie(rgeq.energiel , rgeq.ontaardingen)
    rgcor = cf.CorrelationFunction(None, calc_2rdm = True)
    rgcor.set_rgeq(rgeq, calc_2rdm = True)
    return self.calc_energy(rgcor)

class Psi4_ham(Hamiltoniaan):
  """
  The non-relativistic quantum chemical Hamiltonian. With matrixelements extracted from psi4.
  """
  def __init__(self , filename = None , diagtp = False):
      reader = rp.PsiReader(filename, isbig = False, numorbs = None, read_ints = True)
      self.twoindex , self.fourindex = reader.list_to_matrix()
      if diagtp:
          evals , evecs = self.diagonalize(self.twoindex)
          self.twoindex , self.fourindex = reader.transform_integrals(evecs)
      self.norb = reader.values['norb']
      self.nucrep =  reader.values['nucrep']
      self.nup = reader.values['nalpha']
      self.permarray = [] #when the epsilons of rich red ham are sorted this keeps the permutationarray.

  def calc_energy(self , rgcor):
    som = 0
    #print self.permarray
    for i in range(self.norb):
      som += 2.*(self.twoindex[self.permarray[i],self.permarray[i]]) * rgcor.rdm1[i] #1rdm is diagonal
      #print som ,'  ' ,self.twoindex[self.permarray[i],self.permarray[i]]
    for i in range(self.norb):
      for j in range(self.norb):
        if i != j:
          som += (2.*self.fourindex[self.permarray[i],self.permarray[i],self.permarray[j],self.permarray[j]] - self.fourindex[self.permarray[i],self.permarray[j],self.permarray[j],self.permarray[i]] )*rgcor.get_element_2rdm(i,j,i,j)
          som += self.fourindex[self.permarray[i],self.permarray[j],self.permarray[i],self.permarray[j]] *rgcor.get_element_2rdm(i,i,j,j)
      som += self.fourindex[self.permarray[i],self.permarray[i],self.permarray[i],self.permarray[i]] *rgcor.get_element_2rdm(i,i,i,i)
    som += self.nucrep
    print 'The energy = %f' %som
    return som

  def get_2p_integrals(self):
    return self.fourindex

    
  def get_1p_integrals(self): 
    return self.twoindex

  def diagonalize(self, mat):
      eigval , eigvec = linalg.eigh(mat)
      return (eigval, eigvec)

  def get_rgeq_char(self):
    elev = []
    for i in range(self.norb):
        som = 0
        som += self.twoindex[i,i]
        #for j in range(self.nup): #we run only over the nup lowest occupied orbitals for the HF sp energies, also for the virtual orbitals.
           #som += 1./2.*(self.fourindex[i,i,j,j]*2 - self.fourindex[i,j,j,i]) #the extra factor of 1/2 is to make sure that the energy of the system at g = 0, equals the HF energy so any, lowering of the energy by turning on g will improve on the HF values. (because the sum of the fock energies is not exactly equal to HF energy, up on this factor of 1/2).
           #som += self.fourindex[i,i,j,j]*2 - self.fourindex[i,j,j,i] #This should be equal to the Fock eigenvalues.
        elev.append(som)

    elev = range(1,self.norb+1)
    deg = np.array([2.] * self.norb)
    sen = np.array([0.] * self.norb)
    #self.hfe = 2*sum(elev[:self.nup]) + self.nucrep
    #print 'HF energy: ' , self.hfe
    print elev
    self.permarray = sorted(range(len(elev)), key=lambda k: elev[k])
    #self.permarray = range(len(elev))
    elev.sort()
    return elev , deg , sen , self.nup

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

  def reset_int(self):
    self.eint , self.eri   = self.transforms()
  
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

  def create_dense_one_body(self ,corrgeq):
    #Create dense matrices (to work with horton) from 1rdm matrices of a Richardson-Gaudin ground state
    dob = DenseOneBody(corrgeq.rgeq.alevel)
    for i in range(corrgeq.rgeq.alevel): 
      dob.set_element(i,i,corrgeq.rdm1[i])
    dob.check_symmetry()
    return dob

  def create_dense_two_body(self,corrgeq):
    #Create dense matrices (to work with horton) from 2rdm matrices of a Richardson-Gaudin ground state
    dtwb = DenseTwoBody(corrgeq.rgeq.alevel)
    for i in range(corrgeq.rgeq.alevel): 
      for j in range(i,corrgeq.rgeq.alevel): 
        for k in range(j,corrgeq.rgeq.alevel): 
          for l in range(k,corrgeq.rgeq.alevel): 
            dtwb.set_element(i,j, k , l,corrgeq.get_element_2rdm(i,j,k,l))
    dtwb.check_symmetry()
    return dtwb

  def set_1rdm_rgeq(self , corrgeq , maxiter = 128 , threshold = 1e-8):
    dmrgeq = self.create_dense_one_body(corrgeq)
    lf = self.horton_ham.system.lf
    wfn = self.horton_ham.system.wfn
    overlap = self.horton_ham.system.get_overlap()
    fock = lf.create_one_body()
    self.horton_ham.compute_fock(fock, None)
    wfn.clear()
    wfn.update_exp(fock, overlap,dm_alpha = dmrgeq)
    #converged = False
    #counter = 0
    #while counter < maxiter:
    #  # Construct the Fock operator
    #  fock.clear()
    #  self.horton_ham.compute_fock(fock, None)
    #  # Check for convergence
    #  error = lf.error_eigen(fock, overlap, wfn.exp_alpha)
    #  if error <threshold:
    #    converged = True
    #    break
    #  # agonalize the fock operator
    #  wfn.clear() # discard previous wfn state
    #  wfn.update_exp(fock, overlap)
    #  # t the hamiltonian know that the wavefunction has changed.
    #  self.horton_ham.clear()
    #  #write intermediate results to checkpoint
    #  self.horton_ham.system.update_chk('wfn')
    #  counter += 1
    self.horton_ham.compute()


class General_Pairing_Ham(Hamiltoniaan):
  """
  The paing Hamiltonian with a level dependend interaction constant.
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
    #print 'The energy = %f' %som
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

def calc():
  chemham = Chem_Ham(filename = None ,atoms = [10], positions = [[0,0,0]],  basis = '3-21G') #3-21G
  ehf = chemham.horton_ham.system.extra['energy']
  #efci =-14.5314443939 #Be fci 3-21G
  efci = -127.9189885238 #get this from gaussian fci ne3-21G
  #efciaugpvdz = -14.618287
  myres = -127.822777
  corp = (ehf -myres)/(  ehf - efci)
  print corp

def play_horton():
  chemham = Chem_Ham(filename = None ,atoms = [8,1,1], positions = [[0,0,0],[0,0,3.809934],[3.688913899461,0.,-0.952633889128]],  basis = 'sto-3g') 
  ehf = chemham.horton_ham.system.extra['energy']
  #efci =-14.5314443939 #Be fci 3-21G
  efci = -127.9189885238 #get this from gaussian fci ne3-21G
  #efciaugpvdz = -14.618287
  myres = -127.822777
  corp = (ehf -myres)/(  ehf - efci)
  print corp

def play_psi():
  name = 'h2test2rdm4'
  #chemham = Psi4_ham('../../CIFlow/rundir/results/beh2_sto_3g_symcompc1/output_run/psi0_sto-3g3.00.mout')
  chemham = Psi4_ham('./results/h26-31g/matrixelements/psi0_6-31g0.60.mout')
  print chemham.get_rgeq_char()
  elev, deg , sen, ap  = chemham.get_rgeq_char()
  g = 0.00001 ; elist = []
  nlevel,deg = rgf.checkdegeneratie(elev ,list(deg))
  sen = nlevel* [0.]
  rgeq = vrg.VarSolver(None).create_rgeq(elev,deg, sen , g , ap, tda = 'weak')
  rgeq.energiel , rgeq.ontaardingen ,rgeq.senioriteit ,rgeq.alevel= rgf.uncheckdegeneratie(rgeq.energiel , rgeq.ontaardingen)
  rgeq.energiel = np.array([-1.34976654, -0.51770481, -0.26592237,  0.33503513 ])
  with open('%s.dat' %name , 'w') as d:
    for g in [ 0.001 * i for i in range(1,117)] + [0.119034351059] + [ 0.001 * i for i in range(121,200)] :
      energierg , rgeq =  wp.generating_datak(rgeq, None,'g',0.001,g,rgwrite = True ,tdafilebool = False,exname = '',moviede = False, intofmotion = False, printstep = 30 )       
      rgeqdm = cf.CorrelationFunction(rgeq, calc_2rdm = True)
      #list2 =rgeqdm.get_1rdm() 
      #print list2 , rgeqdm.rgeq.energy
      #for i in range(len(rgeq.energiel)):
          #for j in range(len(rgeq.energiel)):
            #print i , ' ' , i, ' ' , j , ' ' , j ,' ', rgeqdm.get_element_2rdm(i,i,j,j)
            #print i , ' ' , j, ' ' , i , ' ' , j ,' ',rgeqdm.get_element_2rdm(i,j,i,j)
      elist.append(chemham.calc_energy(rgeqdm))
      #corp = (ehf -elist[-1])/(  ehf - efci)
      d.write('%f\t%f\t%f\n' % (g , elist[-1], rgeqdm.rgeq.energy))
  pl.plot([ 0.001 * i for i in range(1,117)] + [0.119034351059] + [ 0.001 * i for i in range(121,200)], elist)
  pl.savefig('%s.png' % name)

if __name__ == "__main__":
  #play_horton()
  play_psi()
  #calc()
  #chemmain()
  #nucmain()
