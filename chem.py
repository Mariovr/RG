# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Mario Van Raemdonck, 2014;
# (c) Ghent University, 2014
# -*- coding: utf-8 -*-
#! /usr/bin/env python 
from horton import *

import correlationfaribault as cf
import richardsongaudin as rg
import rgfunctions as rgf
import writepairing as wp

import pylab as pl
import sys

def create_rgeq(elev ,deg, sen , g , ap, tda = 'strong'):
  rgeq = rg.RichRedBcs(elev,deg,sen,g,ap)
  if tda == 'strong':
    tdastartd = {0:ap }
    rgeq = rg.RichardsonSolver(rgeq).main_solve(tdastartd,xistep = 0.01,xival = 1.,rgwrite = False,plotrgvarpath = False, plotepath =False,xlim = None , ylim = None)
    
  elif tda == 'weak':
    tdastartd = rgf.tdadict_kleinekoppeling(ap,deg,sen)
    rgeq = rg.RichardsonSolver(rgeq).main_solve(tdastartd,xistep = 0.01,xival = 1.,rgwrite = False,plotrgvarpath = False, plotepath =False,xlim = None , ylim = None)
  else:
    g = rgeq.g ; rgeq.g = 0.0001 ; tdastartd = rgf.tdadict_kleinekoppeling(ap,deg,sen)
    energierg , rgeq =  wp.generating_datak(rgeq,tdastartd,'g',0.001,g,rgwrite = True ,tdafilebool = False,exname = '',moviede = False, intofmotion = False, printstep = 30 )       
  return rgeq

def meanfield(filen = None , basis =  '3-21G'):
  #possible bases:  '3-21G' , STO-3G
  if filen != None:
    # Load the coordinates from file
    sys = System.from_file(filen, obasis= basis)
  else:
    #Or define the system by hand
    #sys = System([[0,0,0], [1.4,0,0]], [1,1] ,  obasis= basis)  #one berryllium atom
    sys = System([[0,0,0]], [10] ,  obasis= basis)  #one berryllium atom
  #print sys.get_kinetic()._array
  #print sys.get_overlap()._array
  #print sys.get_nuclear_attraction()._array
  #print sys.get_electron_repulsion()._array
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

def calc_energy2(kin , nucat , erep , nucnuc , rgcor):
  som = 0
  for i in range(kin.nbasis):
    som += (nucat.get_element(i,i) + kin.get_element(i,i)) * rgcor.rdm1[i] #1rdm is diagonal
  for i in range(kin.nbasis):
    for j in range(kin.nbasis):
      som += 1./2. *erep.get_element(i,j,i,j) *rgcor.get_element_2rdm(i,j,i,j)
      som +=  1./2. *erep.get_element(i,i,j,j) *rgcor.get_element_2rdm(i,i,j,j)
  som += nucnuc
  print '------------------------------------------------'
  print 'The total energy = %f' %som
  print '------------------------------------------------'
  return som

def calc_energy(kin , nucat , erep , nucnuc , rgcor):
  som = 0
  for i in range(kin.nbasis):
    som += 2.*( kin.get_element(i,i)) * rgcor.rdm1[i] #1rdm is diagonal
  for i in range(kin.nbasis):
    for j in range(kin.nbasis):
      if i != j:
        som += (2.*erep.get_element(i,j,i,j) - erep.get_element(i,j,j,i) )*rgcor.get_element_2rdm(i,j,i,j)
        som += erep.get_element(i,i,j,j) *rgcor.get_element_2rdm(i,i,j,j)
    som += erep.get_element(i,i,i,i) *rgcor.get_element_2rdm(i,i,i,i)
  som += nucnuc
  print '------------------------------------------------'
  print 'The total energy = %f' %som
  print '------------------------------------------------'
  return som

def create_densematrices(corrgeq):
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

def transforms(ham):
  ham.system.get_kinetic().iadd(ham.system.get_nuclear_attraction() )
  ham.system.get_kinetic()._array = ham.system.get_kinetic().dot(ham.system.wfn.get_exp('alpha').coeffs.T ,ham.system.wfn.get_exp('alpha').coeffs )
  ham.system.get_electron_repulsion()._array = np.tensordot(ham.system.get_electron_repulsion()._array , ham.system.wfn.get_exp('alpha').coeffs ,axes =  ([0],[0]))
  ham.system.get_electron_repulsion()._array = np.tensordot(ham.system.get_electron_repulsion()._array , ham.system.wfn.get_exp('alpha').coeffs ,axes =  ([0],[0]))
  ham.system.get_electron_repulsion()._array = np.tensordot(ham.system.get_electron_repulsion()._array , ham.system.wfn.get_exp('alpha').coeffs ,axes =  ([0],[0]))
  ham.system.get_electron_repulsion()._array = np.tensordot(ham.system.get_electron_repulsion()._array , ham.system.wfn.get_exp('alpha').coeffs ,axes =  ([0],[0]))
  return ham

def tests():
  ham = meanfield()
  transforms(ham)
  
def main(*args , **kwargs):
  name = 'Ne'
  d = open('%s.dat' %name , 'w')
  ham = meanfield()
  ehf = ham.system.extra['energy']
  efcine = -127.9189885238 #get this from gaussian fci
  efcibe =-14.5314443939 
  transforms(ham)
  nucnuc = ham.system.compute_nucnuc()
  for i in ham.system.wfn._iter_expansions():
    elev = i.energies
    ap = int(sum(i.occupations))
  deg = np.array([2.] *len(elev))
  kin = ham.system.get_kinetic(); nucat = ham.system.get_nuclear_attraction() ; ov = ham.system.get_overlap()
  erep = ham.system.get_electron_repulsion()
  g = 0.0001 ; elist = []
  elev = list(elev)
  nlevel,deg = rgf.checkdegeneratie(elev ,list(deg))
  sen = np.array([0.] *len(elev))
  rgeq = create_rgeq(elev,deg, sen , g , ap, tda = 'weak')
  rgeq.energiel , rgeq.ontaardingen ,rgeq.senioriteit= rgf.uncheckdegeneratie(rgeq.energiel , rgeq.ontaardingen)
  rgeq.alevel = len(rgeq.energiel)
  for g in [ 0.001 * i for i in range(1,201)]:
    rgeqdm = cf.CorrelationFunction(rgeq, calc_2rdm = True)
    elist.append(calc_energy(kin , nucat , erep , nucnuc , rgeqdm) )
    corp = (ehf -elist[-1])/(  ehf - efcine)
    d.write('%f\t%f\t%f\n' % (g , elist[-1], corp))
    energierg , rgeq =  wp.generating_datak(rgeqdm.rgeq, None,'g',0.001,g,rgwrite = True ,tdafilebool = False,exname = '',moviede = False, intofmotion = False, printstep = 30 )       
  pl.plot([ 0.001 * i for i in range(1,201)], elist)
  pl.savefig('%s.png' % name)
  d.close()

if __name__ == "__main__":
  main()
  #tests()

