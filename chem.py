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

def corfunc(elev ,deg, sen , g , ap, tda = 'strong'):
  rgeq = rg.RichRedBcs(elev,deg,sen,g,ap)
  if tda == 'strong':
    tdastard = {0:ap }
    rgeq = RichardsonSolver(rgeq).main_solve(tdastartd,xistep = 0.01,xival = 1.,rgwrite = False,plotrgvarpath = False, plotepath =False,xlim = None , ylim = None)
    
  elif tda == 'weak':
    tdastartd = rgf.tdadict_kleinekoppeling(ap,deg,sen)
    rgeq = RichardsonSolver(rgeq).main_solve(tdastartd,xistep = 0.01,xival = 1.,rgwrite = False,plotrgvarpath = False, plotepath =False,xlim = None , ylim = None)
  else:
    g = rgeq.g ; rgeq.g = -0.0001 ; tdastartd = rgf.tdadict_kleinekoppeling(ap,deg,sen)
    energierg , rgeq =  wp.generating_datak(rgeq,tdastartd,'g',-0.001,g,rgwrite = True ,tdafilebool = False,exname = '',moviede = False, intofmotion = False, printstep = 30 )       
  return cf.CorrelationFunction(rgeq)

def main(*args , **kwargs):
  sys = System([[0,0,0]], [4] ,  obasis='3-21G') 
  print sys.coordinates
  kin = sys.get_kinetic()
  for i in range(kin.nbasis):
    for j in range(kin.nbasis):
      print kin.get_element(i,j)
  nucat = sys.get_nuclear_attraction()
  erep = sys.get_electron_repulsion()
  pass

if __name__ == "__main__":
  main()

