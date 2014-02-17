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

import richardsongaudin as rg
import rgfunctions as rgf
import datareader as dr
import varham as vh
import correlationfaribault as cf
import writepairing as wp

import os
from math import log , sqrt
from numpy import identity,array,dot,zeros,argmax , append

class VarSolver(object):
  """
  The base classe for the variatonal RG solver 
  """
  def __init__(self ,ham ):
    self.ham = ham
    self.rgcor = cf.CorrelationFunction(None, calc_2rdm = True)

  def run(self):
    raise Exception("Implement this function in a subclass from Reader")

  def update_rgeq(self,rgeq):
    self.rgcor.set_rgeq(rgeq, calc_2rdm = True)
  
  def fine_tune(self ,rgeq , val, afhvar = 'g' , tol = 1e-6 , tofile = False):
    if rgeq.getvar(afhvar) >= val:
      tol *= -1.
    else:
      pass
    energierg , rgeq =  wp.generating_datak(rgeq,None,afhvar,tol,val,rgwrite = True ,tdafilebool = False,exname = 'finetune',moviede = False, intofmotion = False, printstep = 30 , tofile = tofile)       
    return rgeq

  def create_rgeq(self, elev ,deg, sen , g , ap, tda = 'strong'):
    """
    Function that creates an rgeq (RichRedBcs) 
    """
    elev = list(elev)
    nlevel,deg = rgf.checkdegeneratie(elev ,list(deg))
    sen = nlevel* [0.]
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

class GoldVarRG(VarSolver):
  def __init__(self, ham , g = -0.001 ,step = -0.001 , end = -0.4 , outputfile = True , nameoutput = 'varcal'):
    """
    Variational RG solver for the interaction constant g with a golden search method.
    """
    super(GoldVarRG,self).__init__(ham )
    elev , deg,sen,ap = self.ham.get_rgeq_char() 
    if not outputfile:
      self.create_output_file(nameoutput, elev , deg , sen,ap , g, step ,  end)
    self.create_reader(name = 'plotenergy.dat')

  def g_opt(self,g, var = 'g', tol = 1e-6):
    self.outpr.readrgvars(afhvar = g , linenr = None ,step = 0.001, startrg =3) 
    rgeq = self.outpr.make_rgeq()
    if tol != None:
      rgeq = self.fine_tune(rgeq , g, afhvar = var , tol = tol)
    self.rgcor.set_rgeq(rgeq, calc_2rdm = True)
    return self.ham.calc_energy(self.rgcor)
 
  def create_output_file(self , name, elev , deg , sen, ap ,start , step ,  end):
    """
    Help function for the nucleon() test hereunder, it generates some RG variables for a desired range of the interaction variable
    in a file called 'plotenergy.dat', together with some plots, and at the end it changes to the directory where all the output is saved.
    """
    elev = list(elev)
    nlevel,deg = rgf.checkdegeneratie(elev ,list(deg))
    sen = nlevel* [0.]
    fname = 'input%s.dat' %name
    with open(fname, 'w') as inputfile:
      inputfile.write("""RichardsonEq : RichRedBcs
energylevels: %s
degeneracies: %s
seniorities: %s 
interaction constant: %f
pairs: %g
change : g
""" %(str(elev).translate(None,',\n'),str(deg).translate(None,',\n'),str(sen).translate(None,',\n') , start ,ap) )
    os.system('python writepairing.py -f %s -r k -x %f -y %f -z %f' %(fname , start  ,step ,end))
    dirname = "%srun_%s_p%gl%g" %('RichRedBcs','k',ap,len(elev)) 
    dirname += fname
    os.chdir(dirname)

  def create_reader(self, name = 'plotenergy.dat'):
    self.outpr = dr.ReaderOutput(name)
    self.outpr.elevels , self.outpr.degeneracies , self.outpr.seniorities, self.outpr.npair  = self.ham.get_rgeq_char()

  def bracket(self,f,x1,h):
    """ 
    Finds the brackets (a,b) of a minimum of the energy of self.ham
    The search starts downhill from x1 with a step length h.
    """
    c = 1.61803398875
    f1 = f(x1, tol =None)
    x2 = x1 + h; f2 = f(x2 , tol = None) #no need for finetuning for the bracket determination
    #Determine downhill direction and change sign of h if needed
    if f2 > f1:
      h = -h
      x2 = x1 + h; f2 =f(x2, tol =None)
    #Check if minimum between x1 - h and x1 + h
      if f2 > f1: return x2,x1 - h
    for i in range (100):
      h = c*h
      x3 = x2 + h; f3 =f(x3, tol =None)
      if f3 > f2: return x1,x3
      x1 = x2; x2 = x3
      f1 = f2; f2 = f3
    print 'Bracket did not find a mimimum'

  def search(self,f,a,b,tol=1.0e-7, filestep = 0.001):
    """      
    (see chapter 10 in Kiusalaas_J._Numerical_Methods_in_Engineering_with_Python_(2005)(en)(424s).pdf)
    Golden section method for determining x that minimizes the energy of self.ham
    The minimum has to be in (a,b).
    """
    nIter = -2.078087*log(tol/abs(b-a)) ; tolv = None #we only start to finetune if our xpoints are closer then 2* file stepwidth
    R = 0.61803398875
    C = 1.0 - R
    # First telescoping
    x1 = R*a + C*b; x2 = C*a + R*b
    f1 =f(x1,tol =tolv ); f2 =f(x2,tol =tolv )
    # Main loop
    for i in range(int(nIter)):
      if f1 > f2:
        a = x1
        x1 = x2; f1 = f2
        x2 = C*a + R*b; f2 =f(x2,tol =tolv )
      else:
        b = x2
        x2 = x1; f2 = f1
        x1 = R*a + C*b; f1 = f(x1,tol =tolv )
      if abs(x1 -x2 ) <= 1. * abs(filestep):
        tolv = tol
    if f1 < f2: return x1,f1
    else: return x2,f2
  
  def run(self ,afhvar = 'g' , start= -0.001  ,width = -0.1 , tol = 1.0e-6, filestep = 0.001):
    a,b = self.bracket(self.g_opt ,start, width) # determines an interval that contains a minimum (warning this could be a local minimum)
    xmin , fmin = self.search(self.g_opt , a,b , tol = tol, filestep = filestep)
    print '------------------------------'
    print 'variational minimum: ', fmin , 'at: ', xmin
    print '------------------------------'
    return xmin , fmin

  #REMARK everything hereunder is still very experimental
  def eopt(self,elevg):
    self.rgcor.rgeq.set_energiel(elevg[0:-1])
    self.rgcor.rgeq.g = elevg[-1]
    self.rgcor.re_calc(calc_2rdm = True)
    return self.ham.calc_energy(self.rgcor)

  def optimize(self, start = -0.001, width = -0.1, tol = 1.0e-5, filestep = 0.001, h = 0.001):
    self.run(start = start ,width = width ,  tol = tol , filestep = filestep)
    return self.powell(self.eopt , append(self.rgcor.rgeq.energiel.copy(),self.rgcor.rgeq.g), h =h)

  def powell(self,F,x,h=0.001,tol=1.0e-2):
    """xMin,nCyc = powell(F,x,h=0.1,tol=1.0e-6)
    Powells method of minimizing F(x).
    x= starting point, h= initial search increment used in bracket
    xMin = mimimum point, nCyc = number of cycles
    """
    def f(s,tol =None): return F(x + s*v) # F in direction of v
    n = len(x) # Number of design variables
    df = zeros((n)) # Decreases of F stored here
    u = identity(n)*1.0 # Vectors v stored here by rows
    for j in range(30): # Allow for 30 cycles:
      xOld = x.copy()
      # Save starting point
      fOld = F(xOld)
      # First n line searches record decreases of F
      for i in range(n):
        v = u[i]
        a,b = self.bracket(f,0.0,h)
        s,fMin = self.search(f,a,b)
        df[i] = fOld - fMin
        fOld = fMin
        x = x + s*v
      #Last line search in the cycle
      v = x - xOld
      a,b = self.bracket(f,0.0,h)
      s,fLast = self.search(f,a,b)
      x = x + s*v
      # Check for convergence
      if sqrt(dot(x-xOld,x-xOld)/n) < tol: return x,j+1
      iMax = int(argmax(df))
      for i in range(iMax,n-1):
        u[i] = u[i+1]
        u[n-1] = v
    print 'Powell did not converge'

def main2():
  chemham = vh.Chem_Ham(filename = None ,atoms = [1,1], positions = [[0,0,0], [1.,0,0]],  basis =  '3-21G') #3-21G
  grgvar = GoldVarRG(chemham,g = 0.0001, step = 0.0002 , end = 0.8, outputfile = False, nameoutput = 'varcal')
  #grgvar.run(start = -0.003 , width = 0.05 , tol = 1.0e-7, filestep = 0.0002)
  xmin , fmin = grgvar.optimize(start = 0.0001 ,width = 0.1,  tol = 1e-5, filestep =0.0002, h = 0.0001)
  print xmin , fmin

def main(*args , **kwargs):
  nuc2 = [[-0.24625, -0.16486, -0.14600, -0.18328, -0.23380],
  [-0.16486, -0.23543, -0.19953, -0.36971, -0.22500],
  [-0.14600, -0.19953, -0.72440, -0.17412, -0.24855],
  [-0.18328, -0.36971, -0.17412, -0.17665, -0.17615],
  [-0.23380, -0.22500, -0.24855, -0.17615, -0.20315]] #the pairing Ham changes this matrix automatically to twofold degenerate levels
  enl=  [-6.1210, -5.5080, -3.8910,  -3.7780 , -3.749] ; deg= [8, 6, 2, 12, 4] ; sen= [0, 0, 0, 0, 0 ] ; ap = 8
  nucham = vh.General_Pairing_Ham(enl , deg, sen, nuc2 , ap)
  grgvar = GoldVarRG(nucham, step = -0.001 , end = -0.6 , outputfile = True)
  #xmin , fmin = grgvar.run(width = -0.1, tol = 1e-6, filestep = 0.001)
  xmin , fmin = grgvar.optimize()
  print xmin , fmin

if __name__ == "__main__":
  #main()
  main2()
