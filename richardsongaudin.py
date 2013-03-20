# -*- coding: utf-8 -*-
import math
import numpy as np
from numpy import array, zeros, ones, arange,copy
from numpy import append, delete, polynomial # here are the numpy imports that are at this time 23/08/2012 not compatibel with pypy
import inspect

from tdasolver import *
import newtonraphson as nr
import rgfunctions as rgf
#source code krylov solver is at scipy/optimize/nonlin.py

class RichardsonEq(object): 
  def __init__(self,energiel_,ontaard_,senior_, kop_,apair_, xi = 5.992003000000000027e-6,rgsol = None):
    """
    the start of a general Richardson-Gaudin problem, with XI set in the initialisation to a value that is a bit larger than zero
    (XI = zero corresponds to the tda approximation, XI = one corresponds to the pairing problem that we need to solve)
    """
    self.xi = xi #just a random very small variable so that the near tda approxiamation works
    self.energiel = array(sorted(energiel_))
    self.ontaardingen = array(ontaard_)
    self.senioriteit = array(senior_)
    self.g = kop_
    self.apair = apair_
    self.alevel = len(self.energiel)
    self.rgsolutions = rgsol
    self.energy = None
    assert(len(self.senioriteit) == len(self.ontaardingen) == len(self.energiel))
    
  def __str__(self):
    stringrep = '''#Stringrepresentation of RichardsonEq:
#energylevels: %s
#degeneracys: %s
#senioritys: %s
#the solutions of the set of equations %s
#the interaction constant: %s
#the number of pairs: %s, the xivalue: %s
''' %(str(self.energiel).translate(None,'\n'),str(self.ontaardingen).translate(None,'\n'),str(self.senioriteit).translate(None,'\n'),str(self.rgsolutions).translate(None,'\n'),str(self.g),str(self.apair),str(self.xi))

    return stringrep
  
  #def __del__(self):
    #print 'The following object is destroyed %s' %str(self)
    
  def step_xi(self,step = 1e-2):
    self.xi += step 
    #it's no use to go bigger then one and smaller then 0 because the case that interests us is xi =1
    if self.xi > 1.:
      self.xi = 1.
    if self.xi < 0.:
      self.xi = 0.
  
  def set_g(self,kop):
    self.g = kop
  
  def set_apair(self,ap):
    self.apair = ap
 
  def getvar(self,var):
    #handy to create general programs
    return dict(inspect.getmembers(self))[var]
  
  def setvar(self,var,val):
    if var == 'g':
      self.g = val
    elif var == 'eta':
      assert(isinstance(self,RichFacInt))
      self.eta = val
    elif var == 'xi':
      self.xi = val
    else:
      print 'val is not a member of RichardsonEq'
    
  def test_goodsol(self):
    zerosd = self(self.rgsolutions) 
    print "putting the solutions of the RG variables back in the RG equations gives (needs to be zero):"
    print zerosd
    return sum(abs(test_goodsol))
  
  def solve(self,goodguess = None,tol = 1e-11):
    if goodguess == None:
      self.rgsolutions = nr.solve(self,self.rgsolutions ,tol = tol)
    else:
      self.rgsolutions = nr. solve(self,goodguess,tol = tol)
    self.energy = self.get_energy()
    return self.energy

    
class RichRedBcs(RichardsonEq):
  """
  Class that defines the general Richardson-Gaudin equations produced as bethe ansatz equations
  for the reduced BCS Hamiltonian. (with the pseudo deformation parameter XI that defines the algebra)  
  """ 
  def __call__(self,x):
    rgfunc = zeros(self.apair,complex)
    for i in xrange(self.apair):
      rgfunc[i] = 1.+self.g* sum( (1./2.*self.ontaardingen - self.xi*self.senioriteit)/(2.*self.energiel-x[i]))
      for k in xrange(self.apair):
	if k != i:
	  rgfunc[i] -= 2.*self.g *self.xi* 1./(x[k]-x[i])
    return rgfunc
    
  def jacobian(self,x):
    jac = zeros((self.apair,self.apair),complex)
    for i in xrange(self.apair):
      for j in xrange(self.apair):
	if i == j: 
	  jac[i,i] = -self.g*sum((self.xi*self.senioriteit-1./2.*self.ontaardingen)/(2.*self.energiel-x[i])**2)
	  for k in xrange(self.apair):
	    if k != i:
	      jac[i,i] -= 2.*self.g *self.xi* 1./(x[k]-x[i])**2
	else:  jac[i,j] = 2.*self.g*self.xi/(x[j]-x[i])**2
    return jac  
    
  def get_energy(self):
    return sum(self.rgsolutions.real)
    
  def near_tda(self,tdasolreal,tdadict):
    """
    gives the solution of the RG equations when XI is very small see the notes of S. De Baerdemacker
    by the usage of hermite polynomials ( this is necessary because if the tda start solutions are the same , there are singularities
    in the RG equations. Returns the updated general Richardson-Gaudin variables
    """
    defineupdateE = False
    for key in tdadict.keys():
      if tdadict[key] != 0.:
	tdasol = ones(tdadict[key],int)*tdasolreal[key]
	ai = 1./2.*sum(self.ontaardingen/(2.*self.energiel-tdasol[0])**2.)
	n = len(tdasol)
	coeff = zeros(n+1,float)
	coeff[n] = 1
	rootshermite = polynomial.hermite.hermroots(coeff)
	if defineupdateE == False:
	  updateE = tdasol + 1j*math.sqrt(2.*self.xi/ai)*rootshermite
	  defineupdateE = True
	else:
	  updateE = append(updateE,tdasol+ 1j*math.sqrt(2.*self.xi/ai)*rootshermite)	
    print """the near_tda approximate solutions:
      *******************************************"""
    print updateE
    return updateE
    
  def copy(self):
    #self made copy function because copy.deepcopy() is to damn slow
    d = RichRedBcs(self.energiel,self.ontaardingen,self.senioriteit,self.g,self.apair,xi = self.xi,rgsol = copy(self.rgsolutions))
    d.energy = self.energy
    return d

class RichFacInt(RichardsonEq):
  """
  Class that defines the bethe ansatz equations for a factorisable interaction (see Stijn notes) or the
  article of Stefan Rombouts, Jorge Dukelsky and Gerardo Ortiz; Quantum phase diagram of the integrable px+ipy fermionic superfluid; Phys.Rev.B; (2010)
  """
  def __init__(self,energiel_,ontaard_,senior_, koppelingsconstante,eta_,apair_ , xi = 5.992003000000000027e-6,rgsol = None):
    """
    the start of a general Richardson-Gaudin problem, with XI set in the initialisation to a value that is a bit larger than zero
    (XI = zero corresponds to the tda approximation, XI = one corresponds to the pairing problem that we need to solve)
    """
    self.eta = eta_ #make sure you set first the not inherited variables because the initialisation function of RichardsonEq contains solve and that function makes use of all the variables 
    super(RichFacInt,self).__init__(energiel_,ontaard_,senior_, koppelingsconstante,apair_ ,xi = xi,rgsol = rgsol)
  
  def __call__(self,x):
    rgfunc = zeros(self.apair,complex)
    for i in xrange(self.apair):
      rgfunc[i] = 1.+self.g *sum( (1./2.*self.ontaardingen - self.xi*self.senioriteit)*self.energiel*self.energiel/(self.eta*self.energiel*self.energiel-x[i]))
      for k in xrange(self.apair):
	if k != i:
	  rgfunc[i] -= 2.*self.g *self.xi/self.eta* x[k]/(x[k]-x[i])
    return rgfunc
    
  def jacobian(self,x):
    jac = zeros((self.apair,self.apair),complex)
    for i in xrange(self.apair):
      for j in xrange(self.apair):
	if i == j: 
	  jac[i,i] = -self.g*sum((self.xi*self.senioriteit-1./2.*self.ontaardingen)*self.energiel*self.energiel/(self.eta*self.energiel*self.energiel-x[i])**2)
	  for k in xrange(self.apair):
	    if k != i:
	      jac[i,i] -= 2.*self.g *self.xi/self.eta* x[k]/(x[k]-x[i])**2
	else:  jac[i,j] = 2.*self.g*self.xi/self.eta*x[i]/(x[j]-x[i])**2
    return jac
    
  def get_energy(self):
    return sum(self.rgsolutions.real)-self.eta*sum(self.energiel*self.energiel*(1./4.*self.ontaardingen-1./2.*self.senioriteit))
    
  def near_tda(self,tdasolreal,tdadict):
    """
    gives the solution of the RG equations when XI is very small see the notes of S. De Baerdemacker
    by the usage of hermite polynomials ( this is necessary because if the tda start solutions are the same , there are singularities
    in the RG equations. Returns the updated general Richardson-Gaudin variables
    """
    defineupdateE = False
    for key in tdadict.keys():
      if tdadict[key] != 0.:
	tdas = tdasolreal[key]
	tdasol = ones(tdadict[key],float)*tdas
	ai = 1./2.*sum(self.ontaardingen*self.energiel*self.energiel/(self.eta*self.energiel*self.energiel-tdas)**2.)
	n = len(tdasol)
	coeff = zeros(n+1,float)
	coeff[n] = 1
	rootshermite = polynomial.hermite.hermroots(coeff)
	if defineupdateE == False:
	  if tdas > 0:
	    if self.eta > 0 :
	      updateE = tdasol - 1.j*math.sqrt(2.*self.xi/(ai*self.eta))*np.sqrt(tdasol)*rootshermite
	    else:
	      updateE = tdasol + math.sqrt(2.*self.xi/(ai*self.eta*-1.))*np.sqrt(tdasol)*rootshermite
	  else:
	    if self.eta > 0 :
	      updateE = tdasol + math.sqrt(2.*self.xi/(ai*self.eta))*np.sqrt(abs(tdasol))*rootshermite
	    else:
	      updateE = tdasol + 1.j* math.sqrt(2.*self.xi/(ai*self.eta*-1.))*np.sqrt(abs(tdasol))*rootshermite
	  defineupdateE = True
	else:
	  if tdas > 0:
	    if self.eta >0:
	      updateE = append(updateE,tdasol- 1.j*np.sqrt(2.*self.xi*tdasol/(ai*self.eta))*rootshermite)
	    else:
	      updateE = append(updateE,tdasol+ np.sqrt(2.*self.xi*tdasol/(ai*self.eta*-1.))*rootshermite)
 	  else:
	    if self.eta >0:
	      updateE = append(updateE,tdasol+ np.sqrt(2.*self.xi*abs(tdasol)/(ai*self.eta))*rootshermite)
	    else:
	      updateE = append(updateE,tdasol+ 1.j *np.sqrt(2.*self.xi*abs(tdasol)/(ai*self.eta*-1.))*rootshermite)
	      
    print """the near_tda approximate solutions:
      *******************************************"""
    print updateE
    return updateE  
    
  def __str__(self):
    s = super(RichFacInt,self).__str__()
    s += '\n#And eta is determined by %f' %self.eta
    return s
  
  def copy(self):
    #self made copy function because copy.deepcopy() is to damn slow
    d = RichFacInt(self.energiel,self.ontaardingen,self.senioriteit,self.g,self.eta,self.apair,xi = self.xi,rgsol = copy(self.rgsolutions))
    d.energy = self.energy
    return d
    
class RichardsonSolver(object):
  def __init__(self,richeq):
    '''
    initialisator of solver that changes the Richardson-Gaudin equations
    it's main attributes are a RichardsonEq and a TdaSolver
    '''    
    assert(isinstance(richeq,RichardsonEq))
    self.richeq = richeq.copy()
    if isinstance(richeq,RichFacInt):
      self.tda = TdaFacInt(self.richeq.energiel,self.richeq.ontaardingen,self.richeq.g,self.richeq.eta)
    elif isinstance(richeq,RichRedBcs):
      self.tda = TdaRedBcs(self.richeq.energiel,self.richeq.ontaardingen,self.richeq.g)
    self.xisolutions = []# use as [xi: (E, Rgvar), ...]
    if self.richeq.rgsolutions.dtype != 'object':
      self.set_xisolution(self.richeq.get_energy())
    print 'we created a RichardsonSolver with following RG equations: %s' %str(self.richeq)
    
  def get_tda(self,tdadict):
    return self.tda.bisect_tda(tdadict)
  
  def change_tdadict(self,tdad):
    self.tda.set_tdadict(tdad)
  
  def set_xisolution(self,erg):
    self.xisolutions.append((self.richeq.xi,erg,self.richeq.rgsolutions))
    
  def get_xisolutions(idxi = None):
    if xi == None: 
      return self.xioplossingen
    else: 
      sorted(self.xioplossingen , key = lambda opl: opl[0])
      return self.xioplossingen[idxi]
    
  def deduceTdaDist(self):
    '''
    Function that determines the tda distribution corresponding with the Richardson-Gaudin solution (rgvar) of the Richardson-Gaudin
    variables at Xi = 0
    '''
    tdadict = np.zeros(self.richeq.alevel)
    assert(self.richeq.xi <= 0.01)
    realrgvar = self.richeq.rgsolutions.real
    realrgvar.sort()
    assert(self.tda.g is not 0.)
    i ,j  = 0,0 
    while j != len(realrgvar):
      if self.tda.g < 0 and i is 0:
	x1 = self.tda.singularitys[0]- abs( self.tda.get_extremeboundaries())*2.
	x2 = self.tda.singularitys[0]
      elif i != 0:
	try:
	  x1 = self.tda.singularitys[i-1]
	  x2 = self.tda.singularitys[i]
	except IndexError:  
	  assert(self.tda.g > 0)
	  x1 = self.tda.singularitys[i-1]
	  x2 = self.tda.singularitys[i-1]+ abs(self.tda.get_extremeboundaries()) *2.    
      else: #handles case i ==0 and g>0
	i +=1
	continue
      if x1 < realrgvar[j] and  x2 > realrgvar[j]:
	if self.tda.g > 0:
	  tdadict[i-1] += 1
	else:
	  tdadict[i] += 1
	j += 1
      else:
	i+=1 
    assert(sum(tdadict) == len(realrgvar))
    self.tda.tdadict = tdadict
    return tdadict
    
  def change_xi(self,xiarg,crit = 2.):
    '''
    function that changes the Xi value of an rgequation such as xiarg defines
    the xiarg argument is a dictionary that contains: three keys: 'xibegin','xiend','xistep' it's meaning should be obvious
    '''
    print 'we entered change_xi'
    xistep = xiarg['xistep']
    xiend = xiarg['xiend']
    if xiend > self.richeq.xi and xistep < 0: xistep *= -1.
    elif xiend < self.richeq.xi and xistep > 0: xistep *= -1.
    elif xiend == self.richeq.xi: return self.richeq
    conarray = []
    while(self.richeq.xi != xiend):
      print 'xi', self.richeq.xi 
      rgsol = self.assurecontinuity(conarray,xistep,xiend)
      self.set_xisolution(rgsol) 
    print 'xi is now %f' %self.richeq.xi
    print 'energy is: %f' %rgsol
    return self.richeq
    
  def assurecontinuity(self,carray,xis,xiend):
    '''
    Increases Xi and find the new rg solutions but assure continuity
    REMARK: you can see this function as a wrapper function of both nr.solve and rgf.continuity_check
    '''
    consol = False
    carray2 = list(carray)
    n = 1
    i = 1
    xinewstart = self.richeq.xi
    while(consol == False):
      rgsolsave = self.richeq.rgsolutions
      for a in range(n):
	self.richeq.xi += xis
	if self.richeq.xi+abs(xis)/2. >= xiend and self.richeq.xi - abs(xis)/2. <= xiend:
	  self.richeq.xi = xiend
	  conarray = []
	rgsol = self.richeq.solve()	#MOST IMPORTANT STATEMENT of this function
      if rgsol == None: #to circumvent the special case that rgsol is None (sum(rgsol.real)) so the ValueError get's catched by wrappers of rg_mainsolver
	raise ValueError
      if not rgf.continuity_check(carray2,self.richeq ,crit = 1.5):
	self.richeq.rgsolutions = rgsolsave
	print 'Problems with continuity of the energy with changing xi in main_rgsolver()' 
	self.richeq.xi = xinewstart
	xis /= 10.**i
	n = 10**i
	i += 1
	if n > 201:
	  print 'error we have discoverd a discontinuity in the energy when we variate xi from 0 to 1'
	  sys.exit()
      else:
	consol = True	
    return rgsol 
  
  def circumvent(self,xistep):
    print 'error by changing xi we try to circumvent the problem regio' #TODO make this function more stable and flexible this is just a quick hack
    self.richeq.rgsolutions = self.xisolutions[-1][2]
    self.richeq.xi -= xistep
    self.richeq.g += 1.j *abs(self.richeq.g.real)/100.
    self.richeq.solve()
    self.richeq.xi += xistep
    self.richeq.solve()
    self.richeq.g = self.richeq.g.real
    rgsol = self.richeq.solve()
    assert(isinstance(self.richeq.g,float))
    print self.richeq.rgsolutions, self.richeq.xi, self.richeq.get_energy()
    print '!!!!!!!!!!! circumvented a critical point in xispace'
    return rgsol
  
  def writexipath(self,fname = 'xipath',reverse=False):
    '''
    writes the path in xi space of the rgvars in the complexplane and the energy to a file.
    '''
    sorted(self.xisolutions , key = lambda opl : opl[0])
    fname = '%s%f.dat' % (fname,self.richeq.g)
    xifile = open(fname ,'w')
    xifile.write('#This file contains the energy and rgvariables of a set RG equations with variating xi\n')
    xifile.write('%s\n' % str(self.richeq))
    xifile.write('#Xi\tE' + self.richeq.apair*'\trgvar(real)\trgvar(imag)' + '\n')
    sorted(self.xisolutions , key = lambda opl : opl[0] , reverse = reverse)
    for opl in self.xisolutions:
      xifile.write('%f\t%f\t' %(opl[0] , opl[1]))
      for j in xrange(self.richeq.apair):
	xifile.write('%f\t%f\t' %(opl[2][j].real, opl[2][j].imag ))
      xifile.write('\n')
    xifile.close()
    
  def plotrgvarsxi(self, xiname = 'rgvarxipath' ,xlim= True, ylim = (-20,20)):
    '''
    Function that plots the rgvars in the complexplane as XI changes
    '''
    singularitys = self.tda.calc_singularity()
    pl.figure()
    sorted(self.xisolutions , key = lambda opl : opl[0])
    for i in range(self.richeq.apair):
      try:
        pl.plot(self.xisolutions[0][2][i].real,self.xisolutions[0][2][i].imag,'g.', markersize = 10)
      except KeyError:
	print 'We havent determined the RG solutions with xi =1 yet'
      try:
        pl.plot(self.xisolutions[-1][2][i].real,self.xisolutions[-1][2][i].imag,'r.',mfc = 'None', markersize = 10)
      except KeyError:
	print 'the tdasolutions are not yet determined with change_xi'
    for i in range(self.richeq.apair):
      xw = [opl[2][i].real for opl in self.xisolutions]
      yw = [opl[2][i].imag for opl in self.xisolutions] 
      pl.plot(xw, yw, 'b-')       
    for i in range(self.richeq.alevel):
      pl.axvline(x = singularitys[i] ,c=  'k',linestyle = '--')
    pl.xlabel('real part of rgvars (a.u)')
    pl.ylabel('imaginary part of rgvars (a.u.)')
    pl.title('Richardson-Gaudin variables at g = %f (xi in [0,1])' %(self.richeq.g))
    if xlim != None: pl.xlim((singularitys[0]-(singularitys[-1]-singularitys[0])/2.,singularitys[-1]+(singularitys[-1]-singularitys[0])/2.))
    if ylim != None: pl.ylim(ylim)
    pl.savefig('%s%f.png' %(xiname,self.richeq.g ))
    pl.close()
    
  def plotexi(self,xiname = 'exipath'):
    '''
    Function that plots the energy as function of Xi
    '''
    pl.figure()
    sorted(self.xisolutions , key = lambda opl : opl[0])
    xw = [i[0] for i in self.xisolutions]
    yw = [i[1] for i in self.xisolutions]
    pl.plot(xw , yw)  
    pl.xlabel('xi (a.u.)')
    pl.ylabel('energy (a.u.)')
    pl.title('Energy in function of xi')
    pl.savefig('%s%f.png' %(xiname,self.richeq.g))
    pl.close()
      
  def main_desolve(self,xistep = -0.01,rgwrite = False,plotrgvarpath = False,plotepath = False,xlim = None, ylim = None,xiend = 0.):
    assert(xistep < 0)
    xiarg = {'xiend' : 1e-6, 'xistep' : xistep}
    self.change_xi(xiarg)
    self.richeq.xi = 0. ; print 'xi', self.richeq.xi
    if isinstance(self.richeq,RichRedBcs):
      rgsol = self.richeq.solve() ; self.set_xisolution(rgsol)
    else :
      rgsol = self.richeq.solve(tol = 1e-6) ; self.set_xisolution(rgsol)
    if rgwrite: self.writexipath()
    if plotrgvarpath: self.plotrgvarsxi(xlim = xlim , ylim = ylim)
    if plotepath: self.plotexi()
    return self.deduceTdaDist()

  def main_solve(self,tdadict,xistep = 0.01,xival = 1.,ccrit = 2.,rgwrite = False,plotrgvarpath = False , plotepath = False,xlim = None, ylim = None):
    print 'we entered main_solve'
    assert(xistep > 0 and sum(tdadict.values()) == self.richeq.apair), str(xistep)+str( self.richeq.apair)
    self.richeq.xi = 5.992003000000000027e-6
    xiarg = {'xiend' : xival , 'xistep' : xistep}
    self.richeq.solve(self.richeq.near_tda(self.get_tda(tdadict),tdadict))
    self.change_xi(xiarg,crit = ccrit)
    if rgwrite: self.writexipath()
    if plotrgvarpath: self.plotrgvarsxi(xlim = xlim , ylim = ylim)
    if plotepath: self.plotexi()
    return self.richeq.get_energy(),self.richeq
    
def maintest():
  '''
  test main for main_rgsolver (which is an extraordinary important function) making that function
  more efficient is a huge timewinst. And making it more flexible will cause much better code
  '''
  #eendlev = [9.8696,39.4784,88.8264,157.914,246.74,355.306,483.611,631.655,799.438,986.96]
  eendlev = np.arange(1,13)
  #ontaardingen = [6, 8, 2, 4, 12, 8]
  ontaardingen = np.ones(12,float)*2
  #senioriteit = [3, 1, 1, 0, 2, 1]
  senioriteit = np.zeros(12,float)
  tdastartd = {0:6,1: 0,2:0,3:0,4:0,5:0 ,8:0}
  #tdastartd = {0:6}
  g = -1.
  #apair = 10
  apair = 6
  eta = 1.
  xiarg = {'xistep' : 0.01 , 'xiend': 1.}
  #print ontaardingen,senioriteit,eendlev
  assert(len(ontaardingen) == len(eendlev) and len(ontaardingen) == len(senioriteit))
  rgeq = RichFacInt(eendlev,ontaardingen,senioriteit,g,eta,apair)
  a = RichardsonSolver(rgeq)
  energie,rgeq = a.main_solve(tdastartd,xistep = 0.01,xival = 1.,rgwrite = False,plotrgvarpath = False , plotepath = False,xlim = None , ylim = None)
  tdad  = a.main_desolve(xistep = -0.01,rgwrite = True,plotrgvarpath = True , plotepath = True,xlim = None , ylim = None,xiend = 0.)
  print tdad 

def test_copy():
  eendlev = np.arange(1,13) ; ontaardingen = np.ones(12,float)*2 ; senioriteit = np.zeros(12,float)
  g = -1. ; apair = 6 ; eta = 1.
  rgeq = RichRedBcs(eendlev,ontaardingen,senioriteit,g,apair)
  d = rgeq.copy()
  d.xi = 4.
  print  rgeq.xi ,d.xi 
  d.rgsolutions = np.ones(6,float)
  print rgeq.rgsolutions , d.rgsolutions
  d.energiel = np.zeros(12,float)
  print rgeq.energiel , d.energiel
  d.g = -10.
  print rgeq.g , d.g
  
if __name__ == "__main__":
  maintest()
  #test_copy()
  
  