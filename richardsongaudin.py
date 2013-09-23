# -*- coding: utf-8 -*-
import numpy as np
from numpy import array, zeros, ones, arange,copy
from numpy import append, delete, polynomial #not compatible with pypy 1/10/2012 
from math import sqrt

import copy
from tdasolver import *
import newtonraphson as nr
import rgfunctions as rgf
#source code krylov solver is at scipy/optimize/nonlin.py

class XiError(Exception):
  def __init__(self,xi,energy,rgvars):
    self.xi = xi
    self.energy = energy
    self.rgvars = rgvars
  def __str__():
    strex = 'Exception in RichardsonSolver (XiError) couldn\'t reached the xivalue change_xi asked for\n I posses the following variable: xi (indicates the xi value where breakdown occured), energy (indicates the energy at that xi value), rgvars (indicates the rgvars at that energy)'
    return strex

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
    stringrep = '''#Stringrepresentation of RichardsonEq with subclass: %s
#energylevels: %s
#degeneracys: %s
#senioritys: %s
#the solutions of the set of equations %s
#the interaction constant: %s
#the number of pairs: %s, the xivalue: %s
#the total energy is: %s
''' %(self.__class__.__name__ , str(self.energiel).translate(None,'\n'),str(self.ontaardingen).translate(None,'\n'),str(self.senioriteit).translate(None,'\n'),str(self.rgsolutions).translate(None,'\n'),str(self.g),str(self.apair),str(self.xi),str(self.energy))

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
    return getattr(self,var)

  def get_solutions(self):
    return self.rgsolutions
  
  def setvar(self,var,val):
    setattr(self,var,val)
    
  def test_goodsol(self):
    zerosd = self(self.rgsolutions) 
    print "putting the solutions of the RG variables back in the RG equations gives (needs to be zero):"
    print zerosd
    return sum(abs(test_goodsol))
  
  def solve(self,goodguess = None,tol = 1e-10):
    if goodguess == None:
      self.rgsolutions = nr.solve(self,self.rgsolutions ,tol = tol)
    else:
      self.rgsolutions = nr. solve(self,goodguess,tol = tol)
    self.energy = self.get_energy()
    return self.energy

  def intofmotion(self):
    """
    Calculates the integrals of motion of the Richardson-Gaudin models
    Remark : we need to call the xij function with 2epsilon_j and x_alpha is RG var
    RETURNS: a vector with all the integrals of motion
    """
    di = 1./4.*self.ontaardingen - self.senioriteit *1./2.
    etrans, gtrans = self.iomtrans()
    return di*(-1 + gtrans*array([sum([self.zij(etrans[i],etrans[k])*di[k] for k in range(self.alevel) if k != i]) for i in range(self.alevel)]) + gtrans*array( [sum([self.zij(self.rgsolutions[b],etrans[i]) for b in range(self.apair) ]) for i in range(self.alevel)]))

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
  
  def calc_coef(self, tdasol):
    """
    Used in overlaps.py to calculate the eigenstates but is ugly: ->
    better use the xij and zij functions also to calculate the overlap
    for reasons of simplicity and consistency
    """
    return 1./(2. * self.energiel- tdasol)

  def zij(self, i,j):
    """
    REMARK for the rational Richardson-Gaudin model xij = zij 
    """
    return 1./(i-j)

  def iomtrans(self):
    return 2. * self.energiel , 2.* self.g

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
          updateE = tdasol + 1j*sqrt(2.*self.xi/ai)*rootshermite
          defineupdateE = True
        else:
          updateE = append(updateE,tdasol+ 1j*sqrt(2.*self.xi/ai)*rootshermite)        
    print """the near_tda approximate solutions:
      *******************************************"""
    print updateE
    return updateE
    
  def copy(self):
    #self made copy function because copy.deepcopy() is to damn slow
    if self.rgsolutions == None:
      rgcop = None
    else:
      rgcop = np.copy(self.rgsolutions)
    d = RichRedBcs(self.energiel,self.ontaardingen,self.senioriteit,self.g,self.apair,xi = self.xi,rgsol = rgcop)
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

  def calc_coef(self, tdasol):
    """
    Used in overlaps.py to calculate the eigenstates but is ugly: ->
    better use the xij and zij functions also to calculate the overlap
    for reasons of simplicity and consistency
    """
    return  self.energiel/(self.eta*self.energiel*self.energiel- tdasol)    
    
  def zij(self, i , j):
    return (i+j)/(i-j)

  def iomtrans(self):
    """
    take eta always positive because otherwise some problems arise 
    Transforms the values so they fit in zij and the general formula for the calculation of the integrals of motion.
    """
    return self.energiel*self.energiel*self.eta , 1./(self.eta/self.g + 1 - self.apair + self.alevel / 2. - sum(self.senioriteit)/2.)

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
              updateE = tdasol - 1.j*sqrt(2.*self.xi/(ai*self.eta))*np.sqrt(tdasol)*rootshermite
            else:
              updateE = tdasol + sqrt(2.*self.xi/(ai*self.eta*-1.))*np.sqrt(tdasol)*rootshermite
          else:
            if self.eta > 0 :
              updateE = tdasol + sqrt(2.*self.xi/(ai*self.eta))*np.sqrt(abs(tdasol))*rootshermite
            else:
              updateE = tdasol + 1.j* sqrt(2.*self.xi/(ai*self.eta*-1.))*np.sqrt(abs(tdasol))*rootshermite
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
    s += '#And eta is determined by: %f\n' %self.eta
    return s
  
  def copy(self):
    #self made copy function because copy.deepcopy() is to damn slow
    if self.rgsolutions == None:
      rgcop = None
    else:
      rgcop = np.copy(self.rgsolutions)
    d = RichFacInt(self.energiel,self.ontaardingen,self.senioriteit,self.g,self.eta,self.apair,xi = self.xi,rgsol = rgcop)
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
    if self.richeq.rgsolutions != None:
      self.set_xisolution(self.richeq.get_energy())
    print 'we created a RichardsonSolver with following RG equations: %s' %str(self.richeq)
    
  def get_tda(self,tdadict):
    #if tdadict = None then the tdadict of the self.tda instance is used
    #if that dictionary is not yet changed after initialisation then all tdasolutions are returned
    return self.tda.bisect_tda(tdadict)
  
  def set_tdadict(self,tdad):
    self.tda.set_tdadict(tdad)
  
  def set_xisolution(self,erg):
    self.xisolutions.append((self.richeq.xi,erg,self.richeq.rgsolutions))
    
  def get_xisolutions(idxi = None):
    if xi == None: 
      return self.xisolutions
    else: 
      self.xisorted(self.xisolutions , key = lambda opl: opl[0])
      return self.xisolutions [idxi]
    
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
    #self.tda.tdadict = tdadict
    return tdadict
    
  def change_xi(self,xiarg,crit = 1.8):
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
    testend = False
    while(self.richeq.xi != xiend):
      print 'xi', self.richeq.xi 
      if ((self.richeq.xi + xistep >= xiend and xistep > 0) or (self.richeq.xi +xistep <= xiend and xistep <0 ) ) and testend == False:
        xistep /= 100. ; testend = True ; conarray = []
      rgsol = self.assurecontinuity(conarray,xistep,xiend,crit)
      self.set_xisolution(rgsol) 
    print 'xi is now %f' %self.richeq.xi
    print 'energy is: %f' %rgsol
    return self.richeq
    
  def assurecontinuity(self,carray,xis,xiend,crit):
    '''
    Increases Xi and find the new rg solutions but assure continuity
    REMARK: you can see this function as a wrapper function of both nr.solve and rgf.continuity_check
    '''
    consol = False 
    n = 1
    i = 1
    xisavestep = xis
    smallend = False
    while(consol == False):
      for a in range(n): 
        self.richeq.step_xi(xis)
        if (self.richeq.xi >= xiend and xis > 0) or (self.richeq.xi <= xiend and xis <0 ):
          self.richeq.xi = xiend
        try:
          rgsol = self.richeq.solve()        #MOST IMPORTANT STATEMENT of this function
        except (ValueError, np.linalg.linalg.LinAlgError) as e:
          print e
          self.richeq.energy -= -1e9 #make sure the continuity check fails
      if not rgf.continuity_check(carray,self.richeq ,crit = crit,dvar = 'xi'):
        self.richeq.xi == self.xisolutions[-1][0] ; self.richeq.rgsolutions = self.xisolutions[-1][2]
        print 'Problems with continuity of the energy with changing xi in main_rgsolver()' 
        xis /= 10.**i
        n = 10**i
        i += 1
        if n > 200:
          print 'error we have discoverd a discontinuity in the energy when we variate xi from 0 to 1'
          rgsol = self.circumvent(xisavestep,xiend)
          if not rgf.continuity_check(carray,self.richeq ,crit = crit*5,dvar = 'xi'):
            self.writexipath()
            raise XiError(carray[-1].xi,carray[-1].energy,carray[-1].rgsolutions)
            print 'take another tdasolution'
          else:
            carray.append(copy.deepcopy(self.richeq))
            del(carray[0])
            consol = True
            rgsol = self.richeq.energy
      else:
        consol = True        
    return rgsol 
  
  def circumvent(self,xistep,xiend):
    print 'error by changing xi we try to circumvent the problem regio' #TODO make this function more stable and flexible this is just a quick hack
    gstep = 1.j *abs(self.richeq.xi.real)/1000.   
    cstep = 100
    for i in xrange(cstep):
      print 'making complex %s' %str(self.richeq.xi)
      self.richeq.xi += gstep
      self.richeq.solve()
    for i in xrange(300):
      print 'increasing real part %s'  %str(self.richeq.xi)
      self.richeq.xi += xistep/300.
      if self.richeq.xi.real+abs(xistep)*0.8 >= xiend and self.richeq.xi.real - abs(xistep)*0.8 <= xiend:
        self.richeq.xi = xiend+ gstep*cstep ; self.richeq.solve()
        break
      self.richeq.solve()
    for i in xrange(cstep):
      print 'making back real: %s' %str(self.richeq.xi)
      self.richeq.xi -= gstep
      self.richeq.solve()
    self.richeq.xi = self.richeq.xi.real
    rgsol = self.richeq.solve()
    assert(isinstance(self.richeq.g,float))
    print self.richeq.rgsolutions, self.richeq.xi, self.richeq.get_energy()
    print '!!!!!!!!!!! circumvented a critical point in xispace'
    return rgsol
  
  def writexipath(self,fname = 'xipath',reverse=False):
    '''
    writes the path in xi space of the rgvars in the complexplane and the energy to a file.
    '''
    self.xisolutions = sorted(self.xisolutions , key = lambda opl : opl[0])
    fname = '%s%f%s.dat' % (fname,self.richeq.g,str(self.tda.tdadict).translate(None,' '))
    xifile = open(fname ,'w')
    xifile.write('#This file contains the energy and rgvariables of a set RG equations with variating xi\n')
    xifile.write('%s\n' % str(self.richeq))
    xifile.write('#Xi\tE' + self.richeq.apair*'\trgvar(real)\trgvar(imag)' + '\n')
    if reverse:
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
    self.xisolutions = sorted(self.xisolutions , key = lambda opl : opl[0])
    for i in range(self.richeq.apair):
      try:
        pl.plot(self.xisolutions[0][2][i].real,self.xisolutions[0][2][i].imag,'r.', mfc = 'None' , markersize = 10)
      except KeyError:
        print 'We havent determined the RG solutions with xi =1 yet'
      try:
        pl.plot(self.xisolutions[-1][2][i].real,self.xisolutions[-1][2][i].imag,'g.', markersize = 10)
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
    plotname = '%s%f%s' %(xiname,self.richeq.g,str(self.tda.tdadict).translate(None,' ') )
    if isinstance(self.richeq,RichFacInt):
      plotname += 'eta:'+str(self.richeq.eta)
    pl.savefig(plotname+'.png')
    pl.close()
    
  def plotexi(self,xiname = 'exipath'):
    '''
    Function that plots the energy as function of Xi
    '''
    pl.figure()
    self.xisolutions = sorted(self.xisolutions , key = lambda opl : opl[0])
    xw = [i[0] for i in self.xisolutions]
    yw = [i[1] for i in self.xisolutions]
    pl.plot(xw , yw)  
    pl.xlabel('xi (a.u.)')
    pl.ylabel('energy (a.u.)')
    pl.title('Energy in function of xi')
    pl.savefig('%s%f%s.png' %(xiname,self.richeq.g,str(self.tda.tdadict).translate(None,' ')))
    pl.close()
      
  def main_desolve(self,xistep = -0.01,rgwrite = False,plotrgvarpath = False,plotepath = False,xlim = None, ylim = None,xiend = 0.):
    assert(xistep < 0)
    xiarg = {'xiend' : 1e-6, 'xistep' : xistep} 
    error = False
    self.change_xi(xiarg)
    self.richeq.xi = 0. ; print 'xi', self.richeq.xi
    if isinstance(self.richeq,RichRedBcs):
      rgsol = self.richeq.solve() ; self.set_xisolution(rgsol)
    else :
      rgsol = self.richeq.solve(tol = 1e-4) ; self.set_xisolution(rgsol)
    if rgwrite: self.writexipath()
    if plotrgvarpath: self.plotrgvarsxi(xlim = xlim , ylim = ylim)
    if plotepath: self.plotexi()
    return self.deduceTdaDist() #sets also the tdadict of the self.tda instance to the calculated dictionary of tdasolutions

  def main_solve(self,tdadict,xistep = 0.01 ,xival = 1.,ccrit = 1.7,rgwrite = False,plotrgvarpath = False , plotepath = False,xlim = None, ylim = None):
    print 'we entered main_solve with tdadict: %s' %(str(tdadict))
    assert(xistep > 0 and sum(tdadict.values()) == self.richeq.apair), str(xistep)+str( self.richeq.apair)
    error = False
    self.richeq.xi = 5.992003000000000027e-6
    xiarg = {'xiend' : xival , 'xistep' : xistep}
    rgsol = self.richeq.solve(self.richeq.near_tda(self.get_tda(tdadict),tdadict)) ; self.set_xisolution(rgsol) 
    self.set_tdadict(tdadict)
    try:
      self.change_xi(xiarg,crit = ccrit)
    except (ValueError, np.linalg.linalg.LinAlgError) as e:
      print 'We catched an error in main_solve() of richardsongaudin.py by running self.change_xi(xiarg , crit):' , e
      raise XiError(self.xisolutions[-1][0],self.xisolutions[-1][1],self.xisolutions[-1][2])
    finally:
      if rgwrite: self.writexipath()
      if plotrgvarpath: self.plotrgvarsxi(xlim = xlim , ylim = ylim)
      if plotepath: self.plotexi()
    print  self.richeq.get_energy()
    return self.richeq
    
def maintest():
  '''
  test main for main_rgsolver (which is an extraordinary important function) making that function
  more efficient is a huge timewinst. And making it more flexible will cause much better code
  '''
  #picket fence model
  nlev = 12
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
  g = -1.0000
  #tdastartd = {0:0,1:2,2:0,3:0,4:1,5:0}
  tdastartd = {0:apair }
  alevel = len(eendlev)
  #print ontaardingen,senioriteit,eendlev
  assert(len(ontaardingen) == len(eendlev) and len(ontaardingen) == len(senioriteit))
  rgeq = RichRedBcs(eendlev,ontaardingen,senioriteit,g,apair)
  a = RichardsonSolver(rgeq)
  rgeq = a.main_solve(tdastartd,xistep = 0.01,xival = 1.,rgwrite = True,plotrgvarpath = True , plotepath = True,xlim = None , ylim = None)
  print rgeq.intofmotion()
  #tdad  = a.main_desolve(xistep = -0.01,rgwrite =False,plotrgvarpath = False, plotepath =False,xlim = None , ylim = None,xiend = 0.)
  #print tdad 

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
  
  
