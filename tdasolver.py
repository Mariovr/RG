# -*- coding: utf-8 -*-
#! /bin/env python
import math, sys
import numpy as np
import pylab as pl

class TdaSolver(object):
  '''
  Baseclass that solves and is able to plot all kinds of different tda equations
  '''
  def __init__(self,eendlevi,ontaardingi,gi,tdadic = None):
    '''
    energiel contains the sp levels; ontaardingi contains the degeneracys; gi contains the interaction constant ; tdadic contains 
    a dictionary that determines the solutions that need to be found.
    '''
    self.energiel = np.array(eendlevi)
    self.nlevel = len(eendlevi)
    self.ontaardingen = np.array(ontaardingi)
    self.g = gi
    self.tdasolutions = np.zeros(self.nlevel,np.float)
    if tdadic is None:
      self.tdadict = dict(zip(range(self.nlevel),[1]*self.nlevel))
    else:
      self.tdadict = tdadic
    self.singularitys = self.calc_singularity()
    
  def set_tdadict(self,tdadictionary):
    self.tdadict = tdadictionary
  
  def calc_singularity(self):
    return None

  def getvar(self,var):
    #handy to create general programs
    return getattr(self,var)

  def setvar(self,var,val):
    setattr(self,var,val)

  def get_solutions(self):
    return self.tdasolutions
  
  def set_g(self,gi):
    self.g = gi

  def bisect_tda(self,tdadict = None):
    """
    bisection method to find the zeros of a particular tda equation,
    remark: the poles are 2 times the single particles energies
    remark: we added an option tdadict for efficiency if we want multiple richardson-gaudin solutions for a large system
    that only needs a few tda solutions (one) because of strong interaction constant,
    so only the tda solutions who correspond to the dictionary tdadic are calculated
    remark: this function can only be called from a subclass of tdasolver that implements tda_vgl
    """
    nlevel = len(self.energiel)
    if tdadict == None:
      tdadict = self.tdadict
    self.tdadict = tdadict #to make sure that the tdadict the tdasolve object has is consistent with the tdasolutions vector
    nauwkeurigheid = 1e-14
    verschil =1.
    maxnumbert = 1000 #max number of tries the bisection method searches for a solution in a given interval
    tel = 0
    for key in tdadict.keys():   
      if key == 0 and tdadict[0] is not 0 and self.g < 0:
        # we put x1 equals to a lower bound for the ground state energie solution of the tda equations (see stijn notes)
        x1 = self.singularitys[0]+ self.get_extremeboundaries()
        x2 = self.singularitys[0]
      elif key == nlevel-1 and tdadict[key] is not 0 and self.g > 0: #make sure this works with positive g as well
        #x2 needs to be the value where the tda_vgl is biggest and x1 where the tda_vgl is lowest if g>0 the tdafunction between the poles decreases so x2 becomes the smaller point and x1 the bigger
        x2 = self.singularitys[key]
        x1 = self.singularitys[key]+ self.get_extremeboundaries()
      elif tdadict[key] is not 0 and self.g < 0:
        x1 = self.singularitys[key-1]
        x2 = self.singularitys[key]
      elif tdadict[key] is not 0 and self.g > 0:
        #x2 needs to be the value where the tda_vgl is biggest and x1 where the tda_vgl is lowest if g>0 the tdafunction between the poles decreases so x2 becomes the smaller point and x1 the bigger
        x2 = self.singularitys[key]
        x1 = self.singularitys[key+1]
      else:  
        continue
      numbert = 0
      while(abs(verschil) > nauwkeurigheid):
        midden = (x1+x2)/2.
        verschil = self.tda_vgl(midden)
        if  verschil > 0.:
          x2 = midden
        else:
          x1 = midden
        numbert += 1
        if numbert > maxnumbert:
          print 'tdasolutions didn\'t converge'
          break
      self.tdasolutions[key] = midden
      verschil = 1.
    print 'bisect heeft volgende tda_oplossingen gevonden:%s , bij volgende tdadict: %s' %(str(self.tdasolutions),str(tdadict)) 
    return self.tdasolutions
  
  def get_nauwkeurigheid(self):
    nauwkeurigheid = []
    for i in self.tdasolutions:
      nauwkeurigheid.append(self.tda_vgl(i))
    return nauwkeurigheid
    
  def __str__(self):
    stringrep = '''#stringrepresentatie tdasolver:
#_________________________________________________
#energiewaarden zijn: %s
#ontaardingen zijn: %s
#de koppelingsconstante is: %s
#invullen van volgende oplossingen: %s
#(zou nul moeten geven) geeft: %s''' %(str(self.energiel),str(self.ontaardingen),str(self.g),str(self.tdasolutions),str(self.get_nauwkeurigheid()))
    return stringrep
    
  def plot_tdaopl(self, solutions = True):
    """
    plot the graphical solution method of a particular tda equation
    and the corresponding solutions (if solutions = true). 
    """
    fig1=pl.figure(1)
    #plot first the line -1/g
    hoogtehl = -1./self.g
    plot1 = pl.axhline(y = hoogtehl , c = 'r',label = 'y = -1/g')
    #plot then the summation term of the tda equation relevant for the problem in consideration (xlim and ylim good so we see usefull stuf)
    einde = self.singularitys[-1]+abs(self.get_extremeboundaries())
    if solutions is True: 
      begin = self.tdasolutions[0] - (self.singularitys[-1] - self.singularitys[0])
    else: 
      begin = self.energiel[0] - (self.energiel[-1] - self.energiel[0]) * self.g**2
    xw = np.arange(begin,einde,0.01)
    fxw = [(-self.tda_vgl(x)-1)/self.g for x in xw]
    plot2 = pl.plot(xw,fxw,'b',label = 'sommatieterm TDA sec. vgl.')
    pl.ylim(ymax = -1./(self.g)*3./2. , ymin = 1./(2.*self.g))
    pl.xlabel(r'$E_{\alpha}$ (a.e.)')
    pl.ylabel('de twee delen van de TDA sec. vgl.')
    for i in xrange(self.nlevel):    
      if solutions is True:
        if i == 0:
          plot3 = pl.plot(self.tdasolutions[i],hoogtehl,'g.',markersize=20.0,label ='oplossingen' )
        else:
          pl.plot(self.tdasolutions[i],hoogtehl,'g.',markersize=20.0 )
      pl.axvline(x = self.singularitys[i] ,c=  'k')
      pl.text(self.singularitys[i], self.g*0.5 ,"o:%g"%self.ontaardingen[i])
    
    pl.legend(loc = 'lower left', numpoints=1)
    pl.savefig('TDA.png')
    pl.show()
    
class TdaRedBcs(TdaSolver):
  '''
  Subclass of Tdasolver that generates the tda equations that arise as the bosonic limit of the RG equations of the
  reduced BCS Hamiltonian
  '''
  def tda_vgl(self,x):
    """
    returns the value of the tda eq. defined by 
    the single particle levels("energiel") degeneracy("ontaarding") , interaction constant("g) and x value.
    The solutions of this eqution are the values x that sets the equation underneath to zero.
    """
    opl = -1. -sum(self.g/2.*self.ontaardingen/(2.*self.energiel-x) )
    return opl
    
  def calc_singularity(self):
    return 2.*self.energiel
  
  def get_extremeboundaries(self):
    return self.g/2.*sum(self.ontaardingen)

  def calc_coef(self, tdasol):
    return 1./(2. * self.energiel - tdasol)

class TdaFacInt(TdaSolver):
  '''
  Subclass of Tdasolver that generates and solves the equations that arise as the bosonic limit of the RG equations
  of the factorisable interaction
  '''
  def __init__(self,energieli,ontaardingi,gi,etai,tdadic = None):
    self.eta = etai
    super(TdaFacInt,self).__init__(energieli,ontaardingi,gi,tdadic = tdadic)
 
  def tda_vgl(self,x):
    """
    returns the value of the tda eq. defined by 
    the single particle levels("energiel") degeneracy("ontaarding") , interaction constant("g) and x value.
    The solutions of this eqution are the values x that sets the equation underneath to zero.
    """
    opl = -1.-self.g/2.*sum(self.energiel*self.energiel*self.ontaardingen/(self.eta*self.energiel*self.energiel-x))
    return opl
  
  def calc_singularity(self ):
    return self.eta*self.energiel*self.energiel 
  
  def get_extremeboundaries(self):
    return self.g*sum(self.ontaardingen*self.energiel*self.energiel)/2.

  def calc_coef(self, tdasol):
    return  self.energiel/(self.eta*self.energiel *self.energiel - tdasol)    
  
  def __str__(self):
    s = super(TdaFacInt,self).__str__()
    s += '\n#En eta wordt gegeven door: %f' %self.eta 
    return s

class TdaDicke(TdaSolver):
  '''
  Subclass of Tdasolver that generates and solves the equations that arise as the bosonic limit of the RG equations
  of the Dicke model
  '''
  def __init__(self,energieli,ontaardingi,gi,epsilon0,tdadic = None):
    self.epsilon0 = epsilon0
    super(TdaDicke,self).__init__(energieli,ontaardingi,gi,tdadic = tdadic)
    self.tdasolutions = np.zeros(self.nlevel+1,np.float)
 
  def tda_vgl(self,x):
    """
    returns the value of the tda eq. defined by 
    the single particle levels("energiel" and "epsilon0"), degeneracy("ontaarding"), interaction constant("g") and x value.
    The solutions of this equation are the values x that sets the equation underneath to zero.
    """
    opl = 1./(self.g)*(self.epsilon0-x)-self.g/2.*sum(self.ontaardingen/(self.energiel-x))
    return opl
  
  def calc_singularity(self ):
    return self.energiel 
  
  def get_lowerboundaries(self):
    epsilon_1 = self.energiel[0]
    return 0.5*(-epsilon_1+self.epsilon0) - 0.5*np.sqrt((epsilon_1-self.epsilon0)**2+2*self.g**2*sum(self.ontaardingen))
   
  def get_upperboundaries(self):
     epsilon_m = self.energiel[-1]
     return 0.5*(-epsilon_m+self.epsilon0) + 0.5*np.sqrt((epsilon_m-self.epsilon0)**2 + 2*self.g**2*sum(self.ontaardingen))

  def calc_coef(self, tdasol):
    return  1./(self.energiel - tdasol)  
  
  def __str__(self):
    s = super(TdaDicke,self).__str__()
    s += '\n#En epsilon_0 wordt gegeven door: %f' %self.epsilon0 
    return s

  def bisect_tda(self,tdadict = None):
    """
    bisection method to find the zeros of the Dicke tda equation
    """
    nlevel = len(self.energiel)
    if tdadict == None:
      tdadict = self.tdadict
    self.tdadict = tdadict #to make sure that the tdadict the tdasolve object has is consistent with the tdasolutions vector
    nauwkeurigheid = 1e-14
    verschil =1.
    maxnumbert = 1000 #max number of tries the bisection method searches for a solution in a given interval
    tel = 0
    for key in tdadict.keys():   
      if key == 0 and tdadict[0] is not 0:
      #Lower bound
        x1 = self.singularitys[0] + self.get_lowerboundaries()
        x2 = self.singularitys[0]
      #Extra solution necessary for the Dicke model
      elif key == nlevel and tdadict[key] is not 0:
        x1 = self.singularitys[key-1]
        x2 = self.singularitys[key-1] + self.get_upperboundaries()
      elif tdadict[key] is not 0:
        x1 = self.singularitys[key-1]
        x2 = self.singularitys[key]
      else:  
        continue
      numbert = 0
      while(abs(verschil) > nauwkeurigheid):
        midden = (x1+x2)/2.
        verschil = self.tda_vgl(midden)
        if  verschil > 0.:
          x2 = midden
        else:
          x1 = midden
        numbert += 1
        if numbert > maxnumbert:
          print 'tdasolutions didn\'t converge'
          break
      
      self.tdasolutions[key] = midden
      verschil = 1.
    print 'bisect heeft volgende tda_oplossingen gevonden:%s , bij volgende tdadict: %s' %(str(self.tdasolutions),str(tdadict)) 
    return self.tdasolutions
  


  def plot_tdaopl(self, solutions = True):
    """
    plot the graphical solution method of the tda equation
    in the Dicke model and the corresponding solutions (if solutions = true). 
    """
    fig1=pl.figure(1)
    #plot then the summation term of the tda equation relevant for the problem in consideration (xlim and ylim good so we see usefull stuf)
    einde = self.singularitys[-1] + self.get_upperboundaries()
    begin = self.singularitys[0] + self.get_lowerboundaries()
    xw = np.arange(begin,einde,0.01)
    #plot the linear term in the tda equation
    plot1 = pl.plot(xw,[(self.epsilon0-x)/self.g**2 for x in xw],'r',label = 'lineaire term')
    #plot the summation term
    fxw = [(-self.tda_vgl(x)+1./self.g*(self.epsilon0-x))/self.g for x in xw]
    plot2 = pl.plot(xw,fxw,'b',label = 'sommatieterm TDA sec. vgl.')
    pl.xlim(xmax=einde,xmin=begin)
    pl.ylim(ymax = 1./(self.g)**2*(self.epsilon0-begin) , ymin = 1./(self.g)**2*(self.epsilon0-einde))
    pl.xlabel(r'$E_{\alpha}$ (a.e.)')
    pl.ylabel('de twee delen van de TDA sec. vgl.')
    for i in xrange(self.nlevel+1):    
      if solutions is True:
        if i == 0:
          plot3 = pl.plot(self.tdasolutions[i],1./self.g**2*(self.epsilon0-self.tdasolutions[i]),'g.',markersize=20.0,label ='oplossingen' )
        else:
          pl.plot(self.tdasolutions[i],1./self.g**2*(self.epsilon0-self.tdasolutions[i]),'g.',markersize=20.0 )
      if i != self.nlevel:
        pl.axvline(x = self.singularitys[i] ,c=  'k')
        pl.text(self.singularitys[i], self.g*0.5 ,"o:%g"%self.ontaardingen[i])
    
    pl.legend(loc = 'lower left', numpoints=1)
    pl.savefig(str(self.g)+'-TDA.png')
    pl.show()

    
def user_input():
  """
  asks the user all the necessary variables to determine a particular tda equation
  """
  aantalel = input("Geef het aantal eendeeltjesniveaus")
  #g is negatief voor aantrekkende interactie en positief voor afstotende interactie
  g = input("Geef de koppelingsconstante")
  eendlev = np.zeros(aantalel,np.float)
  ontaardingen = np.zeros(aantalel,np.float)
  x = np.zeros(aantalel,np.float)
  #voer energieniveaus van klein naar groot in nog niet aangepast dat invoer in willekeurige volgorde mogelijk is
  for i in xrange(aantalel):
    a = input("Geef het eendeeltjesniveau?")
    b = input("Geef de bijbehorende ontaarding")
    eendlev[i] = a
    ontaardingen[i] = b
  return g,eendlev,ontaardingen
  
"""
underneath are some tests of the tda equation solvers and newton-raphson method
"""
def main():
  """
  Some tests of the tda solver that is used in the RG code.
  """
  
  #g,eendlev,ontaardingen = user_input()
  eendlev = np.arange(1,13)
  ontaardingen = np.ones(len(eendlev))*2.
  g = -1.
  '''
  eta = 1.
  tdad = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1}
  tda = TdaRedBcs(eendlev,ontaardingen,g)
  '''
  g = -1.
  epsilon0 = 5.5
  tdad = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1}
  tda=TdaDicke(eendlev,ontaardingen,g,epsilon0)
  tda.bisect_tda(tdad)
  print str(tda)
  tda.plot_tdaopl()
  
if __name__ == "__main__":
  main()
