import numpy as np
from numpy import ones, zeros ,array, sort,linalg,empty
import random,sys

import richardsongaudin as rg
import writepairing as wp
import rgfunctions as rf
import plotfunctions as pf

class Continuity(object): #needs to be efficiently written because it will be one of the most used functions
  def continuitycheck(self,rgeq,carray,step,end, ccrit = 2.,pf = None):
    '''
    Increases Xi and find the new rg solutions but assure continuity
    REMARK: you can see this function as a wrapper function of both nr.solve and rgf.continuity_check
    '''
    consol = False
    carray2 = list(carray)
    n = 1
    i = 1
    newstart = self.getvar(rgeq) ; savestep = step
    saverg = rgeq.rgsolutions
    while(consol == False):
      for a in range(n):
	self.changevar(rgeq,step)
	consol = self.endcheck(rgeq,step,end,carray)
	rgsol = rgeq.solve()	#MOST IMPORTANT STATEMENT of this function
      if not consol:
	try:	
	  self.continuity_check(rgeq,step,carray,rgsol, crit = ccrit)
	except rf.DiscontinuityError as e:
	  rgeq.rgsolutions = saverg
	  if pf is not None:
	    pf.write("#We arrived at a discontinuity\n")
	  print 'Problems with continuity of the energy with chance in main_rgsolver()' ,e
	  self.setvar(rgeq,newstart)
	  savestep /= 10.**i
	  n = 10**i
	  i += 1
	  if n > 2000:
	    print 'error we have discoverd a super big discontinuity in the energy when we variate xi from 0 to 1'
	    raise ValueError #so it will handled by making g complex
	else:
	  consol = True
    return rgsol 
    
  def continuity_check(self,rgeq,step,arraysol, nval,crit = 1.5):
    '''
    function that investigates, if the new solution: $nval forfilles the continuity criterium:
    that is that the mean difference of the previous $grootte solutions *$crit is bigger then the difference between $nval
    and the previous solution. If this is the case the last element of $arraysol is deleted and $nval is inserted at the end.
    If this is not the case $arraysol is not changed and the variable $con is set False
    n is a value that says how many steps further the function went
    REMARK: Very important function for (generating_datak) to check if there is no jump to an excited state.
    '''
    grootte = 5
    meand = 0
    if len(arraysol) < grootte:
      arraysol.append((self.getvar(rgeq),rgeq.rgsolutions,nval))
    else:
      adapter = abs((self.getvar(rgeq) - arraysol[-1][0]) /(arraysol[-1][0]-arraysol[-2][0])) #we use this to take into account largerstepsizes
      for i in range(grootte-1): 
	meand += abs(arraysol[i][2]- arraysol[i+1][2])
      verschil = abs(arraysol[-1][2] - nval)      
      if verschil < meand/(grootte-1)* crit*adapter:
	del(arraysol[0])
	arraysol.append((self.getvar(rgeq),rgeq.rgsolutions,nval))
      else:
	raise rf.DiscontinuityError(step,crit,arraysol,nval)
  
  def endcheck(self,rgeq,step,end,carray):
    bool = False
    if self.getvar(rgeq) +abs(step)* 0.8 >= end and self.getvar(rgeq) - abs(step) * 0.8 <= end:
      self.setvar(rgeq,end)
      bool = True
    return bool
  
  def get_back(self,rgeq,carray, i = -1):
    if carray is not []:
      self.setvar(rgeq,carray[i][0])
      rgeq.rgsolutions = carray[i][1]
    return rgeq
  
class Xi_Continuity(Continuity):    
  def changevar(self , rgeq,var):
    rgeq.xi += var
  def setvar(self,rgeq,var):
    rgeq.xi = var
  def getvar(self,rgeq):
    return  rgeq.xi
       
class RG_Continuity(Continuity):
  def changevar(self,rgeq,var):
    rgeq.g  += var
  def setvar(self,rgeq,var):
    if rgeq.g.imag != 0.:
      rgeq.g = var + 1.j*rgeq.g.imag
    else:
      rgeq.g = var
  def getvar(self,rgeq):
    return rgeq.g.real
    
class CG_Continuity(Continuity):
  def changevar(self,rgeq,var):
    rgeq.g += 1.j* var
    if rgeq.g.imag == 0.:
      rgeq.g = rgeq.g.real
  def setvar(self,rgeq,var):
    rgeq.g= rgeq.g.real+ 1.j*var 
  def getvar(self,rgeq):
    return rgeq.g.imag
   
class ELevel_Continuity(Continuity):
  def changevar(self,rgeq,var):
    rgeq.energiel += var #var can be a float or an array with the same length because self.richeq.energiel is an numpy array
  def setvar(self,rgeq,var):
    rgeq.energiel = var
  def getvar(self,rgeq):
    return rgeq.energiel
  def endcheck(self,rgeq,step,end):
    pass

class Change_G(object):
  '''
  This class has as main task the circumvention of a critical point of a Richardson-Gaudin look alike set of equations
  '''
  def __init__(self,rgeq,cg = True):
    self.richeq = rgeq.copy()
    if cg == True: self.concheck = CG_Continuity()
    else: self.concheck = RG_Continuity()
    assert(self.richeq.rgsolutions != None)
    
  def change_g(self,changeg,n = 10.):
    print 'we entered change_g of the class Change_G'
    assert(self.richeq.rgsolutions is not None)
    conarray = []
    gstep = changeg['gstep'] ; gend =  changeg['gend']
    if isinstance(self.concheck,CG_Continuity): 
      self.richeq.g = self.richeq.g+ 1.j* abs(self.richeq.g)/(10000.*n)
      self.richeq.solve()
    while(self.concheck.getvar(self.richeq) != gend or self.richeq.g.imag == self.richeq.g.real/10.):
      print self.richeq.g 
      self.concheck.continuitycheck(self.richeq,conarray,gstep/n,gend)
    print 'we changed g: %s' % str(self.richeq.g)
    return self.richeq
  
  
  def circumvent_point(self,critr,stepr,critc,stepc,n):
    cchange = {'gstep':stepc,'gend':critc}
    rchange = {'gstep':stepr,'gend':critr}
    self.change_g(cchange,n)
    self.change_g(rchange,n)
    cchange = {'gstep':stepc/2.,'gend':0.}
    self.change_g(cchange,n)
    return self.richeq
    
class Evolver(object):
  '''
  This class generates some parallel RG equations, (with different xi , complex g , different xi and complex g)
  '''
  def __init__(self,rgeq,nxi = 3 , ncg = 3):
    self.nc = nxi+2*ncg
    self.gsolutions = []
    self.richeq = rgeq
    if self.richeq.rgsolutions is None:
      self.get_firstsol()
    self.gsolutions.append((self.richeq.g,self.richeq.get_energy(),self.richeq.rgsolutions)) 
    self.xilist = self.get_xiarray(nxi)
    self.compglist = self.get_compglist(ncg, n = 10.)
    self.mixlist = self.get_mixlist()
    self.lists = {'xi':self.xilist , 'cg' : self.compglist , 'mix': self.mixlist}
    
  def get_xiarray(self,nval):
    print 'we start with the creation of the xiarray'
    xilistd = []
    solver = rg.RichardsonSolver(self.richeq)
    for i in reversed(np.linspace(0.1,0.90,nval)):
      xiarg = {'xistep' : -0.01, 'xiend' : i}
      try:
        rgeq = solver.change_xi(xiarg)
      except (ValueError , np.linalg.linalg.LinAlgError) as e:
	print 'error in making xilist'
      xilistd.append(rgeq.copy())
    print 'creating xi list succeeded' 
    return xilistd
    
  def get_compglist(self,nval,complexstep = 20000., n =10.):
    glist = []
    gchanger = Change_G(self.richeq.copy())
    for j in [abs(self.richeq.g/100.) * n**i for i in range(nval)]:    
      cchange = {'gstep': (abs(self.richeq.g)+abs(j)*10.)/complexstep  , 'gend': j }
      try:
        rgeq = gchanger.change_g(cchange,10.)
      except (ValueError , np.linalg.linalg.LinAlgError) as e:
	print 'error in making compglist'
      glist.append(rgeq.copy())
    print 'we have the following glist' 
    return glist
  
  def get_mixlist(self):
    a = np.linspace(0.2,1.,len(self.compglist))
    random.shuffle(a)
    mlist = []
    for i in range(len(self.compglist)):
      xiarg = {'xistep': -0.01 , 'xiend' : a[i]}
      try:
        rgeq = rg.RichardsonSolver(self.compglist[i]).change_xi(xiarg)
      except (ValueError , np.linalg.linalg.LinAlgError) as e:
	pass
      mlist.append(rgeq)
    print'we have the following mixlist' 
    return mlist
    
  def evolve(self,stepg,endg):
    while(self.richeq.g is not endg):
      self.richeq = self.richeq_step(stepg)
      if self.richeq.g.real + abs(stepg)/2. >= endg and self.richeq.g.real - abs(stepg)/2. <= endg:
        self.richeq.g = endg
	self.richeq.solve() ; break
      print self.richeq.g
  
  def richeq_step(self,stepg= None):
    if stepg is not None: 
      self.richeq.g += stepg 
    
    foundsol = False
    try:
      self.richeq.solve()
    except (ValueError , np.linalg.linalg.LinAlgError) as e:
      i = 0
      for klist in sorted(self.lists.keys()): #make sure we probe first the list with complex g , then the mixlist and finally the xilist
        for rgeq in self.lists[klist]: 
          i+= 1
          if rgeq.rgsolutions is not None:
	    changer = {'gstep': stepg , 'gend' : self.richeq.g+stepg}
	    changec = {'gstep': -1.*rgeq.g.imag/300. , 'gend': 0.}
	    changev = {'xistep': 0.01 , 'xiend': 1.}
	    print str(rgeq) ,  changer , changec , changev
	    try:
	      rgeq = Change_G(rgeq, cg= False).change_g(changer)
	      rgeq2 = rg.RichardsonSolver(rgeq).change_xi(changev)
              self.richeq = Change_G(rgeq2).change_g(changer,20)
              foundsol = True ;break
            except:
	      #if i == self.nc: print 'BIG ERROR choose other or more circumventers ' ;sys.exit(1)
	      print 'warning none of the start rg equations can reach the demanded solution'
	      pass
	   
	if foundsol == True: break    
    self.gsolutions.append((self.richeq.g,self.richeq.get_energy(),self.richeq.rgsolutions)) 	  
    return self.richeq
  
  def restore(self,stepg,n):
    for list in self.lists.values():
      for rgeq in list:
	if rgeq.rgsolutions == None:
	  changer = {'gstep': stepg , 'gend': rgeq.g.real}
	  changev = {'xistep': 0.01, 'xiend': rgeq.xi}
	  changec = {'gstep': self.richeq.g/10000 , 'gend': rgeq.g.imag}
	  a = Change_G(rgeq)
	  rgeq = a.change_g(changer)
	  rgeq = a.change_g(changec,20)
	  rgeq = rg.RichardsonSolver(self.richeq).change_xi(changev) 
	
  def delete_none(self):
    for list in self.lists.values():
      for rgeq in list:
	if rgeq.rgsolutions == None:
	  list.remove(rgeq) 
    
    
  def get_rg(self, lis = 'xi',nval = 0):
    return self.lists[lis][nval] 
    
  def get_firstsol(self,pairingdict = None, kop = None):
    d = rf.calculated(self.richeq.energiel,self.richeq.ontaardingen)
    self.richeq = rf.genstartsol(self.richeq,d,self.richeq.g,pairingd = pairingdict , begin = kop)
  
  def write_gpath(self,pairingd,exname = '',rgwrite = True):
    fwrite = open("plotenergy%s.dat" %exname,'w')
    plotenergyfile = open("plotenergy%s.dat" %exname,"w")
    rf.info_1set(plotenergyfile,str(self.richeq) , exinfo = "#g\tcE\tgE\trgvar(real)\trgvar(imag)\t ...\n" ,tdadict = pairingd)
    bb = rf.calcnintgrondtoestand(self.richeq.apair, self.richeq.energiel, self.richeq.ontaardingen, self.richeq.senioriteit)
    a = sorted(self.gsolutions, key = lambda opl : opl[0].real)
    for gopl in a:
      plotenergyfile.write("%f\t%f\t%f" %(gopl[0].real, gopl[1]-bb, gopl[1])) 
      if rgwrite is True:
	for i in range(self.richeq.apair):
	  plotenergyfile.write('\t%f' %( gopl[2][i].real ))
	  plotenergyfile.write('\t%f' %( gopl[2][i].imag ))
      plotenergyfile.write('\n')
    plotenergyfile.close()

def main():
  el = np.arange(1,13)
  deg = ones(12,float) *2.
  sen = zeros(12,float)
  apair = 6
  pd = {0:6}
  g= -1.
  eta = 1.
  rgeq = rg.RichRedBcs(el,deg,sen,g,apair)
  endg = -0.001
  stepg = +0.003
  e,rgeq = rg.RichardsonSolver(rgeq).main_solve(pd)
  tdad = rg.RichardsonSolver(rgeq).main_desolve(plotrgvarpath = True)
  evolution = Evolver(rgeq,nxi = 3, ncg= 3)
  evolution.evolve(stepg,endg)
  evolution.write_gpath(pd)
  pf.plotrgvars(apair,ref = 'plotenergy.dat',afhvar = 'g (a.u.)',namerg = 'rgvar',stop = None,begin = 0,istart = 3)
  
if __name__ == '__main__':
  main()
  
  
    