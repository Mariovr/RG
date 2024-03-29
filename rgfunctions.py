# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Mario Van Raemdonck, 2013;
# (c) Ghent University, 2013
import sys,math , os  , shutil 
import numpy as np
from numpy import ones, zeros ,array, linalg,empty
from itertools import combinations,combinations_with_replacement
from operator import mul
import linecache,copy

import richardsongaudin as rg
import writepairing as wp

#to surf with sockets commandline: ssh -D 32190 molgate0 (32190 is a random number lower then 64000)
#proxy instellen -> manual , sockets : 127.0.0.1 (localhost) vakje er naast(random number dat je in de commandline hebt getypt)

def littleLoop(rgeq,stepg,n,complexstepd = 10000,end = None,backxi = False,xival = 1. , pairingdtje = None ,dvar = 'g'):
  '''
  function makes the interaction constant complex in order to circumvent a critical point
  it makes a little loop around the critical point after the execution we have a rg solution one stepg further
  in realspace of g if the critical point is not yet passed we need to execute this function again until it's passed
  (this function is called when their was an error that signals that their is a critical point in the neighborhood)
  REMARK: g is always negative
  REMARK: you can use this functions on two different ways: 1) to circumvent a critical point at some g then it just increases 
  netto with one step of the interaction constant of the wrapper function that calls this function.
  2) to go from a start interactionconstant 
  '''
  rgeq = copy.deepcopy(rgeq)
  if rgeq.rgsolutions is None:
    rgeq = rg.RichardsonSolver(rgeq).main_solve(pairingdtje, xival = xival)
  print ('######################################')
  print ('Making g complex to circumvent a critical point on rgeq = %s and n = %g, stepg= %s, complexstepd = %s , dvar = %s, end = %s' %(str(rgeq),n,str(stepg),str(complexstepd) ,str(dvar),str(end) ))
  print ('######################################')
  phi = abs(rgeq.g)/(1000.*n)
  if end == None:
    crit = abs(rgeq.g)/100.
    if n > 10:
      crit *= 10.
    phis = abs(rgeq.g)/complexstepd #step of phi the complex variable
  else:
    crit = abs((rgeq.g+end)/10.)/10.
    if n > 10.:
      crit *= 10.
    if n > 1000:
      crit *= 10.
    if n > 10000:
      crit /= 100.
    n = int(abs(end-rgeq.g)/abs(stepg)) + 1
    phis = abs(rgeq.g+end)/(complexstepd*2.) #step of phi the complex variable
  assert(rgeq.rgsolutions is not None)
  rgeq.g += 1.j * phi
  while(rgeq.g.imag < crit):
    #little block to become a better guess for the rgvars
    energierg = rgeq.solve()
    rgeq.g += 1.j*phis
    #print ('distance measure rgeq.g = : %s' % str(rgeq.g))
  print ('we made g complex now we are going to increase the norm of g: %s' %str(rgeq.g))
  argxi = {'xiend': xival}
  if backxi is True:
    argxi = {'xiend': 0.6, 'xistep': -0.01}
    print ('we start by reducing the xi value to xi = %f' %argxi['xiend'])
    rgeq = rg.RichardsonSolver(rgeq).change_xi(argxi)
  nauw= 100. # the bigger the better the accuracy
  n = n *nauw
  stepg /= nauw
  saverg1 = None ; saverg2 = None
  for i in xrange(int(n)):
    if saverg1 is  None:
      saverg1 = rgeq.rgsolutions
    if saverg1 is not None and saverg2 is None :
      saverg2 = rgeq.rgsolutions
    if saverg1 is not None and saverg2 is not None:
      saverg1 = saverg2
      saverg2 = rgeq.rgsolutions
      rgeq.rgsolutions += saverg2-saverg1
    rgeq.setvar(dvar,rgeq.getvar(dvar) + stepg)
    if end is not None:
      if (rgeq.getvar(dvar).real < end and stepg < 0) or (rgeq.getvar(dvar).real > end and stepg > 0):
        rgeq.setvar(dvar,end+ 1.j *  rgeq.getvar(dvar).imag)   
    energierg = rgeq.solve() 
    if end is not None and rgeq.g.real == end:
      break
    #print 'The variable we change in littleLoop : %s is: %s ' %(dvar ,str(rgeq.getvar(dvar)))
  print ('Now are we going to make g back real %s' %str(rgeq.g))
  phis /= 2.
  rgeq2 = copy.deepcopy(rgeq)
  if backxi is True:
    argxi = {'xiend': xival, 'xistep': 0.01}
    print ('now we go back to xi = xival')
    rgeq = rg.RichardsonSolver(rgeq).change_xi(argxi)
  while(rgeq.g.imag >  0.):
    rgeq.g -= 1.j * phis
    if rgeq.g.imag <= 0.:
      rgeq.g = rgeq.g.real
    energierg = rgeq.solve()
    #print ('distance measure rgeq.g = : %s' % str(rgeq.g))
  print ('circumventing the critcal point on %f succeeded !!!!!!!!!!!!!!!!!' % rgeq.g.real)
  return energierg,rgeq,rgeq2
  
def genstartsol(rgeq,d,end,pairingd = None , begin = None):
  '''
  Generates a start solution, if no solution is known in one of the two limits where it's easy to find a solution
  REMARK: make sure that step is positive
  TODO: couple bcs code to this code because my thesis indicated that the critical point at the smallest interaction
  constant was correlated to the point where the BCS gap becomes zero.
  #################
  RETURNS: float(energierg) the energie of the solution of the Richardson-Gaudin variables, np.array(apair,np.complex)  the Richardson-Gaudin
  variables at the corresponding solution, float(begin) the interaction constant of the corresponding solution
  '''
  rgeq.rgsolutions = None ; rgeq.g = begin ;rgeq.xi = 0.
  if isinstance(rgeq, rg.RichRedBcs):
    multiple = 1.
  elif isinstance(rgeq,rg.RichFacInt):
    multiple = 10. * d
  if begin is None: #search for the optimal start g
    if d/(2.*multiple ) * 10. > abs(end):
      if pairingd is None:
        pairingd = tdadict_kleinekoppeling(rgeq.apair,rgeq.ontaardingen,rgeq.senioriteit)
      if end > 0 : #takes into account positive interaction constants and makes sure we increase the interaction constant
        rgeq.g = 0.001/multiple
      else: 
        rgeq.g = -0.001/multiple
      assert(abs(rgeq.g) < abs(end))
    else:
      if pairingd is None:
        pairingd = {}
        pairingd[0] = rgeq.apair
      if end > 0: rgeq.g = abs(math.floor(d)*1.5) 
      else: rgeq.g = -abs(math.floor(d)*1.5 )    
      assert(abs(rgeq.g) > abs(end)) # otherwise it's useless to go to smaller ic (ic > 0)
  print 'in genstartsol the RG equations are now %s' %str(rgeq)
  while(rgeq.xi != 1.0):
    try:
      print "searching for first solution"
      rgeq = rg.RichardsonSolver(rgeq).main_solve(pairingd)
      print 'We found a good startsolution'
    except (ValueError,np.linalg.linalg.LinAlgError) as e:
      if abs(rgeq.g) < abs(end):
        rgeq.g /= 10.
        print "the start interaction constant was not small enough so we reduced it to %f" % rgeq.g
        if abs(rgeq.g) < 1e-5:
          print 'BIG ERROR we dont find a start solution'
          return None , None #use this when using generating_data
      else:
        rgeq.g *= 2.
        print "the start interaction constant was not strong enough so we increased it to %f" % rgeq.g
        if abs(rgeq.g) > 20*abs(d):
          print 'BIG ERROR we didn\'t find a start solution'
          return None , None #use this when using generating_data
  return rgeq
  
def generate_dir(name,infilename,args = None, prefix = ''):
  """ generates a seperate dir for all the information of the
  handled problem
  """  
  dirname = os.path.join('results', "%s" %name)
  print ('#making directory for the results: %s' %dirname)
  dir = dirname #name of the directory
  #check if the directory already exists
  if os.path.isdir(dir): #OEPS the directory exists already
    no = 'l'
    while  no != 'y' and no != 'n':
      no = raw_input('remove dir %s: (y,n)' %name).lower()
    if no == 'y':
      shutil.rmtree(dir)
    else:
      print 'the directory you wanted to create does already exist so we added 1 to the end'
      dir += '1'     
      
  os.makedirs(dir)
  if isinstance(infilename,str): 
    shutil.copy(infilename,dir)
  if args is not None:
    for i in args:
      shutil.copy(i,dir)
  if prefix is not '':
    x = [i for i in os.listdir(os.getcwd()) if prefix in i]
    for a in x:
      shutil.copy(a,dir)
  os.chdir(dir)

  
def calcnintgrondtoestand(rgeq):
  """
  calculates the non-interacting ground energy of the problem by putting the pairs in the lowest sp levels
  VARIABLES: number of pairs, list(or array) of sp levels, list(or array) of the degeneracies
  REMARK: with seniority zero for all sp levels
  """
  som,i = 0  ,0
  #Moet functie bij voor het Dicke model
  apair = rgeq.apair
  while apair > 0:
    som += rgeq.energiel[i]*(rgeq.ontaardingen[i]-rgeq.senioriteit[i]*2)
    apair -= (rgeq.ontaardingen[i]/2-rgeq.senioriteit[i])
    i += 1
  if apair < 0:
    som += rgeq.energiel[i-1]*apair*2
  return som  

#functie die de gemiddelde afstand tussen de eerste alev niveaus, in waarden meegeeft
def calculated(waarden,deg):
  """
  calculates the mean distance between the sp energy levels,
  at this moment without taking into account degeneracies of particular levels
  """
  som = 0.
  for i in xrange(len(waarden)-1):
    som += (waarden[i+1]-waarden[i])
  som /= sum(deg)/2
  return som  


def windowcreation(energielev, deg,sen,npair,hwinu, hwinb):
  '''
  INPUT
  ######
  list (energielev): containing the energielevels, list(deg): containing the degeneracies, list(sen): containing the seniority's ,integer(npair): the total amount of pairs,
  float(hwinu): the parameter that defines the window under the fermilevel, float(hwinb): the parameter that defines the window above the fermilevel
  ########
  REMARKS:
  1) It's a function that reduces the energylevels such that only the energeylevels within an interval [Fermi-hwinu,hwinb+Fermi]
  are kept (REMARK: make sure that hwin is kept positive)
  2) It also calculates the new amount of particles now contained in the interval.
  #######
  OUTPUT:
  list(energielw): the energielevels in the window
  list(deglw): the corresponding degeneracies
  list(senlw): the corresponding seniority's
  integer(npairw): the number of pairs in the window
  integer(nlevelw): the number of levels in the window
  float(extrae): the energycontribution of the electrons in the non interacting sp levels
  '''
  print ('#starting to create a window for the sp levels that interact with each other as in the reduced BCS Hamiltonian')
  fermilev ,pairf = get_fermilevel(npair,energielev,deg)
  print ('#the fermilevel is: %f' %fermilev)
  i = 0
  pairout = 0
  while(fermilev - hwinu > energielev[i]):
    pairout += deg[i]/2
    i += 1
  npairw = npair - pairout 
  print ('npair before :%g ,npair after: %g'  %(npair,npairw ))
  beginelw = i
  i = len(energielev)-1
  while(fermilev + hwinb < energielev[i]):
    i -= 1 
  endelw = i+1
  energielw = energielev[beginelw:endelw]
  deglw = deg[beginelw:endelw]
  senlw = sen[beginelw:endelw]
 
  nlevelw = len(energielw)
  #the calculation of the energycontribution of the electrons in the non interacting sp levels
  ind = energielev.index(energielw[0])
  extrae = calcnintgrondtoestand(pairout,energielev[:ind], deg[:ind],sen[:ind])
  print ('#the energy of the electrons in the not interacting sp levels equals %f'%extrae) 
  return energielw, deglw, senlw, npairw, nlevelw,fermilev , extrae

def get_fermilevel(npair,energielev,deg):
  '''
  function that determines the fermilevel given the amount of pairs, the set of energielevels, the degeneracies of each level
  #######
  OUTPUT:
  a tuple of length 2 : the first argument is the fermilevel, the second argument is the number of pairs in the fermilevel 
  (if this amount is None , the fermilevel is in between 2 sp levels and can't contain any pairs)
  '''
  i = 0
  while(npair > 0):
    npair -= deg[i]/2
    i += 1
   
  if npair < 0:
    fermil = energielev[i-1]
    pairf = deg[i]/2 + npair
  elif npair == 0:
    fermil = (energielev[i-1] + energielev[i])/2.
    pairf = None  
  return fermil,pairf 
  
def get_npair(deg,occ):
  npair = 0
  for i in range(len(occ)):
    npair += deg[i]/2
  return npair
  
def wd_processing(wd2,nlev,npair,waarden):
  """
  function that keeps only the sp levels E_n such that E0 + WD2 > En
  nlev = number of energylevels , npair = number of pairs, waarden = energylevels, wd2 = upper boundary energylevels
  #REMARK this wd2 definition is not physical we need to write another function that defines wd2 around
  the Fermilevel (not only the amount of levels in the interaction band changes but also the amount of pairs need to change then
  as the dependent variabele that characterizes the sp levels changes)
  """
  print 'start wd2 processing'
  nlevel = nlev
  if wd2 > 0:
    for i in xrange (0,len(waarden)-1):
      if waarden[i]- waarden[0] > wd2:
        nlevel = i
        if nlevel < npair:
          print "het aantal energieniveaus bij 2wd = %f : is nlevel = %f dit is kleiner als apair %f dus fatale fout (oplossing: 2wd vergroten)" %(wd2,nlevel,npair)
          sys.exit(1)
        break
  waarden = waarden[0:nlevel]
  return nlevel,waarden

def checkdegeneratie(waarden,deg,nauw = 0.00000000001):
  """
  It determines the degeneracies of an sp spectrum en puts two levels as one if they are degenerate or quasi degenerate (closer then nauw together)
  And then it increases the degeneracie of the level
  the function asks a list of sp levels , and the number of sp levels
  it returns the adjusted number of splevels and an array that contains the degeneracy of the sp levels (the sp level list is adjusted by reference)
  """
  i = 0 ;  nlev = len(waarden)
  while( i < nlev):
    a = 0 # number of pairs a sp level can contain
    c = i
    degsom = 0
    while (waarden[c]- nauw <= waarden[i]):
      a = a + 1
      degsom += deg[c]
      c += 1
      if c >=nlev:
        break
    deg[i] = degsom
    for m in range(a-1):
      del(waarden[i])
      del(deg[i+1]) #because otherwise the just calculated degeneracy of the level is deleted
    nlev -= a-1
    assert(nlev == len(waarden))
    print i,nlev,degsom
    i = i +1 
  assert(nlev == len(deg))
  return nlev, deg

def uncheckdegeneratie(waarden,deg):
  """
  undo checkdegeneratie
  REMARK: returns the new list of energylevels and degeneracies -> so no in place change of the energylevels
  """
  newe = []
  for i , de in enumerate(deg):
    num = int(de/2)
    newe += [waarden[i]] * num
  deg = ones(len(newe)) *2.
  sen = zeros(len(newe)) 
  return array(newe) , deg, sen , len(newe)

def tdadict_kleinekoppeling(npair,degen,sen):
  """
  function that creates the dictionary that relates the pairs to the tda solutions if the interaction constant is really weak.
  ASKS: dictionary, number of pairs in the system , list(or array) of degeneracies
  CHANGES the dictionary in the function returns nothing
  """
  pairingd = {}
  i,up = 0,0
  while up < npair:
    a = int(degen[i]/2)-sen[i]
    pairingd[i] = a
    i += 1
    up += a
  if up > npair:
    pairingd[i-1] -= (up-npair)
  return pairingd
  
def mergefiles(file1,file2,reversed = False):
  f = open(file1,'r')
  g = open(file2,'r')

  mergen = file1+'_'+file2

  m = open(mergen,'w')

  for line in f:
    m.write(line)
  if reversed == True:
    for line in reversed(g.readlines()):
      m.write(line)
  else:
    for line in g:
      m.write(line)
  f.close()
  g.close()
  m.close()
  
#function that generates the degeneracies used for superconductivity (every sp level is 2 times degenarte spin-up and down)
def degeneracies_super(nlevel):
  return ones(nlevel,float)*2.   

def seniority_enhancer_allstates(rgeq,fname,vsen,exewaarde = 0,begin = -0.0005,step = -0.001,exinfo ='#The paired spectrum of the neutrons of Sn' ):
  '''
  INPUT: list(energielev) the sp energielevels, list(deg) a list with the same length as energielev, that contains the degeneracies
  of the corresponding levels, list(sen) the seniority of the corresponding levels, integer(npair) the number of pairs, integer(vsen) the total 
  seniority that the system needs to have
  ########
  REMARKS:
  1)This function calculates all the eigen-energy's and RG-vars from a system at a given seniority
  2)Watch out only pairs that are already contained in the interaction window are considered, so there's no addition 
  of extra particles from non-interacting states. And also no escape of particles to free states above the window of pairing interaction
  3)It creates a file with the total seniority in the filename, that contains all the excited states of the system.
  4)Only works when the seniority of all  the states is initially zero, it then probes all possible seniority distribution
  of (vsen) unpaired electrons
  5)if a level doesn't interact with the other levels through the pairing interaction anymore because the seniority == degeneration/2 of that sp level
  we delete that level out of the calculation
  6) if this functions is enough tested so it runs wihtout FAULTS you can remove the unnecessary wind construct because it's only implemented
  for when a simulation terminated faulty to restart it where it terminated
  #######
  OUTPUT: nothing
  '''
  rgeq.apair -= vsen/2  #the netto amount of pairs is reduced because of the pairbreaking of vsen/2 pairs
  fd = open(fname + 'sen=%g.dat' %vsen , 'w')
  from datareader import info_1set
  info_1set(fd,str(rgeq) , exinfo = '%s total seniority = %f' %(exinfo,vsen) ,contenergy = exewaarde)
  sencombarray = []
  for i in range(rgeq.alevel):
    sencombarray += [i]*(rgeq.ontaardingen[i]/2)
  sen_places_gen = combinations(sencombarray,vsen)
  #Because if you select zero out of a list you have None, we need to take care of seniority zero states seperately.
  if vsen == 0:
    enil = []
    ontaarding = 1
    fd.write('#We are using seniority %s\n' %str(rgeq.senioriteit))
    allstateslowg(rgeq,fd,ontaarding,extrae = enil,exe = exewaarde)
    #allstatesstrongg(rgeq,fd,ontaarding,rgeq.alevel,extrae = enil,exe = exewaarde)
  else: 
    wind = 0 #variable that determines from which states we are gonna start to look for the solutions, because the order
    sw = 0
    #of the instances in sen_places_gen and tdacombinations is always the same
    deg = list(rgeq.ontaardingen) ; energielev = list(rgeq.energiel) 
    for se_pla in sen_places_gen:
      rgeq.ontaardingen = deg ; rgeq.energiel = energielev 
      rgeq.senioriteit = [0]*len(rgeq.ontaardingen)
      ontaarding = 1
      enil = []
      possol = True
      for i in se_pla:
        rgeq.senioriteit[i] += 1
        enil.append(rgeq.energiel[i])
      fd.write('#We are using seniority %s\n' %str(rgeq.senioriteit))
      print ('the seniority is %s' %str(rgeq.senioriteit))
      j = 0
      while j < len(rgeq.senioriteit):
        assert(rgeq.senioriteit[j] <= rgeq.ontaardingen[j]/2.)
        if rgeq.senioriteit[j] == rgeq.ontaardingen[j]/2.:
          del(rgeq.energiel[j]) ; del(rgeq.ontaardingen[j]) ; del(rgeq.senioriteit[j]) #the level is full with unpaired electrons so the level doesn't take part in the pairing interaction
          j -= 1
        j+= 1
      #calculate the degeneracy of the level  
      for i in range(len(rgeq.senioriteit)):
        if rgeq.senioriteit[i] is not 0:        
          ontaarding *= binoml(rgeq.ontaardingen[i]/2.,rgeq.senioriteit[i])
      print rgeq.energiel,rgeq.ontaardingen, rgeq.senioriteit,rgeq.apair, enil
      assert(len(enil) == vsen) #check if we take into account all unpaired elektrons in total energy later  
      if rgeq.apair is not 0:
        wind = allstateslowg(rgeq,fd,ontaarding,extrae = enil,exe = exewaarde,wind = wind,startw = sw,begin = begin,step = step)
      else:
        fd.write("%s\t%f\t%f\t%s\tIm\t%s\n" %(str(dict), sum(enil)+exewaarde,ontaarding,'geen RGvars' ,'geen RGvars'))    
  
  fd.close()
  return 0

  
def allstateslowg(rgeq,fd,ontaarding,extrae = [],step = -0.001,exe = 0,wind = 0,startw = -1,begin = -0.0005):
  '''
  INPUT:Richardson-Gaudinequations(rgeq),filehandler(fd): file where the outputs is going to be written
  list(extrae): list of sp levels in the pairing window where unpaired elektrons are situated
  #############
  Remarks:
  1) runs over the tdastartdistributions that give the eigenstates at small g and calculates the corresponding energy and Richardson-Gaudin variables
  2) And it variates the g then to where 
  3) all the calculated information is written to the file where the filehandler(fd) communicates with
  4) This function is preferable over allstatesstrongg if the interaction constant where we want to know all the states from is close to the 
  small interacting limit.
  j) wind : construction used for debugging and overcoming unovercomable critical points
  ############
  Output: None
  '''
  tdacombinations = combinations_with_replacement(np.arange(0,rgeq.alevel),rgeq.apair)
  afhxas = 'g' ; ende = rgeq.g
  #ende = rgeq.g +step #use this when using generating_datak
  d = calculated(rgeq.energiel, rgeq.ontaardingen)
  for tdadict in tdacombinations:
    #tdastart needs to be a dictionary so we need to convert the list that contains one element of the permutation sequence to a dictionary    
    tdastartd = {}
    goodsol = True
    ontaarding = 1.
    for i in tdadict:
      a = tdadict.count(i)      
      tdastartd[i] = a
      if a*2 + rgeq.senioriteit[i]*2 > rgeq.ontaardingen[i]:
        goodsol = False
    if goodsol == True:
      wind += 1
    if goodsol == True and wind > startw:
      for key,value in tdastartd.iteritems():
        ontaarding *= binoml(int(rgeq.ontaardingen[key]/2.),int(value))
      rgeq.g = begin ; rgeq.rgsolutions = None
      energierg,rgeq = wp.generating_datak(rgeq,tdastartd,afhxas,step,ende , rgwrite = False,exname = '%g' %wind,moviede = False,tdafilebool = False) 
      #littleLoop(rgeq,step,2.,complexstepd = 10000,end = ende,backxi = False,xival = 1.,pairingdtje = tdastartd)
      energierg +=  sum(extrae) +exe  #add the contribution of the unpaired electrons   
      fd.write("%s\t%.12f\t%.12f" %(str(tdastartd), energierg,ontaarding)) 
      for i in range(rgeq.apair):
        fd.write('\t%.12f' %( rgeq.rgsolutions[i].real ))
        fd.write('\t%.12f' %( rgeq.rgsolutions[i].imag ))
      fd.write('\n')
  return wind
  
def allstatesstrongg(rgeq,fd,ontaarding,activelevels,extrae = [],exe= 0 , dataanalyse = {'rgw': False , 'ple' : False , 'plrg' : False , 'ont' : True}): 
  '''
  INPUT:filehandler(fd): file where the outputs is going to be written, list(energielev): list of sp levels, list(degeneration)
  list of degeneracies of the sp levels, float(koppelingsconstante): the interaction constant, int(npair): number of pairs, list(senioriteit) :
  list of the seniority's of the sp levels, array(complex,len(energielev))(rgvar): Richardson-Gaudin variabels to start from (if you have none set on None),
  list(extrae): list of sp levels in the pairing window where unpaired elektrons are situated
  #############
  Remark:
  1) runs over all possible tdastartdistributions and calculates the corresponding energy and Richardson-Gaudin variables
  2) all the calculated information is written to the file where the filehandler(fd) communicates with
  #############
  Output: None
  '''
  #to calculate all the permutations of the pairs in the active sp levels so we choose out a np.arange(nlevel) npair levels where we put our pairs (with repetition)
  rgw = dataanalyse['rgw'] ; ple = dataanalyse['ple'] ; plrg = dataanalyse['plrg'] 
  tdacombinations = combinations_with_replacement(np.arange(activelevels),rgeq.apair)
  for tdadict in tdacombinations:
    dict = {}
  #tdastart needs to be a dictionary; suppose we want to start from a tdadistribution 30 0 0 2 0 3 -> the corresponding tdadict is: {0:30, 3:2, 5:3}
    for a in tdadict:
      dict[a] = tdadict.count(a)
    nosol = False
    print 'we start rg.main_rgsolver with: ', dict#,energielev,koppelingsconstante,npair,sencopy,rgvar,deg
    try:
      rgeqn = rg.RichardsonSolver(rgeq).main_solve(dict,rgwrite = rgw, plotrgvarpath = plrg , plotepath = ple) 
      energierg = rgeqn.solve()+sum(extrae) + exe  #add the contribution of the unpaired electrons
    except rg.XiError as xier:
      nosol = True       
    if nosol == False:
      if dataanalyse['ont'] == True:
        ontaarding = 1.
        for key,value in dict.iteritems():
          ontaarding *= binoml(int(rgeqn.ontaardingen[key]/2),int(value))
      fd.write("%s\t%.12f\t%.12f" %(str(tdastartd), energierg,ontaarding)) 
      for i in range(rgeq.apair):
        fd.write('\t%.12f' %( rgeq.rgsolutions[i].real ))
        fd.write('\t%.12f' %( rgeq.rgsolutions[i].imag ))
      fd.write('\n')
    else:
      fd.write('#%s has no solution we reached xi and energyvalue: %f %f %s\n' %(str(dict),xier.xi,xier.energy,str(xier.rgvars).translate(None,'\n')) )
      pass
  return 0  
  
def create_predefined_rgeq(interactionconstant,eta,npair,nlevel , hamiltonian ,inputname,epsilon0 = 0.0):
  #gathering of the input we need to generate a RichEq object which we shall solve
  if inputname == 'r':
    energy , degeneration, seniority = wp.romboutsprob() #typical interacion constant between -0.075 and -0.0001 for groundstate typical n of pairs is 10 (quarterfilling)
  elif inputname == 's':
    energy , degeneration , seniority = wp.picketfence( nlevel)
  elif inputname == 'd':
    energy , degeneration , seniority , npair = wp.dang(filename = 'Sn120Neutrons' , cutoff = 1e5) #typ int. c. = -0.137
  else:
    energy = [0]*npair; degeneration = [2]*npair; seniority = [0]*npair #is used when runstring = 'f' those are overwritten in generating_data but we need a rgeq object so we create it with some dummy variables
  #Creation of the RichEq object
  if  hamiltonian == 'r':
    rgeq = rg.RichRedBcs(energy, degeneration, seniority , interactionconstant , npair)
  elif hamiltonian == 'f':
    rgeq = rg.RichFacInt(energy , degeneration, seniority , interactionconstant ,eta ,npair)
  elif hamiltonian == 'D':
    rgeq = rg.Dicke(energy, degeneration,interactionconstant, epsilon0,npair)
  else:
    print 'unknown Hamiltonian at the moment only the factorisable interaction (f) and reduced BCS Hamiltonian (r) are implemented'
    sys.exit(1)
  return rgeq

def input_checks(rgeq):
  assert(len(rgeq.ontaardingen) == len(rgeq.senioriteit)) 
  assert(len(rgeq.energiel) == len(rgeq.ontaardingen))
  assert(len(rgeq.energiel) > 1 )   
  assert(rgeq.apair <= sum(rgeq.ontaardingen)/2.)
  assert(rgeq.apair > 0)
  if rgeq.rgsolutions != None:
    rgeq.solve()
  print rgeq
  
def binoml(n,r):
  '''
  gives back the number of possibility's to choose k things out of a set of n things with no repetition and no order
  '''
  #first check what's the lowest r or n-r and set to r.
  if r > n-r:  # for smaller intermediate values
    r = n-r
  return int( reduce( mul, range((n-r+1), n+1), 1) /
    reduce( mul, range(1,r+1), 1) )
  
#def printtda(nlev = 12, apar = 10 , ontaardingen = [4,4,4,8,4,4,8,8,4,8,4,8,12] , senioriteit = none, startw = -1,  endw = 310000):
def printtda(nlev = 12, apar = 6 , ontaardingen = [2,2,2,2,2,2,2,2,2,2,2,2] , senioriteit = None, startw = -1,  endw = 310000):
  '''
  generates a file that contains the connection between the wind parameter (the number of the state as determined by combinations with replacement)
  and the distribution of the pairs over the sp levels at low enough g so the sp levels that contain a pair behave as the corresponding sp levels.
  '''
  tdacombinations = combinations_with_replacement(np.arange(0,nlev),apar)
  if senioriteit == none:
    senioriteit = np.zeros(nlev)
  file = open('correspondentie.dat', 'w')
  wind = 0
  for tdadict in tdacombinations:
    tdastartd = {}
    goodsol = true
    for i in tdadict:
      a = tdadict.count(i)      
      tdastartd[i] = a
      if a*2 + senioriteit[i]*2 > ontaardingen[i]:
        goodsol = false
    if goodsol == true:
      wind += 1
      if wind >= startw and wind < endw:
        print 'the %g th state corresponds to %s ' %( wind,str(tdastartd))
        file.write('%g\t%s\n' %(wind,str(tdastartd)))
  file.close()

def readlevels(infilename , waardeafh,nlevel = None,nauw = None):
  '''
  function that reads the energylevels from file infilename at dependend value wafh
  REMARK: the column at the left is the column with the dependend value that characterizes the energylevels
  at the right of wafh is the corresponding sp spectrum (same column heigth as wafh)
  '''
  waarden = []
  with open( infilename , 'r') as ifile:
    for line in ifile:
      test = line.split()
      if test[0][0] != '#':
        if (waardeafh >= float(test[0])-nauw/2. and waardeafh <= float(test[0])+nauw/2. ) or waardeafh == None:
          waarden = map(float,line.split())
          print 'energie niveaus bij afhankelijke variabele:' , waarden[0] , 'zijn ingelezen'
          break
    if ((len(waarden) < 2)):
      print 'no data retrieved from file'
    if nlevel == None  :
      energielev = waarden[1:]
    else:
      energielev = waarden[1:nlevel+1]
  return energielev      
    
def dangSn(filen,cutoff = 1e5):
  f = open(filen,'r')
  elevel = []
  deg = []
  for line in f:
    if line.startswith('#'):
      continue
    else:
      data = map(float,line.split())
      if data[0] > cutoff:
        break
      elevel.append(data[0])
      #deg.append(data[1]*4)
  deg = [2]*len(elevel)
  f.close()
  return elevel,deg   
    
def readchk(checkfile = 'rgvars.chk' , linen = -1):
  '''
  function that reads the information of the $linen line of the checkfile generated by generating_data (the file rgsolver) and returns the dependend var
  and the rgvars, at the critical point that one wants to investigate
  '''
  a = open(checkfile,'r')
  f = a.readlines()
  ll = f[linen]
  listchk = ll.split()
  begink = listchk[0]
  im = listchk.index('Im')
  rergvar = listchk[1:im]
  imrgvar = listchk[im+1:]
  rgvars = empty(len(rergvar),complex)
  rgvars.real = rergvar
  rgvars.imag = imrgvar
  a.close()
  return begink, rgvars

def continuity_check(arraysol, rgeq,crit = 1.6,dvar = 'g'):
  '''
  function that investigates, if the new solution: $nval forfilles the continuity criterium:
  that is that the mean difference of the previous $grootte solutions *$crit is bigger then the difference between $nval
  and the previous solution. If this is the case the last element of $arraysol is deleted and $nval is inserted at the end.
  If this is not the case $arraysol is not changed and the variable $con is set False
  REMARK: Very important function for (generating_datak) to check if there is no jump to an excited state.
  '''
  grootte = 6
  assert(len(arraysol) <= grootte), 'somewhere the continuity array got populated to much'
  meand = 0 ; adapter = 0 # adapter adapts for bigger steps (important around critical points)
  con = True
  if len(arraysol) < grootte:
    arraysol.append(rgeq.copy())
  else:
    for i in range(grootte-1): 
      adapter += abs(arraysol[i].getvar(dvar)- arraysol[i+1].getvar(dvar)) 
      meand += abs(arraysol[i].energy - arraysol[i+1].energy)+ 0.001 #the + 0.001 is for the case if all elements in carray are equal what happens at very low g and increasing xi
    verschil = abs(arraysol[-1].energy - rgeq.energy)
    adapter = abs(arraysol[-1].getvar(dvar) - rgeq.getvar(dvar))*(grootte-1)/adapter
    if adapter < 1: adapter = 1 
    if verschil < meand* crit/(grootte-1) *(adapter)**1.8: #REMARK 1/(grootte-1) of mean difference but also * (grootte-1) because of adapter
      del(arraysol[0])
      arraysol.append(rgeq.copy())
    else:
      con = False
  return con  
  
def str2tdadict(integerstring):
  return dict(zip( range(len(integerstring)), map(int, list(integerstring))  ))

def remove(name,regexp = False, reva =r'(/d+/./d+)' ):
  '''
  this functions removes all the files that contains a number.number syntax in the current working directory
  '''
  if regexp == True:
    regexp = re.compile(name+reva) 
  dd = [i for i in os.listdir(os.getcwd()) if regexp.match(i)]
  for a in dd:
    os.remove(a)
  
def test():
  '''
  This function is used to test functions that are defined above.
  '''
  el = np.linspace(-12,34,30)
  deg = range(2,7,2)*10
  print deg
  sen = np.zeros(30)
  npair = 15
  fl = get_fermilevel(npair,el,deg)
  print fl
  a,b,c,d,e,f = windowcreation(el,deg,sen,npair,10,10)
  print a,b,c,d,e,f,el
  
if __name__ == "__main__":
  test()
