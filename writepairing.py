# -*- coding: utf-8 -*-
#!/usr/bin/env python
import sys,math , os  , shutil 
import numpy as np
from numpy import ones, zeros ,array, sort, empty
from itertools import izip,combinations, combinations_with_replacement #for permutations in easysolve
import re , copy
#from scoop import futures #for parallel programming

from rgfunctions import *
import rgfunctions as rf 
from plotfunctions import *
import richardsongaudin as rg
import newtonraphson as nr

def main():
  """
  main function that solves the rg equations for general sp levels and degeneracies extracted from a given file 
  """
  #file with the sp energy levels
  infilename = sys.argv[1]
  #independent variable that is characteristic for the sp energy levels
  afhxas = sys.argv[2]
  afhxas = afhxas+ ' (a.u.)'
  try:
    koppelingsconstante = float(sys.argv[3])
  except:
    #negatieve koppelingsconstante invoeren (deeltjes die elkaar aantrekken)
    koppelingsconstante = -400.
  #aan te passen variabelen bij elk verschillend probleem
  probn = 100
  probname = infilename
  npair = 5
  nlevel = 10
  runstring = "f" #string that determines the loop "k" is run over interaction constant, "f" is run over a file of sp levels, 'e' probes the exited states
  typint = 'r' #to make the distinction of the different Hamiltonians we can solve with this technique
  eta = 1. #only necessary if you try to determine the solutions of a factorisable interaction
  #startwaarde afhankelijke variabele only important when runstring = 'f' REMARK: exited states always goes from small interaction constant to the strong interaction regime
  afh = {'start':0. , 'end':1. , 'step':0.01}
  kp = True #if the initial interaction constant lays in the big or small regime for the interaction constant(if we run over a file of splevels)
  #put kleinekoppeling True if the tda distribution is 1 1 1 1 ... (and the interaction constant has little effect) if False the tdastartd is apair 0 0 ...
  kleinekoppeling = False
  #the 2 values underneath are only important when runstring = "k" or runstring = 'e'
  ende = -1.
  wafh = 0.2 # value of the dependend variable where you want to investigate the dependence on the interaction constant of.
  """2wd defines the energy range for the sp levels that interact with each other. If this 
  variable is defined (wd2 > 0 ) then we look at the correlation-energy in function of a changing density of the sp states. If the variable
  is not defined (wd < 0) then the number of sp levels stays constant
  """
  wd2 = -50.
  degeneration = degeneracies_super(nlevel) 
  #set the seniority and degeneracy of the problem (in this case all to zero) REMARK: we can only do this after the definite number of sp levels is known
  senioriteit = zeros(nlevel,float) 
  #elevels,degeneration = dangSn(infilename)
  #nlevel = len(elevels) 
  #nlevel,degeneration = checkdegeneratie(elevels,degeneration,nauw = 0.0011)
  
  #some checks of the initialisation variables
  assert(len(degeneration) == len(senioriteit)) 
  assert(len(degeneration) > 1 )
  assert(koppelingsconstante < 0)
  probname += "l=%gp=%gal4" %(nlevel,npair)
  
  rgeq = None
  if runstring == 'k' or runstring == 'e':  
    elevels = readlevels(infilename,wafh,nlevel)
    nlevel,degeneration = checkdegeneratie(elevels,degeneration) #check if there are any accidental degenrations and reduce the number of sp levels accordingly and return the degeneration dictionary
    #set the seniority and degeneracy of the problem (in this case all to zero) REMARK: we can only do this after the definite number of sp levels is known
    senioriteit = zeros(nlevel,float)  
    if runstring == 'k': name = "runk1geo%s=%fname=%si=%g" %(afhxas , wafh,probname,probn)
    else: name = 'excitedstatesplotenergy1geo%s=%fname=%s' %(afhxas , wafh,probname)  
    
    if typint == 'r':
      rgeq = rg.RichRedBcs(elevels,degeneration,senioriteit,koppelingsconstante,npair,xi = 5.992003000000000027e-6)
    elif typint == 'f':
      rgeq = rg.RichFacInt(elevels,degeneration,senioriteit,koppelingsconstante,eta,npair,xi = 5.992003000000000027e-6)
    else:
      print 'error this interaction is not yet known: put typint one of these: r,f'
      sys.exit(1)
    assert(rgeq is not None) 
  #dictionary that determines the start tda solutions
  pairingd = {}
  #preparing the dictionary according to the limit of the interaction constant we probe
  if kleinekoppeling:    
    pairingd = tdadict_kleinekoppeling(npair,degeneration,senioriteit)
    step = -0.01
  else:
    pairingd[0] = npair
    step = 1.
  #####################################################################################################################2
  ###########################END OF INITIALISATION (start mainwork) ###################################################2
  if runstring == "f":
    #loop over file that contains the sp levels
    name = "nr= %dname=%sk=%fnp=%sL=%swd2=%fvo=100" %( probn , probname, koppelingsconstante,npair, nlevel ,wd2) 
    #generate seperate dir for the solved problem
    generate_dir(name,infilename,None)
    generating_data(infilename,typint,eta,koppelingsconstante,afh,npair,nlevel,wd2,pairingd,degeneration,afhxas,kp,namepf = name+'.dat')
    #generate some nice plots of the created data in the problem directory
    generate_plot(nlevel,npair,koppelingsconstante,afhxas,plotg = True,name = name+'.dat')
    
    
  elif runstring == "k":
    #generate seperate dir for the solved problem
    generate_dir(name,infilename,None)
    generating_datak(rgeq,pairingd,afhxas,step,ende,tdafilebool = False ,exname = name,moviede = False)       
    #generate some nice plots of the created data in the problem directory
    generate_plot(nlevel,npair,wafh,afhxas,plotg = False,name='plotenergy'+name+'.dat')
  
  elif runstring == 'e':
    '''
    because in the non-interacting limit we now some exited states vb: 1111100000 -> 1111010000 we must start the survey of the
    non-interacting states from weak interaction limit to the big limit
    '''
    rgw = True
    mov = False
    tdaf = False
    try:
      assert(kleinekoppeling)
    except AssertionError:
      print 'put kleinekoppeling True because the survey of exited states commands to go from the weak interaction limit to the strong interaction limit'
    probeLowestExitedStates(name,rgeq,pairingd,afhxas,step,ende,rgw,mov,tdaf) 
    #probeAllExitedStates(name,elevels,koppelingsconstante,npair,nlevel,afhxas,wafh,pairingd,degeneration,senioriteit,step,ende,rgw,mov,tdaf) :
  

def generating_data(infilename,typint,eta,koppelingsconstante,afhw,npair,nlevels,wd2,pairingdict,deg,afhxas,kp = False,namepf = 'plotenergy.dat',checkfile = True,rgwrite = True,exname = ''):
  """
  function that generates the file with all the data 
  does the main job
  """
  #startwaarde afhankelijke variabele
  afh = afhw['start']
  end = afhw['end']
  step = afhw['step']
  saveopl = 123456.
  discstep = 200. #discontinuity step if two succeeding solutions differ by more then this we search the ground state from the beginning
  #create the data files
  plotenergyfile = open(namepf,"w")
  rf.info_1set(plotenergyfile,'#interactionconstant = %f\n#deg=%s\n'%(koppelingsconstante,str(deg)) , exinfo = "#at (afh = %s)\n#g\tcE\tgE\tnig\td\tnlevels\trgvar(real)\trgvar(imag)\t ...\n" %(afhxas),tdadict = pairingdict)
  ifile = open( infilename , 'r')
  #if rgvar is None we haven't determined any rg variables so the first solution of the file has to be determined from the corresponding tda solutions (xi = 0 -> xi = 1)
  #REMARK after the first solution we have a good guess for the next sol of the file so we don't need to start from tda but can directly
  #start from the previous solution if the stepwidth of the dependent variable of the file is low enough (WATCH OUT for critical points)
  rgeq = None
  for line in ifile:
    if (line[0] is '#'):
      continue
    waarden = map(float,line.split())
    #REMARK remember waarden[0] consists of the dependent variable that characterizes the sp levels
    if math.fabs(waarden[0])+0.000001 > math.fabs(afh) and math.fabs(waarden[0])-0.000001 < math.fabs(afh) :	
      print "***************************"
      print "afhwaarde = %f" %afh
      print "***************************" 
      defvar = waarden[0]
      energielev = waarden[1:]
      nlevels,energielev = wd_processing(wd2,nlevels,npair,energielev)
      nlevel,deg = checkdegeneratie(energielev,deg) 
      sen = np.zeros(nlevel,float)
      energielev.sort()
      energielev = array(energielev)
      d = calculated(energielev,deg) 
      bb = calcnintgrondtoestand(npair,energielev,deg,sen)
      print "%f is de waarde van d bij afh = %f" %(d,defvar)
      """
      If the interaction constant is to small we can't find any solution with this pairing dictionary, therefore we call the generateRGsolKoppeling
      a function that gives back the RG energy of the groundstate and the RG variables for arbitrary interaction constants
      """      
      try:
	if rgeq is None: #We just started so we need to find a startsolution but if kp is True we win some time by immediately raising a valueError this timewinst is significant for large systems
	  if typint == 'r':
	    rgeq = rg.RichRedBcs(energielev,deg,sen,koppelingsconstante,npair,xi = 5.992003000000000027e-6)
	  elif typint == 'f':
	    rgeq = rg.RichFacInt(energielev,deg,sen,koppelingsconstante,eta,npair,xi = 5.992003000000000027e-6)
	  else:
	    print 'error this interaction is not yet known: put typint one of these: r,f'
	    sys.exit(1)
	  if kp is True:
	    raise ValueError
          energierg,rgeq = rg.RichardsonSolver(rgeq).main_solve(pairingdict,plotrgvarpath = True,xlim = None , ylim = None)	  
	else:
	  rgeq.energiel = energielev
          energierg = rgeq.solve()
      except (ValueError, np.linalg.linalg.LinAlgError) as e:
	plotenergyfile.write('# we reached a critical point\n') ; print 'we reached a critical point'
	if checkfile == True and rgeq.rgsolutions is not None:
	  chk = open('rgvars.chk','a')
	  chk.write('%f\t%s\tIm\t%s\n' %(afh,' '.join(map(str,rgeq.rgsolutions.real)),' '.join(map(str,rgeq.rgsolutions.imag))))	
	  chk.close()
	if kp is True:
	  send = rgeq.g
          rgeq = genstartsol(rgeq,d,send,begin = None,pairingd = None)
          a = False ; gg =[(0.5,10.),(0.2,4.),(3.,5.),(5.,10.),(1.,2.),(20.,10.)] ; i =0
          while(a is False and i < len(gg)):
            try:
	      energierg, rgeq,rgeq2 = littleLoop(rgeq,(send - rgeq.g)/gg[i][0],gg[i][1],complexstepd = 10000,end = send ,backxi = False,xival = 1.)
	      a = True
	    except (ValueError,np.linalg.linalg.LinAlgError) as e:
	      i += 1
	   
	else:
	  try:
	    energierg,rgeq = rg.RichardsonSolver(rgeq).main_solve(pairingdict)
	  except:
	    print '######################################################'
	    print 'ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR'
	    print ' ####################################################'
      #little discontinuity check
      if (energierg + discstep < saveopl or energierg - discstep > saveopl) and saveopl is not 123456.:
	print '######################################################################################'
	print 'we arrived at a discontinuity'
	print '######################################################################################'
	if kp is True:
	  send = rgeq.g 
          rgeq = genstartsol(rgeq,d,send,begin = None,pairingd = None)
          a = False ; gg =[(0.5,10.),(2.,5.),(3.,20.),(5,10.),(1.,2.)] ; i =0
          while(a is False):
            stepg = (send - rgeq.g)/gg[i][0]
            print send , str(rgeq),stepg
            try:
	      energierg, rgeq,rgeq2 = littleLoop(rgeq,stepg,gg[i][1],complexstepd = 10000,end = send ,backxi = False,xival = 1.)
	      a = True
	    except (ValueError,np.linalg.linalg.LinAlgError) as e:
	      i += 1
        else:
	  energierg,rgeq = rg.RichardsonSolver(rgeq).main_solve(pairingdict)
      saveopl = energierg
      if energierg is not None:
	plotenergyfile.write("%f\t%f\t%f\t%f\t%f\t%f" %(defvar, energierg-bb, energierg,bb , d,nlevel))
      if rgwrite is True:
	for i in range(npair):
	  plotenergyfile.write('\t%f' %( rgeq.rgsolutions[i].real ))
	  plotenergyfile.write('\t%f' %( rgeq.rgsolutions[i].imag ))
      plotenergyfile.write('\n')
      afh = afh + step
      if math.fabs(afh)+0.00001 > math.fabs(end) and math.fabs(afh)-0.00001 < math.fabs(end):
	break
  plotenergyfile.close()
  if rgwrite is True:
    plotrgvars(rgeq.apair,ref = namepf,afhvar = afhxas,namerg = 'rgvar%f%s' %(koppelingsconstante,exname),istart = 6 )
  ifile.close()

def generating_datak(rgeq,pairingd,dvar,step,end ,xival = 1.,rgwrite = True,exname = '',moviede = False,tdafilebool = False):
  conarray = [] # array that is used to check continuity
  plotenergyfile = open("plotenergy%s.dat" %exname,"w") #opening of the file that's going to contain the results
  if tdafilebool is True:
    tdafile = open('tdafile.dat','w')
    dirname = "moviedir%s" %exname
    print ('#making directory for the movie of the rgvars: %s' %dirname)
    dir = dirname
    if os.path.isdir(dir):
      shutil.rmtree(dir)
    os.mkdir(dir)
    olddir = os.getcwd() # remember the old directory
    os.chdir(dir) #change to the moviedirectory   
  d = rf.calculated(rgeq.energiel,rgeq.ontaardingen)  #Calculation of the mean distance between the used sp levels
  bb = rf.calcnintgrondtoestand(rgeq.apair, rgeq.energiel,rgeq.ontaardingen,rgeq.senioriteit) 
  rf.info_1set(plotenergyfile,str(rgeq) , exinfo = "#The variable we change is %s d = %f \n#and the noninteracting groundstate = %f\n#g\tcE\tgE\trgvar(real)\trgvar(imag)\t ...\n" %(dvar, d,bb),tdadict = pairingd)
  lastkp = False #boolean to know if the last g value circumvented a critical point
  #if rgeq.rgsolutions is None we haven't determined any rg variables so the first solution has to be determined from the corresponding tda solutions (xi = 0 -> xi = 1)
  #REMARK after the first solution we have a good guess for the next sol of the file so we don't need to start from tda but can directly
  #start from the previous solution if the stepwidth of the dependent variable of the file is low enough (WATCH OUT for critical points)
  if rgeq.rgsolutions is None: 
    energierg,rgeq = rg.RichardsonSolver(rgeq).main_solve(pairingd,xival = xival)   
  rgeqsaveback = None
  while (rgeq.getvar(dvar) != end ):  
    rgeq.setvar(dvar,rgeq.getvar(dvar) + step)
    if (rgeq.getvar(dvar)- abs(step)*0.8 < end and rgeq.getvar(dvar) + abs(step)*0.8 > end):
      rgeq.setvar(dvar,end)    
    print 'The variable we change: %s is: %s ' %(dvar ,str(rgeq.getvar(dvar)))
    try:
      energierg = rgeq.solve()
      if not rf.continuity_check(conarray, rgeq,crit = 1.8,dvar = dvar):
	raise ValueError
    except (ValueError, np.linalg.linalg.LinAlgError) as e:
      """
      we arrived at a critical point so now we circumvent it by making g complex.
      """
      assert(rgeq.rgsolutions is not None)
      n = 1 #integer that determines how many steps g need to increase to circumvent the critical point
      # if you use the generating_datak wrapper probeExitedStates make sure n = 1 because it takes the difference of exited states with the groundstate and
      #if the critical points of the exited state are different in contrast to those of the ground state then you'll have noisy graphs.
      plotenergyfile.write("# %s kritisch punt \n" %str(rgeq.getvar(dvar)))
      savestep = step
      send = None
      complexstep = 10000 ; extremecp = False
      while lastkp == False:	
	if n > 1:
	  step += step /3.
	elif (rgeq.getvar(dvar) - abs(step)*1.5 < end and rgeq.getvar(dvar) + abs(step) *1.5 > end):
	  send = end
        try:
          #make sure that you divide step by a dividor from step 
          savergeq = rgeq.copy()
          if rgeqsaveback == None:
	    send = rgeq.getvar(dvar) + savestep
	    energierg,rgeq,rgeqsaveback = rf.littleLoop(conarray[-2].copy(),step/10.,n*10,complexstepd = complexstep,end = send,dvar = dvar)
	  elif n < 1000:  
            energierg,rgeq,rgeqsaveback = rf.littleLoop(rgeq,step/2.,n*2,complexstepd = complexstep,end = send,dvar = dvar)
          elif n > 1000:
	    send = rgeq.getvar(dvar)+ savestep+savestep/2.
	    complexstep = 100000
	    energierg,rgeq, rgeqsaveback = rf.littleLoop(rgeqsaveback,savestep*10,2.,complexstepd = complexstep,end = send,dvar = dvar)
            lastkp = True
            step = savestep 
          if not rf.continuity_check(conarray, rgeq,crit = 1.5,dvar = dvar):
	    raise ValueError
          lastkp = True
          step = savestep
        except (ValueError, np.linalg.linalg.LinAlgError) as e:
	  rgeq = savergeq
	  if abs(rgeq.getvar(dvar)) > 1e-3:
	    complexstep *= 2.
	  else:
	    complexstep /= 10.
	  if rgeq.getvar(dvar) - abs(step)*1.5 < end and rgeq.getvar(dvar) + abs(step) *1.5 > end:
	    step /= 2.	    
	    if abs(step) < 1e-7:
	      n = 100000
	  else:  
	    step /= 10.
	    n *= 10
	  print 'couldn\'t circumvent the critical point at %s = %f, try to make the step of the interaction constant bigger in complexspace' %(dvar,rgeq.getvar(dvar))
	  print 'circumventing a critical point at %g steps and step in complexspace is %f' %(n, complexstep)       
    finally:
      if tdafilebool is True:
	try:
	  pairingdict  =  rg.RichardsonSolver(rgeq).main_desolve(xistep = -0.01,rgwrite = tdafilebool, plotrgvarpath = moviede,plotepath = False)
	  tdafile.write('%f\t%s\n'  %(rgeq.getvar(dvar),' '.join(map(str,pairingdict))))
        except (ValueError, np.linalg.linalg.LinAlgError,NameError) as e:
	  print 'problem in going back in xispace to xi = o to determine the corresponding tdadistribution'
	  tdafile.write('#%f\t%s\n'  %(rgeq.getvar(dvar),'We couldn\'t find any solutions because their occured an error in desolving the Richardson-Gaudin solutions from XI =1 to XI = 0')) 
      	  
    plotenergyfile.write("%f\t%f\t%f" %(rgeq.getvar(dvar), energierg-bb, energierg)) 
    if rgwrite is True:
      for i in range(rgeq.apair):
	plotenergyfile.write('\t%f' %( rgeq.rgsolutions[i].real ))
	plotenergyfile.write('\t%f' %( rgeq.rgsolutions[i].imag ))
    plotenergyfile.write('\n')
    lastkp = False   ; extremecp = False
  #end calculation underneath is just some cleaning up and closing files  
  plotenergyfile.close()  
  if tdafilebool is True:
    tdafile.close()
    if moviede is True:
      rf.makemovie()
    os.chdir(olddir) #because we were in the movie dir but we want our other plots in the parentdir
  if rgwrite == True:
    plotrgvars(rgeq.apair, ref = "plotenergy%s.dat" %exname, afhvar = dvar,namerg = 'rgvar%s' %exname,istart = 3) 
  return energierg,rgeq  
    
def probeLowestExitedStates(dirname,infilename,rgeq,pairingd,afhxas,step,ende,rgw,mov,tdaf) :
  '''
  function that investigates excited states
  '''
  dir = dirname
  if os.path.isdir(dir):
    shutil.rmtree(dir)
  os.mkdir(dir)
  if isinstance(infilename,str):
    shutil.copy(infilename,dir)
  os.chdir(dir)  
  pairingd = tdadict_kleinekoppeling(rgeq.apair,rgeq.ontaardingen,rgeq.senioriteit)
  generating_datak(rgeq,pairingd,afhxas,step, ende,exname = 'grondtoestand',rgwrite = rgw, moviede =mov , tdafilebool = tdaf)
  generate_plot(rgeq.alevel,rgeq.apair,rgeq.g,afhxas,name = 'plotenergygrondtoestand.dat', plotg = False)
  for j in xrange(rgeq.apair):
    del(pairingd[rgeq.apair-j-1])
    for i in xrange(rgeq.alevel-rgeq.apair):
      pairingd[rgeq.apair+i] = 1.
      extraname = str(rgeq.apair-j) +'_'+str(rgeq.apair +i+1)
      generating_datak(rgeq,pairingd,afhxas,step,ende,rgwrite = rgw,exname = extraname,moviede = mov,tdafilebool = tdaf)
      del(pairingd[rgeq.apair+i])
      generate_plot(rgeq.alevel,rgeq.apair,rgeq.g,afhxas,name = 'plotenergy%s.dat' %extraname , plotg = False)
    pairingd[rgeq.apair-j-1] = 1.  
    break #uncomment if you only want the lowest exited states
  generatePlotExited(rgeq.alevel,rgeq.apair,wafh,afhxas)

def probeAllExitedStates(dirname,rgeq,pairingd,afhxas,wafh,step,ende,rgw,mov,tdaf) :
  '''
  function that is meant to generate all the exited state of a particular reduced BCS Hamiltonian
  but first generate the groundstate so we can compare
  REMARK: make sure step is negative and koppelingsconstante is very small so we can make the interaction constant stronger
  and probe excited states (we know already the distribution of tdapairs for the excited states at small interaction constant)
  REMARK: only works when we can start from the small interaction limit
  '''
  #generate seperate dir for the solved problem
  generate_dir('allexcited%s' % dirname,None,None)
  #generate special dir for the groundstate info
  generate_dir('grondtoestand',None,None)
  pairingd = tdadict_kleinekoppeling(npair,degeneration,seniority)
  generating_datak(rgeq,pairingd,afhxas,step, ende,exname = 'grondtoestand',rgwrite = rgw,moviede = mov,tdafilebool = tdaf)
  generate_plot(rgeq.alevel,rgeq.apair, rgeq.g,afhxas,name = 'plotenergygrondtoestand.dat', plotg = False)
  os.chdir(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir)))
  #to calculate all the permutations of 1 1 0 0 ... so we choose out a np.arange(alevel) apair levels where we put our pairs (without repetition)
  tdacombinations = combinations(np.arange(0,rgeq.alevel),rgeq.apair)
  onezipper = np.ones(rgeq.apair)
  tdacor = open('tdacor.dat','w')
  tdacor.write('#This file contains the correspondence between the directorys and the start tda distributions \n #The first column is the directory number and on the same line is the tda start distribution written \n')
  i = 0
  for tdadict in tdacombinations:
    generate_dir('%g' %i,None,None) #only handy for a small amount of levels
    #tdastart needs to be a dictionary so we need to convert the list that contains one element of the permutation sequence to a dictionary    
    tdastartd = dict(zip(tdadict,onezipper))
    tdacor.write('%g\ttdadict= %s\n' %(i,' '.join(map(str,tdastartd.values()))))
    print 'we start generating_datak with: ', tdastartd
    generating_datak(rgeq,tdastartd,afhxas,step,ende  ,rgvars = None,rgwrite = rgw,exname = '',moviede = mov,tdafilebool = tdaf)
    generate_plot(rgeq.alevel,rgeq.apair,waardeafh,afhxas,plotg = False)
    os.chdir(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir)))
    i += 1
  tdacor.close()
  
  
def testcircumvent():
  '''
  test for the rgsolutions when xi is not one and/ or g is imaginary
  '''
  eendlev = np.arange(1,13)
  seniority = zeros(len(eendlev),float)
  degeneration = np.ones(len(eendlev))*2.
  alevel = len(eendlev)
  apair = alevel/2
  tdastartd = {}
  tdastartd = {0:6}
  enddatak = -0.001
  stepg = 0.001
  rgvar = None
  afhxas = 'g'
  waardeafh = 0.
  #generate seperate dir for the solved problem
  generate_dir('testcircumventwithlowxicomplexg3',None,None)
  print tdastartd
  xival = 0.1
  rgeq = RichRedBcs(eendlev,degeneration,seniority,g,apair)
  while xival <= 1.:
    g = -1.0007
    energierg,rgeq = rg.RichardsonSolver(rgeq).main_solve(tdastartd,gwrite = False,xistep = 0.01,xival=xival)   
    while g.imag < abs(g.real/10.):
      rgeq.g += 1j*0.01
      energierg,rgvar = rgeq.solve()
      print rgeq.g
    try:
      generating_datak(rgeq,tdastartd,afhxas,stepg,enddatak ,rgvars = rgvar,rgwrite = True,exname = '%f' %xival,moviede = False,tdafilebool = False,xival = xival)
      generate_plot(rgeq.alevel,rgeq.apair,waardeafh,afhxas,plotg = False,name = 'plotenergy%f.dat' %xival)
    except:
      pass   
    xival += 0.1

def testeasysolve():
  '''
  calculates all the eigenvalues of the pairingsHamiltonian of a system with double degenerate, equidistant levels and zero seniority
  '''
  eendlev = np.arange(1,13)
  eendlev = array([0.04,0.08,0.16,0.20,0.32,0.36,0.40,0.52,0.64,0.68,0.72,0.8,1]) 
  eendlev = np.sqrt(eendlev)
  seniority = zeros(len(eendlev),float)
  degeneration = np.ones(len(eendlev))*2
  degeneration = [4,4,4,8,4,4,8,8,4,8,4,8,12]
  alevel = len(eendlev)
  apair = alevel/2
  apair = 10
  afhxas = 'g'
  waardeafh = 0.
  eta =1.
  generate_dir('romboutstest',None,None)
  #to calculate all the permutations of 1 1 0 0 ... so we choose out a np.arange(alevel) apair levels where we put our pairs (without repetition)
  tdacombinations = combinations_with_replacement(np.arange(alevel),apair)
  onezipper = np.ones(apair)
  tdacor = open('tdacor.dat','w')
  tdacor.write('#This file contains the correspondence between the directorys and the start tda distributions \n #The first column is the directory number and on the same line is the tda start distribution written \n')
  i = 0
  for a in [0]:
    if a == 0:
      g = -0.0001 ; enddatak = -0.075 ;   stepg = -0.0001
    else:
      g = 0.0001 ; enddatak = 1. ; stepg = 0.003
    rgeq = rg.RichFacInt(eendlev,degeneration,seniority,g,eta,apair)
    for tdadict in tdacombinations:
      tdastartd = {}
      goodsol = True
      for j in tdadict:
        a = tdadict.count(j)      
        tdastartd[j] = a
        if a*2 + rgeq.senioriteit[j]*2 > rgeq.ontaardingen[j]:
	  goodsol = False
      if goodsol == False:
	continue
      if i > 691:
	generate_dir('%g' %i,None,None) #only handy for a small amount of levels
	#tdastart needs to be a dictionary so we need to convert the list that contains one element of the permutation sequence to a dictionary    
	tdacor.write('%g\ttdadict= %s\n' %(i,' '.join(map(str,tdadict))))
	print 'we start generating_datak with: ', tdastartd
	generating_datak(rgeq,tdastartd,afhxas,stepg,enddatak ,rgwrite = True,exname = '',moviede = False,tdafilebool = False)
	generate_plot(alevel,apair,waardeafh,afhxas,plotg = False)
	os.chdir(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir)))
      i += 1
    tdacor.close()

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
  

def dangmain():
  '''
  main to describe the tin pairing spectrum recieved by dang: http://ribf.riken.go.jp/~dang/publication.html
  '''
  #initialisation variables
  filename = sys.argv[1]
  cutoff = 1e5
  koppelingsconstante = -0.137
  energielev,degeneration = dangSn(filename,cutoff)
  nlevel,degeneration = checkdegeneratie(energielev,degeneration,nauw = 0.0011)
  npair = 35
  rgvard = None
  senioriteit = [0]*len(energielev)  
  pairingd = {0:35}
  energielev,degeneration,senioriteit,npair,nlevel,fermil,extrae = windowcreation(energielev,degeneration,senioriteit,npair,5,5)
  #senioriteit = [3,3,0,0,0,0]
  pairingd = tdadict_kleinekoppeling(npair,degeneration,senioriteit)
  #apair = 7
  print energielev,degeneration,senioriteit,npair,nlevel,fermil
  #pairingd = {4:5}
  #[-10.4576,-8.4804,-7.6512,-7.7025] #ijzer
  #rgvard = readrgvars(-0.310001,'rgvarDang.dat')
  generate_dir('dangsentest','Sn120Neutrons',None)  #generate seperate dir for the solved problem
  step = -0.001
  ende = -10
  afhxas = 'g' 
  tel = 0
  eta = 1.
  wafh = 120
  d = calculated(energielev,degeneration)
  rgeq = RichRedBcs(energielev,degeneration,senioriteit,koppelingsconstante,npair)
  #while(energielev[-1] >= 0):
  #generating_datak(rgeq,pairingd,afhxas,step,ende,rgvars = rgvard,tdafilebool = False ,exname = 'Dang%g' %tel)       
  #generate_plot(len(energielev),npair,wafh,afhxas,plotg = False,name='plotenergyDang%g.dat' %tel)
  """
    tel += 1
    del(energielev[-1])
    del(senioriteit[-1])
    del(degeneration[-1])
  """ 
  #we make the problem smaller by using a cutoff energy
  
  #totsen = np.arange(8.,npair*2+1,2.)
  totsen = [12.]
  for vsen in totsen: 
    seniority_enhancer_allstates(rgeq,'DangSn120neutronwindow(5,5)',vsen,exewaarde = extrae,begin = -0.001 ,step = -0.0001)

def facintmain():
  #artikel Stefan test
  eendlev = array([0.04,0.08,0.16,0.20,0.32,0.36,0.40,0.52,0.64,0.68,0.72,0.8,1])
  eendlev = np.sqrt(eendlev)
  seniority = zeros(len(eendlev),float)
  degeneration = [4,4,4,8,4,4,8,8,4,8,4,8,12]
  alevel = len(eendlev)
  apair = 10
  g = -0.075
  eta = 1.
  tdastartd = {}
  tdastartd = {0:10,1:0,2:0}
  enddatak = -0.0002
  stepg = 0.0001
  dvar = 'g'
  rgeq = rg.RichFacInt(eendlev,degeneration,seniority,g,eta,apair)
  generate_dir('stefantda',None,None)
  generating_datak(rgeq,tdastartd,dvar,stepg,enddatak ,rgwrite = True,exname = '',moviede =True,tdafilebool = True)
  generate_plot(alevel,apair,0.,dvar,plotg = False)  

if __name__ == "__main__":
  #test()
  #main()
  #dangmain()
  #testcircumvent()
  testeasysolve()
  #addlevel()
  #facintmain()