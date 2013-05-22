# -*- coding: utf-8 -*-
#!/usr/bin/env python
import sys,math , os  , shutil 
import numpy as np
from numpy import ones, zeros ,array, sort, empty
from itertools import izip,combinations, combinations_with_replacement #for permutations in easysolve
import re , copy
from optparse import OptionParser
#from scoop import futures #for parallel programming

from rgfunctions import *
from plotfunctions import *
import richardsongaudin as rg
import newtonraphson as nr

def main():
  """
  main function that solves the rg equations for general sp levels and degeneracies extracted from a given file 
  """
  kp =True #for runstring = 'f' if the initial interaction constant lays in the big or small regime for the interaction constant(if we run over a file of splevels)
  #startwaarde afhankelijke variabele only important when runstring = 'f' REMARK: exited states always goes from small interaction constant to the strong interaction regime
  probname = 'testchanges' #adds to the directory name where all the output data is going to be stored
  afh = {'start':0. , 'end':1. , 'step':0.01} ; spkar = 'xo' #characterizes the sp levels in the file as run = 'f' and input is a filename
  step = {'n': -0.0001,'p': 0.003}; ende = {'p' : 10. , 'n' : -1.0750} #for runstring = e or k , the sign operator determines if the interaction constant is negative (n) or positive (p)
  rgw = True ; mov = False ; tdaf = False
  sign = 'n' #for runstring = k determines which instances of step and ende it receives
  npair =6#if you use a filename as input or gives the number of pairs in the rombout system
  nlevel = 12#if you use a filename as input or gives the number of levels in the picketfence model
  eta = 1.
  seniority = zeros(nlevel,float)  #if you use filename as input
  degeneration = degeneracies_super(nlevel) #if you use filename as input
  #set the seniority and degeneracy of the problem (in this case all to zero) REMARK: we can only do this after the definite number of sp levels is known
  kleinekoppeling =True#put kleinekoppeling True if the tda distribution is 1 1 1 1 ... (and the interaction constant has little effect) if False the tdastartd is apair 0 0 ...
  wd2 = -50. #(wd2 > 0 ) then we look at the correlation-energy in function of a changing density of the sp states. If the variable is not defined (wd < 0) then the number of sp levels stays constant
  #handling of input variables
  usage = "python writepairing.py [options] arguments \n(take also a look at which function is the main function in the bottom of the writepairing.py file)"
  parser = OptionParser(usage = usage)
  parser.add_option('-f','--filename',dest = 'filename' ,default = None,type = str, help = 'File of which each line consists of a dependend variable and the corresponding sp energylevels (seperated by tabs)')
  parser.add_option('-v','--variable',dest = 'depvar',default = 'g' , help = 'Give the string representation of the independend variable that is going to be changed')
  parser.add_option('-i','--interactionc',dest = 'interactionconstant',default = -0.0001,type = float,help = 'Give the pairing strength')
  parser.add_option('-r','--run', dest = 'runstring' , default = 'k' , help = 'Give the mainrun you want the program to execute')
  parser.add_option('-H' , '--hamiltonian' , dest = 'hamiltonian' , default = 'r' , help = 'Give the hamiltonian you want to solve: \'r\' is the red. bcs. Hamiltonian, \'f\' is the fac. int. Hamiltonian')
  parser.add_option('-n', '--inputname', dest = 'inputname' , default = None , help = 'If you don\'t use a file to get the sp levels choose some predefined sets: r (rombouts:2010) , s (sambataro:2008)' )
  parser.add_option('-w','--afhwaarde' , dest = 'wafh' , default = 0. , help = 'If you use a file containing sp energy levels as input, wafh is the value that characterizes the sp levels you want to investigate')
  #parser.add_option('-p','--apair' , dest = 'apair' , default = 10 , help ='the number of pairs ', type = int)
  (options , args) = parser.parse_args(sys.argv[1:])
  filename = options.filename
  depvar = options.depvar
  interactionconstant = options.interactionconstant #negatieve koppelingsconstante invoeren (deeltjes die elkaar aantrekken)
  runstring = options.runstring
  typeint = options.hamiltonian
  inputname = options.inputname
  wafh = options.wafh
  #creation of directory name where we save the output of the run
  if args: 
    probname = args[0] +'prob'
  inputlist = ['r','s','d']
  assert(filename != None or inputname in inputlist)
  if filename is not None: inp = filename
  else: inp = inputname
  name = "wp:run=%stype=%sinput=%sp=%gl=%gvar=%sinteract=%fname=%s" %(runstring,typeint,inp,npair,nlevel,depvar,interactionconstant,probname) 
  #creation of the Richardson-Equations we are going to solve
  if filename != None:
    elevels = readlevels(filename,wafh,nlevel)
    nlevel,degeneration = checkdegeneratie(elevels,degeneration) #check if there are any accidental degenrations and reduce the number of sp levels accordingly and return the degeneration dictionary
    if typeint == 'r':
     rgeq = rg.RichRedBcs(elevels,degeneration,seniority,interactionconstant,npair,xi = 5.992003000000000027e-6)
    elif typeint == 'f':
      rgeq = rg.RichFacInt(elevels,degeneration,seniority,interactionconstant,eta,npair,xi = 5.992003000000000027e-6)
    else:
      print 'error this interaction is not yet known: put typeint one of these: r,f'
      sys.exit(1)
  elif inputname in inputlist:
    print inputname
    if inputname == 'r':
      rgeq = romboutsprob(g = interactionconstant, apair = npair) #typical interacion constant between -0.075 and -0.0001 for groundstate
    elif inputname == 's':
      rgeq = picketfence(hamiltonian = options.hamiltonian, alevel = nlevel,g = interactionconstant,eta = eta , apair = npair)
    elif inputname == 'd':
      rgeq = dang(filename = 'Sn120Neutrons' , cutoff = 1e5, koppelingsconstante = interactionconstant) #typ int. c. = -0.137
  else:
    print 'bad input filename = %s and inputname = %s' %(filename , inputname)
  
  #dictionary that determines the start tda solutions only import for runstring = 'f' (when kp is True) or 'k' 
  pairingd = {}
  #preparing the dictionary according to the limit of the interaction constant we probe
  if kleinekoppeling:    
    pairingd = tdadict_kleinekoppeling(rgeq.apair,rgeq.ontaardingen, rgeq.senioriteit)
  else:
    pairingd[0] = rgeq.apair
  
  #some checks of the initialisation variables
  assert(len(rgeq.ontaardingen) == len(rgeq.senioriteit)) 
  assert(len(rgeq.energiel) > 1 )   
  #####################################################################################################################
  ###########################END OF INITIALISATION (start mainwork) ####################################################  
  #generate seperate dir for the solved problem
  assert(rgeq is not None) 
  generate_dir(name,filename,None) #if filename is None nothing will be copied in the directory see the implementation of generate_dir in rgfunctions
  
  if runstring == "f":
    #loop over file that contains the sp levels
    generating_data(rgeq,filename,afh,wd2,pairingd,spkar,kp,namepf = name+'.dat')
    #generate some nice plots of the created data in the problem directory
    generate_plot(nlevel,npair,depvar,plotg = True,name = name+'.dat')
        
  elif runstring == "k":
    generating_datak(rgeq,pairingd,depvar,step[sign],ende[sign],rgwrite = rgw ,tdafilebool = tdaf,exname = name,moviede = mov)       
    #generate some nice plots of the created data in the problem directory
    generate_plot(nlevel,npair,depvar,plotg = False,name='plotenergy'+name+'.dat')
  
  elif runstring == 'e':
    '''
    because in the non-interacting limit we know some exited states vb: 1111100000 -> 1111010000 we must start the survey of the
    non-interacting states from weak interaction limit to the big limit
    '''
    try:
      assert(kleinekoppeling)
    except AssertionError:
      print 'put kleinekoppeling True because the survey of exited states commands to go from the weak interaction limit to the strong interaction limit'
    allstatesgenerating_datak(rgeq,depvar,step,ende,rgw,mov,tdaf)
    #probeLowestExitedStates(rgeq,depvar,step,ende,rgw,mov,tdaf) 
  else:
    print 'runstring: %s is not a valid runstring' %runstring
   
def generating_data(rgeq,infilename,afhw,wd2,pairingdict,afhxas,kp = False,namepf = 'plotenergy.dat',checkfile = True,rgwrite = True,exname = ''):
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
  info_1set(plotenergyfile,'#interactionconstant = %f\n#deg=%s\n'%(rgeq.g,str(rgeq.ontaardingen)) , exinfo = "#at (afh = %s)\n#g\tcE\tgE\tnig\td\tnlevels\trgvar(real)\trgvar(imag)\t ...\n" %(afhxas),tdadict = pairingdict)
  ifile = open( infilename , 'r')
  #if rgvar is None we haven't determined any rg variables so the first solution of the file has to be determined from the corresponding tda solutions (xi = 0 -> xi = 1)
  #REMARK after the first solution we have a good guess for the next sol of the file so we don't need to start from tda but can directly
  #start from the previous solution if the stepwidth of the dependent variable of the file is low enough (WATCH OUT for critical points)
  for line in ifile:
    if (line[0] is '#'):
      continue
    waarden = map(float,line.split())
    #REMARK remember waarden[0] consists of the dependend variable that characterizes the sp levels
    if math.fabs(waarden[0])+0.000001 > math.fabs(afh) and math.fabs(waarden[0])-0.000001 < math.fabs(afh) :        
      print "***************************"
      print "afhwaarde = %f" %afh
      print "***************************" 
      defvar = waarden[0]
      energielev = waarden[1:]
      nlevels,energielev = wd_processing(wd2,rgeq.alevel,rgeq.apair,energielev)
      nlevel,deg = checkdegeneratie(energielev,rgeq.ontaardingen) 
      energielev.sort()
      energielev = array(energielev)
      rgeq.energiel = energielev ; rgeq.ontaardingen = deg
      d = calculated(rgeq.energiel,rgeq.ontaardingen) 
      bb = calcnintgrondtoestand(rgeq)
      print "%f is de waarde van d bij afh = %f" %(d,defvar)
      """
      If the interaction constant is to small we can't find any solution with this pairing dictionary, therefore we call the generateRGsolKoppeling
      a function that gives back the RG energy of the groundstate and the RG variables for arbitrary interaction constants
      """      
      try:
        if rgeq.energy is None: #We just started so we need to find a startsolution but if kp is True we win some time by immediately raising a valueError this timewinst is significant for large systems
          energierg,rgeq = rg.RichardsonSolver(rgeq).main_solve(pairingdict,plotrgvarpath = False,xlim = None , ylim = None)          
        else:
          rgeq.energiel = energielev
          energierg = rgeq.solve()
      except (ValueError, np.linalg.linalg.LinAlgError) as e:
        plotenergyfile.write('# we reached a critical point\n') ; print 'we reached a critical point'
        if checkfile == True and rgeq.rgsolutions is not None:
          chk = open('rgvars.chk','a')
          chk.write('%f\t%s\tIm\t%s\n' %(afh-step,' '.join(map(str,rgeq.rgsolutions.real)),' '.join(map(str,rgeq.rgsolutions.imag))))        
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
        plotenergyfile.write('#we arrived at a discontinuity\n')
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
        for i in range(rgeq.apair):
          plotenergyfile.write('\t%f' %( rgeq.rgsolutions[i].real ))
          plotenergyfile.write('\t%f' %( rgeq.rgsolutions[i].imag ))
      plotenergyfile.write('\n')
      afh = afh + step
      if math.fabs(afh)+0.00001 > math.fabs(end) and math.fabs(afh)-0.00001 < math.fabs(end):
        break
  plotenergyfile.close()
  if rgwrite is True:
    plotrgvars(rgeq.apair,ref = namepf,afhvar = afhxas,namerg = 'rgvar%f%s' %(rgeq.g,exname),istart = 6 )
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
  d = calculated(rgeq.energiel,rgeq.ontaardingen)  #Calculation of the mean distance between the used sp levels
  bb = calcnintgrondtoestand(rgeq) 
  info_1set(plotenergyfile,str(rgeq) , exinfo = "#The variable we change is %s d = %f \n#and the noninteracting groundstate = %f\n#g\tcE\tgE\trgvar(real)\trgvar(imag)\t ...\n" %(dvar, d,bb),tdadict = pairingd)
  lastkp = False #boolean to know if the last g value circumvented a critical point
  #if rgeq.rgsolutions is None we haven't determined any rg variables so the first solution has to be determined from the corresponding tda solutions (xi = 0 -> xi = 1)
  #REMARK after the first solution we have a good guess for the next sol of the file so we don't need to start from tda but can directly
  #start from the previous solution if the stepwidth of the dependent variable of the file is low enough (WATCH OUT for critical points)
  if rgeq.rgsolutions is None: 
    energierg,rgeq = rg.RichardsonSolver(rgeq).main_solve(pairingd,xival = xival)   
  rgeqsaveback = None
  while (rgeq.getvar(dvar) != end ):  
    savergsolutions = np.copy(rgeq.rgsolutions) ; savedepvar = rgeq.getvar(dvar) #important if continuity check fails
    rgeq.setvar(dvar,rgeq.getvar(dvar) + step)
    if (rgeq.getvar(dvar)- abs(step)*0.8 < end and rgeq.getvar(dvar) + abs(step)*0.8 > end):
      rgeq.setvar(dvar,end)    
    print 'The variable we change: %s is: %s ' %(dvar ,str(rgeq.getvar(dvar)))
    try:
      energierg = rgeq.solve()
      if not continuity_check(conarray, rgeq,crit = 1.8,dvar = dvar):
        plotenergyfile.write('#we arrived at a discontinuity\n')
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
      send = savesend = rgeq.getvar(dvar) 
      rgeq.setvar(dvar,rgeq.getvar(dvar) - step) ; rgeq.rgsolutions = savergsolutions  #go back to the rgeq where we know the solution
      complexstep = 10000 ; extremecp = False
      while lastkp == False:        
        if (rgeq.getvar(dvar) - abs(step)*1.5 < end and rgeq.getvar(dvar) + abs(step) *1.5 > end):
          send = end
          step /= 2.
        try:
          #make sure that you divide step by a divisor from step 
          if n > 1:
            send += savestep/2.
            step += savestep/4.
          if rgeqsaveback == None:
            send += step
            energierg,rgeq,rgeqsaveback = littleLoop(conarray[-2],step/2.,n*2,complexstepd = complexstep,end = send,dvar = dvar)
          elif abs(rgeqsaveback.getvar(dvar)-send)/abs(savestep)   > 4 and n < 100: 
            energierg,rgeq,rgeqsaveback = littleLoop(rgeq,step/2.,n*2,complexstepd = complexstep,end = send,dvar = dvar)
          else:
            if n < 310:
              send = savesend 
              complexstep = 10000
              step = savestep
            send += savestep/2.
            assert(isinstance(rgeqsaveback.g,complex))
            energierg,rgeq, rgeqsaveback = littleLoop(rgeqsaveback,step/2.,n*2,complexstepd = complexstep,end = send,dvar = dvar)
            lastkp = True
            if n > 1000:
              print 'fatal error in program to many iterations in generating_datak' ; sys.exit(1) 
          if not continuity_check(conarray, rgeq,crit = 1.9 * n % 4 ,dvar = dvar):
            print'problems with continuity of the found solutions %s with the following richardsoneq %s' %(str(conarray),str(rgeq))
            plotenergyfile.write('# discontinuity at %s  , with n = %g , send = %f , complexstep = %f , step = %f' %(str(rgeq.getvar(dvar)), n , send , complexstep , step))
            rgeq.rgsolutions = savergsolutions
            rgeq.setvar(dvar,savedepvar) #go back to the rgeq where we know the solution
            raise ValueError
          lastkp = True
          step = savestep
        except (ValueError, np.linalg.linalg.LinAlgError) as e:
          n *= 4
          if abs(rgeq.getvar(dvar)) > 1e-3:
            if complexstep <= 4e5:
              complexstep *= 2. #how bigger complexstep how smaller step in complexspace
          else:
            complexstep /= 10.
          print 'circumventing a critical point at %g steps and step in complexspace is %f , at %s = %f ' %(n, complexstep,dvar,rgeq.getvar(dvar))       
    finally:
      if tdafilebool is True:
        rgsolver = rg.RichardsonSolver(rgeq)
        try:
          pairingdict  =  rgsolver.main_desolve(xistep = -0.01,rgwrite = tdafilebool, plotrgvarpath = moviede,plotepath = False)
          tdasol = rgsolver.get_tda(None) 
          tdafile.write('%f\t%s\t%s\n'  %(rgeq.getvar(dvar),' '.join(map(str,pairingdict)), ' '.join(map(str,(tdasol)))))
        except (ValueError, np.linalg.linalg.LinAlgError,NameError , rg.XiError) as e:
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
      makemovie()
    os.chdir(olddir) #because we were in the movie dir but we want our other plots in the parentdir
  if rgwrite == True:
    plotrgvars(rgeq.apair, ref = "plotenergy%s.dat" %exname, afhvar = dvar,namerg = 'rgvar%s' %exname,istart = 3) 
  return energierg,rgeq  
    
def probeLowestExitedStates(rgeq,afhxas,step,ende,rgw,mov,tdaf,sign = 'n') :
  '''
  function that investigates excited states
  '''
  print rgeq
  pairingd = tdadict_kleinekoppeling(rgeq.apair,rgeq.ontaardingen,rgeq.senioriteit)
  generating_datak(rgeq,pairingd,afhxas,step[sign], ende[sign],exname = 'grondtoestand',rgwrite = rgw, moviede =mov , tdafilebool = tdaf)
  generate_plot(rgeq.alevel,rgeq.apair,afhxas,name = 'plotenergygrondtoestand.dat', plotg = False)
  nullevel =  list(pairingd)[-1] + 1 #nullevel is the first not occupied level in pairingd (remember tdadict_kleinekoppeling gives us a dict back that only contains the lowest  possibley occupied levels)
  for j in xrange(nullevel):
    pairingd[nullevel-j-1] -= 1
    for i in xrange(rgeq.alevel-nullevel+1):
      print 'we are calculating the following excitation: pair in level %g -> level %g ' %(j,nullevel+i)
      pairingd[nullevel+i] = 1
      extraname = str(nullevel-j-1) +'_'+str(nullevel +i)
      generating_datak(rgeq,pairingd,afhxas,step[sign],ende[sign],rgwrite = rgw,exname = extraname,moviede = mov,tdafilebool = tdaf)
      pairingd[nullevel+i] -= 1
      generate_plot(rgeq.alevel,rgeq.apair,afhxas,name = 'plotenergy%s.dat' %extraname , plotg = False)
    pairingd[nullevel-j-1] += 1.  
    break #uncomment if you only want the lowest exited states
  generatePlotExited(rgeq.alevel,rgeq.apair)

def allstatesgenerating_datak(rgeq,afhxas,step,ende,rgw,mov,tdaf):
  '''
  calculates all the eigenvalues of the pairingsHamiltonian
  '''
  tdacombinations = combinations_with_replacement(np.arange(rgeq.alevel),rgeq.apair)
  onezipper = np.ones(rgeq.apair)
  tdacor = open('tdacor.dat','w')
  tdacor.write('#This file contains the correspondence between the directorys and the start tda distributions \n #The first column is the directory number and on the same line is the tda start distribution written \n')
  i = 0
  for a in [0,1]:
    if a == 0:
      enddatak = ende['n'];   stepg = step['n']
      rgeq.g = -0.0001
    else:
      enddatak = ende['p'] ; stepg = step['p']
      rgeq.g = 0.0001
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
      if i > -1: #parameter to restart easely after failure
        generate_dir('%g' %i,None,None) #only handy for a small amount of levels
        #tdastart needs to be a dictionary so we need to convert the list that contains one element of the permutation sequence to a dictionary    
        tdacor.write('%g\ttdadict= %s\n' %(i,' '.join(map(str,tdadict))))
        print 'we start generating_datak with: ', tdastartd
        generating_datak(rgeq,tdastartd,afhxas,stepg,enddatak ,rgwrite = rgw,exname = '',moviede = mov,tdafilebool = tdaf)
        generate_plot(rgeq.alevel,rgeq.apair,afhxas,plotg = False)
        os.chdir(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir)))
      i += 1
      print i
    tdacor.close()

def romboutsprob(g = -0.075, apair =10):
  #definition of problem of artikel stefan:2010 
  eendlev = array([0.04,0.08,0.16,0.20,0.32,0.36,0.40,0.52,0.64,0.68,0.72,0.8,1])
  eendlev = np.sqrt(eendlev)
  seniority = zeros(len(eendlev),float)
  degeneration = [4,4,4,8,4,4,8,8,4,8,4,8,12]
  eta = 1.
  rgeq = rg.RichFacInt(eendlev,degeneration,seniority,g,eta,apair)  
  return rgeq
  
def picketfence(hamiltonian = 'r' , alevel = 12,g = -1.,eta = 1., apair = None):
  #little toy problem introduced by sambataro:2008
  eendlev = np.arange(1,alevel+1)
  seniority = zeros(alevel,float)
  degeneration = np.ones(alevel)*2
  if apair == None:
    apair = alevel/2
  if hamiltonian == 'r':
    rgeq = rg.RichRedBcs(eendlev,degeneration, seniority,g,apair)
  elif hamiltonian == 'f':
    rgeq = rg.RichFacInt(eendlev,degeneration, seniority,g,eta,apair)
  return rgeq

def dang(filename = 'Sn120Neutrons' , cutoff = 1e5, koppelingsconstante = -0.137):
  energielev,degeneration = dangSn(filename,cutoff)
  nlevel,degeneration = checkdegeneratie(energielev,degeneration,nauw = 0.0011)
  senioriteit = [0]*len(energielev)  
  npair = 35 
  energielev,degeneration,senioriteit,npair,nlevel,fermil,extrae = windowcreation(energielev,degeneration,senioriteit,npair,5,5)
  if ham == 'r':
    rgeq = RichRedBcs(eendlev,degeneration, seniority,g,apair)
  elif ham == 'f':
    eta = 1.
    rgeq = RichFacInt(energielev,degeneration, senioriteit,koppelingsconstante,eta,npair)
  return rgeq  
  
def facintmain():
  #definition of problem of artikel stefan:2010 
  rgeq = romboutsprob()
  tdastartd = {0:10,1:0,2:0}
  enddatak = -0.0002
  stepg = 0.0001
  dvar = 'g'
  generate_dir('stefantda',None,None)
  generating_datak(rgeq,tdastartd,dvar,stepg,enddatak ,rgwrite = True,exname = '',moviede =True,tdafilebool = True)
  generate_plot(rgeq.alevel,rgeq.apair,dvar,plotg = False)  

def allstatesoneg(npair = 3):
  rgeq = romboutsprob(g = -0.1)
  rgeq = picketfence(hamiltonian = 'r' , alevel = 6,g = -2.00,eta = 1.)
  generate_dir('picketfencefac6lev3redbcs-2.00',None,None)  #generate seperate dir for the solved problem
  activelevels = 6
  '''
  totsen = np.arange(0,npair*2+1,2,int)
  for vsen in totsen: 
    seniority_enhancer_allstates(rgeq,'romboutsall',vsen,exewaarde = 0.,begin = -0.0001 ,step = -0.0001)  
  '''
  fd = open('allstatesactivelevels=%g' %activelevels,'w')
  fd.write(str(rgeq))
  ontaarding = 1.
  dataanalyse = {'rgw': True , 'ple': False, 'plrg': False , 'ont': False }
  allstatesstrongg(rgeq,fd,ontaarding,activelevels,extrae = [],exe= 0 , dataanalyse = dataanalyse) 
  plot_spectrumxichange('.','.dat') 

def testcircumvent():
  '''
  test for the rgsolutions when xi is not one and/ or g is imaginary
  '''
  rgeq = picketfence(hamiltonian = 'r' , alevel = 12,g = -1.,eta = 1.)
  tdastartd = {0:6}
  enddatak = -0.001
  stepg = 0.001
  rgvar = None
  afhxas = 'g'
  #generate seperate dir for the solved problem
  generate_dir('testcircumventwithlowxicomplexg3',None,None)
  print tdastartd
  xival = 0.1
  while xival <= 1.:
    g = -1.0007
    energierg,rgeq = rg.RichardsonSolver(rgeq).main_solve(tdastartd,gwrite = False,xistep = 0.01,xival=xival)   
    while g.imag < abs(g.real/10.):
      rgeq.g += 1j*0.01
      energierg,rgvar = rgeq.solve()
      print rgeq.g
    try:
      generating_datak(rgeq,tdastartd,afhxas,stepg,enddatak ,rgvars = rgvar,rgwrite = True,exname = '%f' %xival,moviede = False,tdafilebool = False,xival = xival)
      generate_plot(rgeq.alevel,rgeq.apair,afhxas,plotg = False,name = 'plotenergy%f.dat' %xival)
    except:
      pass   
    xival += 0.1  

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
  #generate_plot(len(energielev),npair,afhxas,plotg = False,name='plotenergyDang%g.dat' %tel)
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


def testrestart():
  dfile = 'plotenergytest.dat'
  afhvar = -0.3601
  rgvars =  readrgvarsplote(afhvar, name = dfile)
  rgeq = picketfence(hamiltonian = 'f' , alevel = 12,g = afhvar,eta = 1.) 
  rgeq.rgsolutions = rgvars ; rgeq.xi = 1. ; rgeq.solve()
  print rgeq
  step =  0.0001
  end = -0.0001
  dvar = 'g'
  name = ''
  pairingd = {0:rgeq.apair}
  generating_datak(rgeq,pairingd,dvar,step,end ,xival = 1.,rgwrite = True,exname = '',moviede = False,tdafilebool = False)
  generate_plot(rgeq.alevel,rgeq.apair,dvar,plotg = False,name='plotenergy'+name+'.dat')



if __name__ == "__main__":
  #testrestart()
  main()
  #dangmain()
  #testcircumvent()
  #addlevel() #function in rgfunctions
  #facintmain()
  #allstatesoneg()
