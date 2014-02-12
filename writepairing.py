# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Mario Van Raemdonck, 2013;
# (c) Ghent University, 2013
# -*- coding: utf-8 -*-
#!/usr/bin/env python
import sys,math,os , shutil 
import numpy as np
from numpy import ones, zeros ,array, sort, empty
from itertools import izip,combinations, combinations_with_replacement #for permutations in easysolve
import re , copy
from optparse import OptionParser
#from scoop import futures #for parallel programming

from rgfunctions import *
import rgfunctions as rgf
from plotfunctions import *
import richardsongaudin as rg
import datareader as dr

def parse_commandline():
  #handling of commandline variables #negatieve koppelingsconstante invoeren (deeltjes die elkaar aantrekken) 
  usage = "python writepairing.py [options] arguments \nOr use an inputfile (the header of outputfiles is also a valid inputformat, to know more about the structure of the input file check the datareader.py file (especially the Reader classes))"
  parser = OptionParser(usage = usage)
  parser.add_option('-f','--filename',dest = 'filename' ,default = None,type = str, help = 'File of which each line consists of a dependend variable and the corresponding sp energylevels (seperated by tabs)')
  parser.add_option('-v','--variable',dest = 'depvar',default = 'g' , help = 'Give the string representation of the independend variable that is going to be changed')
  parser.add_option('-i','--interactionc',dest = 'interactionconstant',default = -0.001,type = float,help = 'Give the pairing strength')
  parser.add_option('-e','--eta',dest = 'eta',default = 1.,type = float,help = 'Give the weight of the single particle part of the factorisable interacion Hamiltonian')
  parser.add_option('-l','--levels', dest = 'levels' , default = None , help = 'Give the number of levels', type = int)
  parser.add_option('-p','--pairs',dest = 'pairs',default = None,help = 'Give the number of pairs', type = int)
  parser.add_option('-r','--run', dest = 'runstring' , default = 'k' , help = 'Give the mainrun you want the program to execute (f file , k one parameter changes , e excited states), r (restart from an outputfile or a set of outputfiles under current dir)')
  parser.add_option('-H' , '--hamiltonian' , dest = 'hamiltonian' , default = 'r' , help = 'Give the hamiltonian you want to solve: \'r\' is the red. bcs. Hamiltonian, \'f\' is the fac. int. Hamiltonian')
  parser.add_option('-n', '--inputname', dest = 'inputname' , default = None , help = 'If you don\'t use an inputfile you can only use some predefined model problems such as: r (rombouts:2010) , s (sambataro:2008), d (dang inputfile of Sn120 (it is called Sn120Neutrons))' )
  parser.add_option('-t', '--tdadict', dest = 'tdadict' , default = 'weak', help = 'Give the start_tda distribution if weak -> start tdadict of weak interacting regime, if strong -> start tdadict of strong interacting regime, else integer that contains the desired tdadict: 210401 -> {0:2 , 1:1 , 2:0 ,3:4 ,4:0, 5:1}')
  parser.add_option('-x','--start', dest = 'start' , default = -0.0001, help = 'Give the start value of the independend variable' , type = float)
  parser.add_option('-y','--step', dest = 'step' , default = -0.0001, help = 'Give the step of the independend variable', type = float)
  parser.add_option('-z','--end', dest = 'end' , default = -1., help = 'Give the end value of the independend variable', type = float)
  (options , args) = parser.parse_args(sys.argv[1:]) #args is list of positional arguments that remains after the processing of the parsing options
  return options.filename,  options.depvar,  options.interactionconstant, options.eta,options.levels, options.pairs, options.runstring,  options.hamiltonian,  options.inputname, options.tdadict, options.start, options.step , options.end, args

def main():
  """
  main function that solves the rg equations for general sp levels and degeneracies extracted from a given file 
  """
  rgw = True ; mov = False ; tdaf = False ; intm = True #rgw: writes Rg vars to outputfile , tdaf: for each step go back to xi = 0 and determine start tdasol , intm: write integrals of motion to outputfile, mov: if tdaf = True plot all xipaths and make a movie
  wd2 = -50. #(wd2 > 0 ) then we look at the correlation-energy in function of a changing density of the sp states. If the variable is not defined (wd < 0) then the number of sp levels stays constant
  filename, depvar, interactionconstant,eta , nlevel, npair, runstring, hamiltonian, inputname,tdadict,start , step , ende, args =  parse_commandline()

  if runstring == 'r': #handle the special case of a restart (this is special because we don't need to generate tdadicts or richardsoneq objects because all the start information is already contained in some outputfile)
    restart_somestates(depvar,step,ende,rgw,mov,tdaf ,afhvar = None ,linenr = -1 , intm = True , tw = 'a', exname = '')
    sys.exit(0) # we end the program normally
      
  #create RichardsonEq object
  if filename == None or runstring == 'f':
    rgeq = create_predefined_rgeq(interactionconstant,eta,npair,nlevel ,hamiltonian, inputname)
  else:
    rgeq = dr.ReaderOutput(filename, inputline= '', comment = '%').make_rgeq()

  #initialise the start tda dictionary 
  if tdadict == 'weak':    
    tdadict = tdadict_kleinekoppeling(rgeq.apair,rgeq.ontaardingen, rgeq.senioriteit)
  elif tdadict == 'strong':
    tdadict = {0 : rgeq.apair}
  else:
    tdadict = str2tdadict(tdadict)

  #creation of directory name where we save the output of the run
  name = "%srun_%s_p%gl%g" %(rgeq.__class__.__name__,runstring,rgeq.apair,rgeq.alevel) 
  if filename != None: name += filename
  else: name += 'inputname_' + inputname 
  if args: 
    name = args[0] + name 
  name.translate(None,'}{\\[];:.') 

  #generate seperate dir for the outputfiles if this is not a restart
  generate_dir(name,filename,None) #if filename is None nothing will be copied in the directory see the implementation of generate_dir in rgfunctions

  ###########################END OF INITIALISATION (start mainwork) ####################################################  
  if runstring == "f":
    kp =False # If the interaction constant is in the critical regime set kp = True
    #loop over file that contains the sp levels
    generating_data(rgeq,nlevel,filename,start,step,ende,wd2,tdadict,depvar,kp,namepf = 'plotenergygeo.dat' , rgwrite = rgw, intofmotion = intm)
    #generate some nice plots of the created data in the problem directory
    Plot_Geo_File('plotenergygeo.dat').standard_plot(rgw,intm) 
    
  elif runstring == "k":
    generating_datak(rgeq,tdadict,depvar,step,ende,rgwrite = rgw ,tdafilebool = tdaf,exname = '',moviede = mov, intofmotion = intm, printstep = 30 )       
    #generate some nice plots of the created data in the problem directory
    Plot_Data_File("plotenergy.dat").standard_plot(rgw , intm)
  
  elif runstring == 'e':
    allstatesgenerating_datak(rgeq,depvar,step,ende,rgw,mov,tdaf , intm = intm)
    #probeLowestExitedStates(rgeq,depvar,step,ende,rgw,mov,tdaf , intm = intm) 

  else:
    Exception('runstring: %s is not a valid runstring' %runstring)
   
def generating_data(rgeq,nlevel,infilename,afh,step,end,wd2,pairingdict,afhxas,kp = False,namepf = 'plotenergy.dat',checkfile =False,rgwrite = True, intofmotion = True , exname = ''):
  """
  function that generates the file with all the data 
  does the main job
  """
  saveopl = 123456.
  discstep = 1000#discontinuity step if two succeeding solutions differ by more then this we search the ground state from the beginning
  #create the data files
  plotenergyfile = open(namepf,"w")
  dr.info_1set(plotenergyfile,str(rgeq), exinfo = "#The variable we change is: %s\n#at (afh = %s)\tcE\tgE\tnig\td\tnlevels\trgvar(real)\trgvar(imag)\t ...\n" %(afhxas,afhxas),tdadict = pairingdict)
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
      energielev = waarden[1:nlevel+1]
      energielev.sort()
      nlevels,energielev = wd_processing(wd2,nlevel,rgeq.apair,energielev)
      nlevela,dega = checkdegeneratie(energielev,nlevels*[2]) #set degeneration standard to 2 of all the levels if you want higher degeneration of a level put the same energyvalue multiple times in the energyfile on the same line until the desired degeneracy (This function sets the degeneracy higher if energy values are close enough and treats them as one energylevel with increased degeneracy) 
      rgeq.energiel = array(energielev) ; rgeq.ontaardingen = np.array(dega); rgeq.senioriteit = np.zeros(nlevela)
      d = calculated(rgeq.energiel,rgeq.ontaardingen) 
      bb = calcnintgrondtoestand(rgeq)
      print "%f is de waarde van d bij afh = %f" %(d,defvar)
      """
      If the interaction constant is to small we can't find any solution with this pairing dictionary, therefore we call the generateRGsolKoppeling
      a function that gives back the RG energy of the groundstate and the RG variables for arbitrary interaction constants
      """      
      try:
        if rgeq.energy is None: #We just started so we need to find a startsolution but if kp is True we win some time by immediately raising a valueError this timewinst is significant for large systems
          rgeq = rg.RichardsonSolver(rgeq).main_solve(pairingdict,plotrgvarpath = False,xlim = None , ylim = None)          
          energierg = rgeq.get_energy()
        else:
          rgeq.energiel = array(energielev)
          energierg = rgeq.solve()
      except (ValueError, np.linalg.linalg.LinAlgError , rg.XiError) as e:
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
            rgeq = rg.RichardsonSolver(rgeq).main_solve(pairingdict)
            energierg = rgeq.get_energy()
          except:
            print '######################################################'
            print 'ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR'
            print ' ####################################################'
      #little discontinuity check
      if (rgeq.get_energy() + discstep < saveopl or rgeq.get_energy() - discstep > saveopl) and saveopl is not 123456.:
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
          rgeq = rg.RichardsonSolver(rgeq).main_solve(pairingdict)
          energierg = rgeq.get_energy()
      saveopl = energierg
      if energierg is not None:
        plotenergyfile.write("%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f" %(defvar, energierg-bb, energierg,bb , d,nlevel))
      if rgwrite is True:
        for i in range(rgeq.apair):
          plotenergyfile.write('\t%.12f' %( rgeq.rgsolutions[i].real ))
          plotenergyfile.write('\t%.12f' %( rgeq.rgsolutions[i].imag ))
      if intofmotion == True:
        integralsofm = rgeq.intofmotion()
        for i in range(rgeq.alevel):
          plotenergyfile.write('\t%.12f' %(integralsofm[i].real))
      plotenergyfile.write('\n')
      afh = afh + step
      if math.fabs(afh)+0.00001 > math.fabs(end) and math.fabs(afh)-0.00001 < math.fabs(end):
        break
  plotenergyfile.close()
  ifile.close()

def generating_datak(rgeq,pairingd,dvar,step,end ,xival = 1.,rgwrite = True,exname = '',moviede = False,tdafilebool = False , intofmotion = True,tw= "w" , conarray = [], printstep = 30):
  plotenergyfile = open("plotenergy%s.dat" %exname,tw) #opening of the file that's going to contain the results if you want to append data to existing file put howwrite = "a"
  olddir = os.getcwd() # remember the current directory, at the end of this routine we change back to the current dir so the final plots are contained in this dir
  if tdafilebool is True:
    tdafile = open('tdafile.dat','w')
    dirname = "moviedir%s" %exname
    rgf.generate_dir(dirname,None) #we use the function from rgfunctions.py
  d = rgf.calculated(rgeq.energiel,rgeq.ontaardingen)  #Calculation of the mean distance between the used sp levels
  bb = rgf.calcnintgrondtoestand(rgeq) 
  if tw == 'w': dr.info_1set(plotenergyfile,str(rgeq) , exinfo = "#The variable we change is: %s \n#d:  %f \n#and the noninteracting groundstate = %f\n#g\tcE\tgE\trgvar(real)\trgvar(imag)\t ...\n" %(dvar, d,bb),tdadict = pairingd)
  lastkp = False #boolean to know if the last g value circumvented a critical point
  #if rgeq.rgsolutions is None we haven't determined any rg variables so the first solution has to be determined from the corresponding tda solutions (xi = 0 -> xi = 1)
  #REMARK after the first solution we have a good guess for the next sol of the file so we don't need to start from tda but can directly
  #start from the previous solution if the stepwidth of the dependent variable of the file is low enough (WATCH OUT for critical points)
  if rgeq.rgsolutions is None and pairingd is not None: 
    rgeq = rg.RichardsonSolver(rgeq).main_solve(pairingd,xival = xival)   
    energierg = rgeq.get_energy()
  rgeqsaveback = None
  rgeq.solve()
  while (rgeq.getvar(dvar) != end ):  
    savergsolutions = np.copy(rgeq.rgsolutions) ; savedepvar = rgeq.getvar(dvar) #important if continuity check fails
    rgeq.setvar(dvar,savedepvar + step)
    if (rgeq.getvar(dvar)- abs(step)*0.8 < end and rgeq.getvar(dvar) + abs(step)*0.8 > end):
      rgeq.setvar(dvar,end)    
    if int(abs(savedepvar- end)/abs(step)) % printstep == 0 : print 'The variable we change: %s is: %s we end at: %f' %(dvar ,str(rgeq.getvar(dvar)), end) #print each 50 steps
    try:
      energierg = rgeq.solve()
      if not rgf.continuity_check(conarray, rgeq,crit = 1.8,dvar = dvar):
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
            print rgeq
            energierg,rgeq,rgeqsaveback = rgf.littleLoop(conarray[-2],step/2.,n*2,complexstepd = complexstep,end = send,dvar = dvar)
          elif abs(rgeqsaveback.getvar(dvar)-send)/abs(savestep)   > 4 and n < 100: 
            energierg,rgeq,rgeqsaveback = rgf.littleLoop(rgeq,step/2.,n*2,complexstepd = complexstep,end = send,dvar = dvar)
          else:
            if n < 310:
              send = savesend 
              complexstep = 10000
              step = savestep
            send += savestep/2.
            if n < 1000:
              step /= 2.
            if n > 1000:
              step *= 5.
            assert(isinstance(rgeqsaveback.g,complex))
            energierg,rgeq, rgeqsaveback = rgf.littleLoop(rgeqsaveback,step/2.,n*2,complexstepd = complexstep,end = send,dvar = dvar)
            if n > 1000:
              print 'n is starting to get large: %g  maybe you should consider changing the parameters that determine the circumvention of critical points' % int(n)
          
          rgeq.setvar(dvar, rgeq.getvar(dvar) + savestep - rgeq.getvar(dvar) % savestep) #To make sure that we solve the desired points of the independend variable
          rgeq.solve()
          if not rgf.continuity_check(conarray, rgeq,crit = 1.9 * n % 4 ,dvar = dvar):
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
          pairingdict  =  rgsolver.main_desolve(xistep = -0.01,rgwrite = False, plotrgvarpath = moviede,plotepath = False)
          tdasol = rgsolver.get_tda(None) 
          tdafile.write('%f\t%s\t%s\n'  %(rgeq.getvar(dvar),' '.join(map(str,pairingdict)), ' '.join(map(str,(tdasol)))))
        except (ValueError, np.linalg.linalg.LinAlgError,NameError , rg.XiError) as e:
          print 'problem in going back in xispace to xi = o to determine the corresponding tdadistribution'
          tdafile.write('#%f\t%s\n'  %(rgeq.getvar(dvar),'We couldn\'t find any solutions because their occured an error in desolving the Richardson-Gaudin solutions from XI =1 to XI = 0')) 
    plotenergyfile.write("%.12f\t%.12f\t%.12f" %(rgeq.getvar(dvar), energierg-bb, energierg)) 
    if (rgeq.getvar(dvar)- abs(step)*0.8 < end and rgeq.getvar(dvar) + abs(step)*0.8 > end):
      rgeq.setvar(dvar,end)    
    if rgwrite is True:
      for i in range(rgeq.apair):
        plotenergyfile.write('\t%.12f' %( rgeq.rgsolutions[i].real ))
        plotenergyfile.write('\t%.12f' %( rgeq.rgsolutions[i].imag ))
    if intofmotion == True:
      integralsofm = rgeq.intofmotion()
      for i in range(rgeq.alevel):
        plotenergyfile.write('\t%.12f' %(integralsofm[i].real))
    plotenergyfile.write('\n')
    lastkp = False   ; extremecp = False
    if rgeq.getvar(dvar) < end and np.sign(step) < 0 :
      rgeq.setvar(dvar,end)
    elif rgeq.getvar(dvar) > end and np.sign(step) > 0:
      rgeq.setvar(dvar,end)
  #end calculation underneath is just some cleaning up and closing files  
  plotenergyfile.close()  
  if tdafilebool is True:
    tdafile.close()
    if moviede is True:
      makemovie()
    os.chdir(olddir) #because we were in the movie dir but we want our other plots in the parentdir
  return energierg,rgeq  
    
def probeLowestExitedStates(rgeq,afhxas,step,ende,rgw,mov,tdaf , intm = True) :
  '''
  function that investigates excited states
  '''
  pairingd = tdadict_kleinekoppeling(rgeq.apair,rgeq.ontaardingen,rgeq.senioriteit)
  generating_datak(rgeq,pairingd,afhxas,step, ende,exname = 'groundtoestand',rgwrite = rgw, moviede =mov , tdafilebool = tdaf)
  nullevel =  list(pairingd)[-1] + 1 #nullevel is the first not occupied level in pairingd (remember tdadict_kleinekoppeling gives us a dict back that only contains the lowest  possibley occupied levels)
  for j in xrange(nullevel):
    pairingd[nullevel-j-1] -= 1
    for i in xrange(rgeq.alevel-nullevel):
      print 'we are calculating the following excitation: pair in level %g -> level %g ' %(j,nullevel+i)
      pairingd[nullevel+i] = 1
      extraname = str(nullevel-j-1) +'_'+str(nullevel +i)
      generating_datak(rgeq,pairingd,afhxas,step,ende,rgwrite = rgw,exname = extraname,moviede = mov,tdafilebool = tdaf)
      pairingd[nullevel+i] -= 1
    pairingd[nullevel-j-1] += 1.  
    #break #uncomment if you only want the lowest exited states
  Plot_Data_File().plot_spectrum(name = 'spectrum' , rgw = rgw , intm = intm)

def restart_somestates(afhxas,step,ende,rgw,mov,tdaf , afhvar = None , linenr = None , intm = True, tw = 'a' , exname = '' ):
  filecollection = File_Collector('.', r'\bplotenergy' , notsearch = r'\.swp|\.png' , notdir = 'xvysf')
  print filecollection.plotfiles
  #uncomment the three lines hereunder if you want a restart of a restart (need to implement the index where the run broke down manually in this case ('923'))
  #indexl = [filecollection.plotfiles.index(x) for x in filecollection.plotfiles if '245' in x] 
  #print indexl
  #filecollection.plotfiles = filecollection.plotfiles[indexl[0]+1:]
  cwd = os.getcwd()
  for file in filecollection.plotfiles:
    conarray = []
    for lnr in range(linenr-4,linenr+1):
      conarray.append(dr.get_restart_rgeq(file,afhvar = afhvar, linenr = lnr))
    print 'we start generating_datak with: ', file , conarray[-1]
    os.chdir(os.path.split(file)[0])
    generating_datak(conarray[-1].copy(),None,afhxas,step,ende ,rgwrite = rgw,exname = '',moviede = mov,tdafilebool = tdaf,tw = tw , conarray = conarray)
    os.chdir(cwd)
  Plot_Data_File().plot_spectrum(name = 'spectrum' , rgw = rgw , intm = intm)

def get_tdadictiteratorlowestgroup(npair, nlevel):
  """
  Remark quick hack, implement this for general degeneracies
  This is an example of a generator (as is the output from combinations_with_replacement) you can only run once with a for loop over a generator because it doesn't keep the values in memory it uses one value then forgets and go to next,... ; as opposed to iterators which keep all values in memory so you can run over an iterable object multiple times (list,string,file,...) (sometimes not so good if you have to many values ;) solution -> generators)
  """
  o = range(npair)
  v = range(npair , nlevel)
  for i in range(0,-1*npair,-1):
    if i == 0:
      yield o
    else:
      tdad = o[:i] + v[:-1*i]
      yield tdad

  for j in range(npair-1):
    oc = list(o)
    del(oc[j])
    yield oc + [npair]

def get_facgap(npair,nlevel):
  yield range(npair) 
  yield range(1,npair+1)

def allstatesgenerating_datak(rgeq,afhxas,step,ende,rgw,mov,tdaf , intm = True):
  '''
  Calculates all the eigenvalues of the pairingsHamiltonian at positive or negative interaction constant.
  If you want both write a wrapper or use the wrappermainsfile ;). 
  '''
  #tdacombinations = combinations_with_replacement(np.arange(rgeq.alevel),rgeq.apair)
  tdacombinations = get_facgap(rgeq.apair,rgeq.alevel)
  tdacor = open('tdacor.dat','w')
  tdacor.write('#This file contains the correspondence between the directories and the start tda distributions \n#The first column is the directory number and on the same line is the tda start distribution written \n')
  i = 0 #the i parameter makes it easy to restart after a failure just but i in the condition bigger as the i at failure
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
    if i == 0: exname = 'ground' #important for the plot of the spectrum
    else: exname = ''
    if i > -1: #parameter to restart easely after failure
      generate_dir('%g' %i,None,None) #only handy for a small amount of levels
      #tdastart needs to be a dictionary so we need to convert the list that contains one element of the permutation sequence to a dictionary    
      tdacor.write('%g\ttdadict= %s\n' %(i,' '.join(map(str,tdadict))))
      print 'we start generating_datak with: ', tdastartd
      generating_datak(rgeq,tdastartd,afhxas,step,ende ,rgwrite = rgw,exname = exname,moviede = mov,tdafilebool = tdaf)
      os.chdir(os.path.abspath(os.path.join(os.getcwd(), os.path.pardir)))
    i += 1
    print i
  tdacor.close()
  Plot_Data_File().plot_spectrum(name = 'spectrum' , rgw = rgw , intm = intm)
    
def romboutsprob():
  #definition of problem of artikel stefan:2010 (endg is typically -0.075 , and apair = 10)
  eendlev = array([0.04,0.08,0.16,0.20,0.32,0.36,0.40,0.52,0.64,0.68,0.72,0.8,1])
  eendlev = np.sqrt(eendlev)
  seniority = zeros(len(eendlev),float)
  degeneration = [4,4,4,8,4,4,8,8,4,8,4,8,12]
  return eendlev , degeneration ,seniority 
  
def picketfence(alevel):
  #little toy problem introduced by sambataro:2008
  eendlev = np.arange(1,alevel+1)
  seniority = zeros(alevel,float)
  degeneration = np.ones(alevel)*2
  return  eendlev , degeneration ,seniority 

def dang(filename = 'Sn120Neutrons' , cutoff = 1e5, apair = 35):
  #typical number of pairs of Sn is 35 interacting  neutron pairs , typical interaction constant is -0.137
  eendlev ,degeneration = dangSn(filename,cutoff)
  nlevel,degeneration = checkdegeneratie(eendlev ,degeneration,nauw = 0.0011)
  seniority = [0]*len(eendlev )  
  eendlev ,degeneration,senioriteit,npair,nlevel,fermil,extrae = windowcreation(eendlev ,degeneration,senioriteit,apair,5,5)
  return  eendlev , degeneration ,seniority , apair
  
def facintmain():
  #definition of problem of artikel stefan:2010 
  energy , deg , sen = romboutsprob()
  rgeq = rg.RichFacInt(energy , deg, sen, -0.075 , 1., 10)
  tdastartd = {0:10,1:0,2:0}
  enddatak = -0.0002
  stepg = 0.0001
  dvar = 'g'
  generate_dir('stefantda',None,None)
  generating_datak(rgeq,tdastartd,dvar,stepg,enddatak ,rgwrite = True,exname = '',moviede =True,tdafilebool = True)
  Plot_Data_File("plotenergy.dat").standard_plot(True , True)

def allstatesoneg(npair = 3):
  energy , deg ,sen = romboutsprob()
  energy , deg , sen = picketfence(12) 
  rgeq = rg.RichFacInt(energy , deg, sen, -0.075 , 1., 10)
  rgeq = rg.RichRedBcs(energy , deg ,sen , -1, 6)
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
  plotterxi = Plot_Xi_File(rgeq.g)
  plotterxi.procesfiles('.', 'xipath')    
  plotterxi.plot_spectrumxichange()

def testcircumvent():
  '''
  test for the rgsolutions when xi is not one and/ or g is imaginary
  '''
  energy , deg , sen = picketfence(alevel = 12) 
  g = -1.0001
  rgeq = rg.RichFacInt(energy , deg ,sen , g,1., 6)
  tdastartd = {0:6}
  #tdastartd = tdadict_kleinekoppeling(rgeq.apair,rgeq.ontaardingen, rgeq.senioriteit)
  enddatak = -0.0001
  stepg = 0.001
  rgvar = None
  afhxas = 'g'
  #generate seperate dir for the solved problem
  generate_dir('testcircumventlowxifacintrest2',None,None)
  print tdastartd
  xival =1.
  while xival >= 0.:
    try:
      rgeq2 = rg.RichardsonSolver(rgeq.copy()).main_solve(tdastartd,xistep = 0.01,xival=xival)   
    except rg.XiError:
      pass
    energierg = rgeq2.get_energy()
    while rgeq2.g.imag < abs(rgeq2.g.real/10.):
      rgeq2.g += 1j*0.001
      energierg = rgeq2.solve()
      print rgeq2.g 
    generating_datak(rgeq2,tdastartd,afhxas,stepg,enddatak ,rgwrite = True,exname = '%f' %xival,moviede = False,tdafilebool = False,xival = xival)
    Plot_Data_File("plotenergy%f.dat" %xival).standard_plot(True , True)
    xival -= 0.1  

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

def stijnijzer():
  rgeq = dr.get_restart_rgeq('plotenergy1.dat', afhvar = -0.9233, linenr = None, startrg = 3)
  print rgeq
  exname = 'nauwkeurigomgekeerd'
  generating_datak(rgeq,None, 'g' ,-0.000001 , -9.23435 , tdafilebool = True , exname = exname)
  Plot_Data_File("plotenergy%s.dat" %exname).standard_plot(True , True)

def behaviour_conpoint():
  reader = dr.ReaderOutput('plotenergy.dat')
  conp = reader.npair #condensed pairs -> npair == Moore-Read point
  nconp = reader.npair - conp #uncondensed pairs
  conpoint = reader.eta/(2.*nconp +conp -1-2.*np.sum(np.array(reader.degeneracies)/4. -np.array(reader.seniorities)))
  reader.readrgvars(afhvar = conpoint + 0.005, linenr = None, startrg = 3)
  rgeq = reader.make_rgeq(types = 'RichFacInt', depvar = None)
  print rgeq
  exname = 'nauwkeurigconpoint'
  generating_datak(rgeq,None, 'g' ,-0.00001 ,conpoint -0.005 , tdafilebool = False , exname = exname)
  Plot_Data_File("plotenergy%s.dat" %exname).standard_plot(True , True)

def dirrunreadgreen():
  for name in os.listdir('.'):
    if os.path.isdir(name):
      os.chdir(name)
      print "changed to dir: %s " %name
      mainreadgreen()
      os.chdir('../..')

def mainreadgreen():
  plotter = Plot_Data_File().plot_spectrum(xlim = (-0.6,0), ylim = None,search = 'plotenergy', rgw = True, intm = True, name = 'specreadgreen' , readgreen = True , standard = False)
  file =  open('readgreenfile.dat', 'r')
  for line in file:
    if 'filename' in line:
      match = re.search(r'/(\d+)/', line)
      print os.getcwd()
      os.chdir(str(match.group(1)) )
      #readgreentda()
  file.close()
  
def readgreentda(): 
  reader = dr.ReaderOutput('plotenergy.dat')
  readgreenpoint = reader.eta/(2.*(reader.npair-1) -2.*np.sum(np.array(reader.degeneracies)/4. -np.array(reader.seniorities)))
  reader.readrgvars(afhvar = readgreenpoint, linenr = None, startrg = 3)
  rgeq = reader.make_rgeq(types = 'RichFacInt', depvar = None)
  print rgeq
  if readgreenpoint < rgeq.g :
    generating_datak(rgeq, None, 'g' ,-0.0002 , readgreenpoint  , tdafilebool = True, exname = 'readgreen', moviede = True)
  elif readgreenpoint > rgeq.g :
    generating_datak(rgeq, None, 'g' ,0.0002 , readgreenpoint , tdafilebool = True, exname = 'readgreen' , moviede = True)
  else:   
    rgsolver = rg.RichardsonSolver(rgeq)
    pairingdict  =  rgsolver.main_desolve(xistep = -0.01,rgwrite = True, plotrgvarpath = True,plotepath = False)
    print pairingdict
    Plot_Xi_File('.',"xipath").plotrgvarsxi(name = 'rgvxi' ,xlim = None , ylim = None)

def stijnd():
  readdata = dr.ReaderInp('pairing-parameters.inp', comment = '*')
  rgeq = readdata.make_rgeq()
  tdadict = tdadict_kleinekoppeling(rgeq.apair,rgeq.ontaardingen, rgeq.senioriteit)
  rgeq.g = -0.0001
  #tdadict = {0:rgeq.apair}
  exname = 'stijnconf3'
  #rgeq = dr.get_restart_rgeq('plotenergystijnconf2.dat', afhvar = -0.319600, linenr = None)
  generating_datak(rgeq, tdadict , 'g' ,-0.0005 , -0.050100	, tdafilebool = False, exname = exname)
  Plot_Data_File("plotenergy%s.dat" %exname).standard_plot(True , True)

def test_critical():
  for i in range(1,8):
    elev = [1,2] ; ont = [2*i,2] ; sen = [0,0] ; g = -0.0001
    rgeq = rg.RichRedBcs(elev,ont,sen,g,1+i)
    tdadict =tdadict_kleinekoppeling(rgeq.apair, rgeq.ontaardingen, rgeq.senioriteit)
    dir = '%gcorner' %(i+1)
    #os.mkdir(dir)
    os.chdir(dir)
    generating_datak(rgeq,tdadict,'g',-0.001, -1.1,exname = '%g' %(i+1))
    os.chdir('..')

if __name__ == "__main__":
  main()
  #dirrunreadgreen()
  #behaviour_conpoint()
  #parse_commandline()
  #dangmain()
  #testcircumvent()
  #addlevel() #function in rgfunctions
  #facintmain()
  #allstatesoneg()
  #stijnd()
  #stijnijzer()
  #test_critical()
