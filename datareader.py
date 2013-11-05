#! /bin/env python
from string import maketrans
import re
from numpy import empty

import richardsongaudin as rg

def get_restart_rgeq(filename, afhvar = None, linenr = None , startrg = 3 , types = None , depvar = 'xi'):
  #remark last two options of this function are not necessary, but they handle old outputfiles that didn't contain type and independend variable information in the header
  outputreader = ReaderOutput(filename, comment = '#')
  outputreader.readrgvars(afhvar = afhvar, linenr = linenr , startrg = startrg) #if ahvar is not None we read rgvars when depvar[1] == afhvar otherwise we read rgvars at linenumber = linenr, remark xifile has startrg = 2, plotgeofile has startrg = 6 , plotenergyfile has startrg = 3 
  return outputreader.make_rgeq(types = types , depvar = depvar)
  
def matchlist(line):
  match = re.search(r'\[(?P<list>[+\-\d\s.eE]*)\]' , line) 
  if match:
    return map(float,match.group('list').split())
  else:
    Exception("No match for a list was found in the following line: %s", line)
    
def matchdict(line):
  match = re.search(r'\{(?P<dict>[+\-\d\s.eE:,]*)\}' , line) 
  if match:
    return dict(map(lambda x: tuple(map(lambda x : int(float(x)),x.split(':'))) ,match.group('dict').split(','))) #little trick at second lambda function to make strings '3.0' transformable to integer 3
  else:
    Exception("No match for a list was found in the following line: %s", line)

def info_1set(filehandler,rgeqstring, exinfo = '#',tdadict = None,contenergy = None):
  """
  Writes header to outputfile very important for datareader because it uses this to read the info back in.
  """
  filehandler.write(rgeqstring)
  filehandler.write('%s' %exinfo)
  if tdadict is not None:
    filehandler.write('#the tda start distribution is: %s\n' %(str(tdadict)))
  if contenergy is not None:
    filehandler.write('#The energy from the electrons that don\'t participate in the pairing is %f\n' % contenergy)

def make_newstyle(file , comment = '#', rgindex = 3):
  """
  Function that changes the order of the Richardson-Gaudin variables instead of the old style reval reval reval ... Im imval imval imval ...
  to the new style reval imval reval imval .....
  """
  filehandler = open(file , 'r')
  savenew = open('newstyle' + file , 'w')
  for line in filehandler:
    if line[0] == comment:
      savenew.write(line)
      continue
    data = line.split('\t')
    del(data[rgindex+1]) #delete 'Im'
    rerg = data[rgindex].split()
    imrg =data[rgindex+1].split()
    newrgstr = data[0:rgindex] + [re + '\t' + im for re, im in zip(rerg , imrg)]
    savenew.write('\t'.join(newrgstr) + '\n')
  filehandler.close()
  savenew.close()

class Reader(object):
  """
  Abstract data class that reads in data from files
  """
  def __init__(self, filename , comment = '#'):
    self.filename = filename
    self.comment = comment
    self.nlevel , self.npair , self.g , self.rgvars, self.eta = None , None , None, None, None
    self.elevels, self.degeneracies , self.seniorities = [],  [] , []
    self.depvar = {'depvar' : "",'depval': None} #This attribute contains a string, that defines the independend variable that is chanced under the header (only relevant for outputfiles, set by read (if it apply's )) , and the second element is the value of the independend variable where the rgvars are read in (set by readrgvars)
    self.tdadist = None
    self.kind = None
    self.read()

  def read(self):
    """
    Function that reads the variables that define the problem at hand
    """
    raise Exception("Implement this function in a subclass from Reader")

  def readrgvars(self):
    """
    Function that reads the corresponding solutions to the system defined by the variables obtained by self.read
    It should set self.rgvars and update self.depvar['depval']
    """
    raise Exception("Implement this function in a subclass from Reader")
  
  def make_rgeq(self , types = None , depvar = 'g'):
    """
    returns a RichardsonEq object, if the rgvars are read in with xi = 1 otherwise with low xi and no energy
    """
    if types == None: types = self.kind
    if types == 'RichRedBcs':
      rgeq = rg.RichRedBcs(self.elevels, self.degeneracies , self.seniorities, self.g, self.npair , rgsol = self.rgvars)
    elif types == 'RichFacInt' :
      rgeq = rg.RichFacInt(self.elevels, self.degeneracies, self.seniorities,self.g,self.eta,self.npair, rgsol = self.rgvars)
    else:
      raise Exception("%s is not a subclass of RichardsonEq (see the file: richardsongaudin.py)" %types )
    if rgeq.rgsolutions != None:
      rgeq.xi = 1.
      try: rgeq.setvar(self.depvar['depvar'],self.depvar['depval'])
      except : rgeq.setvar(depvar,self.depvar['depval'])
      rgeq.solve() #sets the energy

    return rgeq

class ReaderOutput(Reader):
  def read(self):
    """
    Reads in the header of the outputfiles of generating_datak , generating_data of which the name is typically: "plotenergy.dat"
    """
    with open(self.filename, 'r') as file:
      for line in file:
        if line.startswith(self.comment):
          if 'RichardsonEq' in line:
            try: self.kind = re.search(r':\s*(\w+)',line).group(1)
            except: print 'a kind was not found probably old outputformat you can fix it by submitting an extra parameter -> make_rgeq(types = RichRedBcs)'
          elif 'energylevels' in line:
            self.elevels = matchlist(line)
          elif 'degeneracies' in line or 'degeneracys' in line:
            self.degeneracies = matchlist(line)
          elif 'seniorities' in line or 'senioritys' in line:
            self.seniorities = matchlist(line)
          elif 'interaction' in line:
            self.g = float(re.search(r':\s*([\-+.\deE]+)' , line).group(1))
          elif 'pairs' in line:
            self.npair = int(re.search(r'pairs[:=\s]*(\d+)',line).group(1))
          elif 'change' in line:
            try: self.depvar['depvar'] = re.search(r'change[\s\w]*:\s*(\w+)',line).group(1)
            except: print 'a indep var. was not found probably old outputformat you can fix it by manually giving a parameter to make_rgeq(depvar = \'xi\')'
          elif 'eta' in line:
            match = re.search(r':*\s*([+\-\d.]+[eE+\-\d]*)', line) 
            if match: 
              self.eta = float(match.group(1))
          elif 'tda' in line:
            self.tdadist = matchdict(line) 
        else:
          break #header is finished
    self.nlevel = len(self.elevels)
    print 'we have read the information from the header of ', self.filename

  def readrgvars(self,afhvar = None , linenr = None , startrg = 3):
    """
    Function that reads back in the rgvars out of an plotenergy.dat file generated by generating_datak.
    If linenr is given it reads the rgvars at linenr in the file with name: name
    If linenr is not given it reads the rgvars when the float in the first column = afhvar
    If the file is from the old output format of generating_datak (re re re ... Im im im im ...) use the file converter -> make_newstyle and then use this function to read the rgvariables
    """
    with open(self.filename,'r') as ref: #neat trick in python closes the file automatically when you leave the with block
      if linenr != None:
        self.rgvars, self.depvar['depval'] = self.get_rgvars(ref.readlines()[linenr].split(), startrg)
      else:
        for line in ref: #make sure afhvar is not None if linenr is None
          if line[0]== '#':
            continue
          data = line.split()
          try:
            float(data[0])
          except:
            continue
          if float(data[0]) == afhvar:
            self.rgvars, self.depvar['depval'] = self.get_rgvars(data , startrg)
    print 'we have read the Richardson-Gaudin variables:', self.rgvars , ' at :', self.depvar

  def get_rgvars(self,datalist , start):
    mixdata = map(float,datalist)[start:]
    redata = mixdata[:2*self.npair:2]
    imdata = mixdata[1:2*self.npair:2]
    assert(len(redata) == len(imdata))
    rgvars = empty(len(redata),complex)
    rgvars.real = redata
    rgvars.imag = imdata
    return rgvars , float(datalist[0])

class ReaderInp(Reader):
  """
  class to read in data files from stijns program
  """
  def __init__(self, filename , comment = '#' , transf = "dD" , transa = "ee"):
    self.trans = maketrans(transf,transa)
    super(ReaderInp,self).__init__(filename,comment)

  def read(self): 
    with open(self.filename,"r") as file:
      for line in file:
        if line[0] == self.comment:
          continue
        data = line.split()
        if len(data) == 0: break
        if data[0] == 'nlev':
          self.nlevel = int(data[1])
        elif data[0] == 'Npair':
          self.npair = int(data[1])
        elif data[0] == 'G':
          data = data[1].translate(self.trans)
          self.g = float(data)
        elif data[0].startswith('lev'):
          self.elevels.append(float(data[1].translate(self.trans)))
          self.degeneracies.append(int(data[2]))
          self.seniorities.append(int(data[3]))

  def maketranslationtable(self,first = "dD",after = "ee"): 
    return maketrans(first, after)

def main():
  d = ReaderInp('pairing-parameters.inp', comment = '*')
  print  d.make_rgeq(types = 'RichRedBcs')
  d = ReaderOutput('xipath-0.039095.dat')
  print d.make_rgeq(types = 'RichRedBcs') #is not the type of the solutions but no problem because we haven't read in the solutions (we need to give it manually)
  d.readrgvars(linenr = -6, startrg = 2) #reading of solutions: remark xifile has startrg = 2, plotgeofile has startrg = 6 , plotenergyfile has startrg = 3 
  rgeq = d.make_rgeq(types = 'RichFacInt', depvar = 'xi')
  print str(rgeq)

if __name__ == "__main__":
  main()
  #make_newstyle( 'Dang120neutronwin5_5sen2.dat', comment = '#', rgindex = 3)
