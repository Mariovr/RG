#! /bin/env python
from string import maketrans

import richardsongaudin as rg

class Reader(object):
  """
  class to read in data files from stijns program
  """
  def __init__(self, filename , comment = '#' , transf = "dD" , transa = "ee"):

    self.file = open(filename , 'r')
    self.comment = comment
    self.nlevel = None
    self.npair = None
    self.g = None
    self.elevels = []
    self.degeneracies = []
    self.senioritys = []
    self.trans = maketrans(transf,transa)
    self.read()

  def read(self): 
    for line in self.file:
      if line[0] == self.comment:
        continue
      data = line.split()
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
        self.senioritys.append(int(data[3]))

  def maketranslationtable(self,first = "dD",after = "ee"): 
    return maketrans(first, after)

  def make_rgeq(self , types = 'RichRedBcs' , eta =None):
    if types == 'RichRedBcs':
      return rg.RichRedBcs(self.elevels, self.degeneracies , self.senioritys, self.g, self.npair)
    if types == 'RichFacInt' :
      return rg.RichFacInt(self.elevels, self.degeneracies, self.senioritys,self.g,
          eta,self.npair)


def main():
  d = Reader('pairing-parameters.inp', comment = '*')
  rgeq = d.make_rgeq()
  print str(rgeq)

if __name__ == "__main__":
  main()
