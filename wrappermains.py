import math
import os

from rgfunctions import *

def mainpairchange():
  name = 'groundstateallpairsrombouts2010' 
  generate_dir(name,None,None) #if filename is None nothing will be copied in the directory see the implementation of generate_dir in rgfunctions
  for ap in range(1,40):
    cmd = 'python ../writepairing.py -n r -r k -p %g' % ap
    print 'we execute %s' %cmd
    os.system(cmd)
    
def main():
  name = 'balktosquaredifferentgfacint1'
  fname = 'vierkantrechthoekconstantV'
  generate_dir(name,fname,None)
  for interaction in  [-10000]:
    cmd = 'python ../writepairing.py -r f -f %s -H r -i %f' %(fname , interaction )
    print 'we execute %s' %cmd
    os.system(cmd)

if __name__ == "__main__":
  main()

