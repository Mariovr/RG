import math
import os

from rgfunctions import *

def main():
  name = 'groundstateallpairsrombouts2010' 
  generate_dir(name,None,None) #if filename is None nothing will be copied in the directory see the implementation of generate_dir in rgfunctions
  for ap in range(1,40):
    cmd = 'python ../writepairing.py -n r -r k -p %g' % ap
    print 'we execute %s' %cmd
    os.system(cmd)
    

if __name__ == "__main__":
  main()
