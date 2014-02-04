import os
import writepairing as wp
import rgfunctions as rgf

"""
In this file I collect my wrapperfunctions for the main program in writepairing.py.
Typically these are functions that contain a loop over some input variables which are getting 
feed to the mainprogram in writepairing.py by calling os.system(python path_to_writepairing.py input_variables)
"""

def mainpairchange():
  nlevel = 8
  name = '%gp' %nlevel 
  rgf.generate_dir(name,None,None) #if filename is None nothing will be copied in the directory see the implementation of generate_dir in rgfunctions
  for ap in range(1,nlevel/2 +1):
    cmd = 'python ../../RG/writepairing.py -n s -r e -p %g -l %g -H f -x -0.001 -y -0.001 -z -0.6' % (ap, nlevel)
    print 'we execute %s' %cmd
    os.system(cmd)
  wp.dirrunreadgreen() #if an analyses of the read-green line is wanted
    
def main():
  name = 'balktosquaredifferentgfacint1'
  fname = 'vierkantrechthoekconstantV'
  rgf.generate_dir(name,fname,None)
  for interaction in  [-10000]:
    cmd = 'python ../writepairing.py -r f -f %s -H r -i %f' %(fname , interaction )
    print 'we execute %s' %cmd
    os.system(cmd)

def gapsurvey():
  """
  Make sure that the generator that determines the tdastartdistributions for the allstatesgenerating_datak function is set to the generator: get_facgap() because standard the generator is the generator which generates all states.
  """
  scale = 4 #multiplication factor for the sp levels if the number of pairs is given
  for ap in [1,2,3,4,5,10,15,20,25,30,50,75,100]:
    cmd = 'python writepairing.py -n s -r e -p %g -l %g -H f -x -0.001 -y -0.0001 -z -1. exname_%g' % (ap, scale * ap, ap)
    print 'we execute %s' %cmd
    os.system(cmd)
  
if __name__ == "__main__":
  #main()
  #mainpairchange()
  gapsurvey()

