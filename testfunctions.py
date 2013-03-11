import numpy as np
from numpy import array, zeros, ones, arange
from numpy import append, delete, polynomial

def generateRGsolKoppelingcomplex(rgeq,end,d,pairingd = None,step = 1.,begin = None,klein = True,gdata = True):
  """
  function that generates the RG solution of a particular system for the interaction constant "interactionc" , this because the interaction
  constant is not compatible with "apair 0 0 0 ... ". we determine what is the closest d-interactionc or interactionc to determine
  if we start from a strong interaction constant or a small. If there happens to be a critical point we circumvent it by
  making the interaction constant complex (using litteLoop from rgfunctions.py)
  -------------------------------------------------------------------INPUT--------------------------------
  interactionc = the interaction constant
  npair = the amount of pairs that are contained in the system
  nlevel = the amount of single particle levels
  wd2 = if neg, not used, if pos, we only use the sp levels in the interval groundstate+wd2
  elev = list of sp levels
  startI = start interaction constant that goes down until the interaction constant were we are searching the solution's for is reached, standard set on -0.1
  step = the step that determines the variation of the interaction constant; standard set on 1
  pairingd(start tda distribution) (OPTIONAL argument)
  ---------------OUTPUT------------------------------------------------------------------------
  the groundstate energy and the Richardson-Gaudin variables at the given interaction constant
  """  
  #if rgvar is None we haven't determined any rg variables so the first solution of the file has to be determined from the corresponding tda solutions (xi = 0 -> xi = 1)
  #REMARK after the first solution we have a good guess for the next sol of the file so we don't need to start from tda but can directly
  #start from the previous solution if the stepwidth of the dependent variable of the file is low enough (WATCH OUT for critical points)
  if rgeq.rgsolutions is None:
    rgeq = genstartsol(rgeq,d,end,step,begin = begin,pairingd = pairingd, gdata = gdata, klein = klein)
  assert(rgeq.rgsolutions is not None)  
  """
  we arrived at a critical point so now we circumvent it by making g complex.
  """
  goodsol = False
  i = 0
  while goodsol == False:
    cs = 100000
    try:
      savergeq = copy.deepcopy(rgeq)
      energierg,rgeq = getRGatg(rgeq,pairingd,end,step/4. ,complexstepd = cs)
      goodsol = True
    except (ValueError, np.linalg.linalg.LinAlgError) as e:
      print 'nosolution so we make the step in complex space and real space smaller'
      rgeq = savergeq
      cs *= 20
      i+= 1
      if i == 2:
        print'In generateRGsolKoppelingcomplex fatal error we can\'t reach the end interaction constant even after we made the step of the interaction constant in complex and real space 100 times smaller'
        sys.exit(1)
               
  print 'The interaction constant is: %s ' % str(rgeq.g) 
  return energierg,rgvar




def testhermite(n):
  """
  test contains the coefficients of the hermite polynomial that gives us the roots.
  the length is (n+1) because there exists also a zero Hermite function who has no roots, we want the n roots
  of the n^e Hermite function, zo the length of the coefficient matrix need to be n+1
  OKE the graph generated with this functions corresponds to the graph in the personal notes of S. De Baerdemacker
  """
  test = zeros(n+1,float)
  test[n] = 1
  pl.figure(1)
  for i in xrange(n):
    rootsi = polynomial.hermite.hermroots(test)
    print rootsi
    for j in xrange(len(rootsi)):
      pl.plot(n-i,rootsi[j],'go',markersize=5 )
    test = delete(test,[0],None)
  pl.show()
  pl.savefig("hermitepol.png")

def fill_hermitematrix(n):
  hermitematrix = zeros((n,n),float)
  hermitematrix[arange(n-1)+1,arange(n-1)] = [math.sqrt(2*(a+1)) for a in range(n-1)]
  hermitematrix[arange(n-1),arange(n-1)+1] = [math.sqrt(2*(a+1)) for a in range(n-1)]
  return hermitematrix  
  
def main2():
  """
  test main for hermite polynomials
  """
  testhermite(20)   
  
if __init__ == '__main__':
  main2()
  
  
