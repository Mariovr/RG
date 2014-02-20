# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Mario Van Raemdonck, 2014;
# (c) Ghent University, 2014
# -*- coding: utf-8 -*-
#! /usr/bin/env python 

## module powell
from numpy import identity,array,dot,zeros,argmax , asarray
from math import sqrt , log
from scipy.optimize import fmin_bfgs , fmin_cg

def bracket(f,x1,h):
  """ 
  Finds the brackets (a,b) of a minimum of the energy of self.ham
  The search starts downhill from x1 with a step length h.
  """
  c = 1.61803398875
  f1 = f(x1)
  x2 = x1 + h; f2 = f(x2)
  #Determine downhill direction and change sign of h if needed
  if f2 > f1:
    h = -h
    x2 = x1 + h; f2 = f(x2)
  #Check if minimum between x1 - h and x1 + h
    if f2 > f1: return x2,x1 - h
  #Search loop
  for i in range (100):
    h = c*h
    x3 = x2 + h; f3 = f(x3)
    if f3 > f2: return x1,x3
    x1 = x2; x2 = x3
    f1 = f2; f2 = f3
  print 'Bracket did not find a mimimum'

def search(f,a,b,tol=1.0e-7):
  """      
  (see chapter 10 in Kiusalaas_J._Numerical_Methods_in_Engineering_with_Python_(2005)(en)(424s).pdf)
  Golden section method for determining x that minimizes the energy of self.ham
  The minimum has to be in (a,b).
  """
  nIter = -2.078087*log(tol/abs(b-a)) 
  R = 0.61803398875
  C = 1.0 - R
  # First telescoping
  x1 = R*a + C*b; x2 = C*a + R*b
  f1 = f(x1); f2 = f(x2)
  # Main loop
  for i in range(int(nIter)):
    if f1 > f2:
      a = x1
      x1 = x2; f1 = f2
      x2 = C*a + R*b; f2 = f(x2)
    else:
      b = x2
      x2 = x1; f2 = f1
      x1 = R*a + C*b; f1 = f(x1)
  if f1 < f2: return x1,f1
  else: return x2,f2

def powell(F,x,h=0.1,tol=1.0e-6):
  """xMin,nCyc = powell(F,x,h=0.1,tol=1.0e-6)
  Powells method of minimizing F(x).
  x= starting point, h= initial search increment used in bracket
  xMin = mimimum point, nCyc = number of cycles
  """
  def f(s): return F(x + s*v) # F in direction of v
  n = len(x) # Number of design variables
  df = zeros((n)) # Decreases of F stored here
  u = identity(n)*1.0 # Vectors v stored here by rows
  for j in range(30): # Allow for 30 cycles:
    xOld = x.copy()
    # Save starting point
    fOld = F(xOld)
    # First n line searches record decreases of F
    for i in range(n):
      v = u[i]
      a,b = bracket(f,0.0,h)
      s,fMin = search(f,a,b)
      df[i] = fOld - fMin
      fOld = fMin
      x = x + s*v
    # Last line search in the cycle
    v = x - xOld
    a,b = bracket(f,0.0,h)
    s,fLast = search(f,a,b)
    x = x + s*v
    # Check for convergence
    if sqrt(dot(x-xOld,x-xOld)/n) < tol: return x,j+1
    iMax = int(argmax(df))
    for i in range(iMax,n-1):
      u[i] = u[i+1]
      u[n-1] = v
  print 'Powell did not converge'

def f(x, *args):
    u, v = x
    a, b, c, d, e, f = args
    return a*u**2 + b*u*v + c*v**2 + d*u + e*v + f
def gradf(x, *args):
    u, v = x
    a, b, c, d, e, f = args
    gu = 2*a*u + b*v + d     # u-component of the gradient
    gv = b*u + 2*c*v + e     # v-component of the gradient
    return asarray((gu, gv))

def main2():
  """test for the scipy optimization functions"""
  args = (2, 3, 7, 8, 9, 10)  # parameter values
  x0 = asarray((0, 0))  # Initial guess.
  from scipy import optimize
  res1 = optimize.fmin_cg(f, x0,  args=args)
  #res1 = optimize.fmin_bfgs(f, x0,  args=args)
  print 'res1 = ', res1

def main(*args , **kwargs):
  def F(x): return 100*(x[1] - x[0]**2)**2+ (1-x[0])**2
  xstart = array([-1,1])
  print F(xstart)
  xmin , niter = powell(F , xstart)
  print xmin , niter
  print F(xmin)

if __name__ == "__main__":
  main2()

