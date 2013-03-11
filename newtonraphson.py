# -*- coding: utf-8 -*-
from numpy import zeros,dot,sqrt
from numpy import *
import sys
import math 
from scipy import linalg

def swapRows(v,i,j):
  if len(v) == 1:
    v[i],v[j] = v[j],v[i]
  else:
    temp = v[i].copy()
    v[i] = v[j]
    v[j] = temp
def swapCols(v,i,j):
  temp = v[:,j].copy()
  v[:,j] = v[:,i]
  v[:,i] = temp


def err(string):
  print string
  raw_input("Press return to exit")
  sys.exit()


def gaussPivot(a,b,tol=1.0e-9):
  """ x = gaussPivot(a,b,tol=1.0e-9).
  Solves [a]{ x} = { b} by Gauss elimination with
  scaled row pivoting
  """
  n = len(b)
  # Set up scale factors
  s = zeros((n))
  for i in xrange(n):
    s[i] = max(abs(a[i,:]))
  for k in xrange(0,n-1):
    # Row interchange, if needed
    p = int(argmax(abs(a[k:n,k])/s[k:n])) + k
    if abs(a[p,k]) < tol:
      error.err("Matrix is singular")
    if p != k:
      swapRows(b,k,p)
      swapRows(s,k,p)
      swapRows(a,k,p)
    # Elimination
    for i in xrange(k+1,n):
      if a[i,k] != 0.0:
	lam = a[i,k]/a[k,k]
	a[i,k+1:n] = a [i,k+1:n] - lam*a[k,k+1:n]
	b[i] = b[i] - lam*b[k]
  if abs(a[n-1,n-1]) < tol:
    error.err("Matrix is singular")
  # Back substitution
  for k in xrange(n-1,-1,-1):
    b[k] = (b[k] - dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
  return b

def solve(f,x , tol=1.0e-10):
  """
  newton-raphson method for a system of nonlineair equations.
  f is the system
  x is the first guess
  tol is the accuracy of the solution ,standard tol = 1.0e-8
  """
  def jacobian(f,x ):
    """
    if our callabal function class f has a jacobian method we use this as the jacobian
    otherwise we approximate the jacobian by a simple formule for the derivative
    """
    if hasattr(f,"jacobian"):
      return f.jacobian(x),f(x)
    h = 1.0e-8
    n = len(x)
    jac = zeros((n,n),complex)
    f0 = f(x)
    for i in xrange(n):
      temp = x[i]
      x[i] = temp + h
      f1 = f(x)
      x[i] = temp
      jac[:,i] = (f1 - f0)/h
    return jac,f0   
  for i in xrange(200):
    jac,f0 = jacobian(f,x)
    nauw = abs(f0)
    if sum(nauw) < tol: return x
    #dx = gaussPivot(jac,-f0) #we commented our self written linalg solver out becaus the scipy solver is quicker and does the same thing
    dx = linalg.solve(jac,-f0)
    x = x + dx
    #print "iteration: %g" %i , "|f(x)|= " , nauw , x 
  print "Too many iterations"
  return x
  
def main():
  def test_system(x):
    f = zeros((len(x)),float)
    f[0] = sin(x[0]) + x[1]**2 + log(x[2]) - 7.0
    f[1] = 3.0*x[0] + 2.0**x[1] - x[2]**3 + 1.0
    f[2] = x[0] + x[1] + x[2] - 5.0
    return f
  guess = array([1., 1.,1.])
  print solve(test_system,guess)
  #solution of this test system is : [ 0.59905376  2.3959314   2.00501484]
  
if __name__ == "__main__":
  main()