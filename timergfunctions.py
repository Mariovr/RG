import sys
import os , unittest
import time

def timer(func , args = [] , kwargs = {} , repetitions = 10 , comment = ''):
  """
  Function that tests the speed of a given function.
  """
  t0 = time.time() ; c0 = time.clock()
  for i in range(repetitions):
    func(*args,**kwargs)
  cpu_time = time.clock()-c0
  elapsed_time = time.time() - t0
  try: #instance method ?
    name = func.im_class.__name__ + '.' + func.__name__
  except: #ordinary funtion
    name = func.__name__
  print '%s %s (%d calls): elapsed time = %g , CPU = %g' % (comment , name , repetitions , elapsed_time , cpu_time)
  return cpu_time/float(repetitions)


#Demonstration of the functionality of the above function to choose an implementation to calculate the binomialcoefficients. One uses an elementary implementation and the other makes extensive use of python library functions with an underlying C implementation (inherited from functional languages).
#Conclusion binomialCoeff is faster for low n and k comb is faster for big n and k -> use comb in our richardson-gaudin solver program because for low n and k time difference is not important but for big n and k it is.

def binomialCoeff(n, k):
  result = 1
  for i in range(1, k+1):
      result = result * (n-i+1) / i
  return result


from operator import mul
def comb(n,r):
  if r > n-r:  # for smaller intermediate values
      r = n-r
  return int( reduce( mul, range((n-r+1), n+1), 1) /
    reduce( mul, range(1,r+1), 1) )

def main():
  args = [1000 , 500]
  timer(binomialCoeff , args = args, repetitions = 50 )
  timer(comb , args = args, repetitions = 50 )

if __name__ == "__main__":
  main()
