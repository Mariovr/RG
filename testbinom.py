
from datetime import datetime

def binomialCoeff(n, k):
  result = 1
  for i in range(1, k+1):
      result = result * (n-i+1) / i
  return result


from operator import mul
def comb(n,r):
  ''' calculate nCr - the binomial coefficient
  >>> comb(3,2)
  3
  >>> comb(9,4)
  126
  >>> comb(9,6)
  84
  >>> comb(20,14)
  38760
  '''

  if r > n-r:  # for smaller intermediate values
      r = n-r
  return int( reduce( mul, range((n-r+1), n+1), 1) /
    reduce( mul, range(1,r+1), 1) )

 
if __name__ == "__main__":
  number = 4
  a = datetime.now()
  for i in range(number):
    binomialCoeff(number, i)
  print(binomialCoeff(number, 2))
  b = datetime.now()
  print b-a
  a = datetime.now()
  for i in range(number):
    comb(number,i)
  print(comb(number, 2))
  b = datetime.now()
  print b-a
  
  