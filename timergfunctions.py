import sys
import os , unittest
import time

def timer(func , args = [] , kwargs = {} , repetitions = 10 , comment = ''):
  t0 = time.time() ; c0 = time.clock
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



if __name__ == "__main__":
  main():
