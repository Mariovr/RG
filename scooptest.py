#
#SCOOP is free software: you can redistribute it and/or modify
#it under the terms of the GNU Lesser General Public License as
#published by the Free Software Foundation, either version 3 of
#the License, or (at your option) any later version.
#
#SCOOP is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU Lesser General Public License for more details.
#
#You should have received a copy of the GNU Lesser General Public
#License along with SCOOP. If not, see <http://www.gnu.org/licenses/>.
#
"""
A simple example showing how to resolve a full balanced tree with multiples
techniques using SCOOP.
"""
import datetime
import scoop
from scoop import futures

def func4(n):
  # Example of calculus you may want to perform
  result = n*n
  #print scoop.WORKER_NAME ,' returned this result ', result
  #print 'at this moment we have %g workers' % scoop.SIZE
  return result
  
def func3(n):
  # This call results in a generator function
  result = futures.map(func4, [i+1 for i in range(n)])
  # The results are evaluated here when they are accessed.
  return sum(result)

def func2(n):
  launches = [futures.submit(func3, i + 1) for i in range(n)]
  # Spawn a generator for each completion, unordered
  result = (a.result() for a in futures.as_completed(launches))
  return sum(result)

def func1(n):
  # To force an immediate evaluation, you can wrap your map in a list such as:
  result = list(futures.map(func2, [i+1 for i in range(n)]))
  return sum(result)
  
def func0(n):
  # Task submission is asynchronous; It will return immediately.
  task = futures.submit(func1, n)
  # The call blocks here until it gets the result
  result = task.result()
  return result

def main(n):
  #creates a tree of all different workers (a worker goes down and generates other workers who generate other workers)
  n =n
  task = futures.submit(func0, n)
  # You can wait for a result before continuing computing
  futures.wait([task], return_when=futures.ALL_COMPLETED)
  result = task.result()
  print(result)
  return result


def main2(n):
  # This call results in a generator function
  result = futures.map(func4, [i+1 for i in range(n)])
  print result
  # The results are evaluated here when they are accessed.
  d = sum(result)
  print d
  return d


def mainserial(n):
  n =n
  result = sum(map(func4, [i+1 for i in range(n)] ))
  print result
  return result


if __name__ == "__main__":
  n = 100000
  start = datetime.datetime.now() 
  presult = main2(n)
  end1 = datetime.datetime.now()
  sresult = mainserial(n)
  end2 = datetime.datetime.now()
  print 'the total execution time of the parallel program was: %f with result = %f' %( (end1-start).total_seconds() , presult)
  print 'the total execution time of the serial program was: %f with result = %f' %( (end2-end1).total_seconds(), sresult )


