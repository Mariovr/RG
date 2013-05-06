from itertools import combinations,combinations_with_replacement
import numpy as np
tdacombinations = combinations_with_replacement(np.arange(0,13),10)
senioriteit = np.zeros(13)
ontaardingen = [4,4,4,8,4,4,8,8,4,8,4,8,12]
file = open('correspondentie.dat', 'w')
wind = 0
startw = -1
endw = 310000

for tdadict in tdacombinations:
  tdastartd = {}
  goodsol = True
  for i in tdadict:
    a = tdadict.count(i)      
    tdastartd[i] = a
    if a*2 + senioriteit[i]*2 > ontaardingen[i]:
      goodsol = False
  if goodsol == True:
    wind += 1
    if wind >= startw and wind < endw:
      print 'the %g th state corresponds to %s ' %( wind,str(tdastartd))
      file.write('%g\t%s\n' %(wind,str(tdastartd)))
