#/usr/bin/python3
import random
import os
import subprocess

def nextInt():
  return random.randint(1, 40)

if __name__ == "__main__":
  cnt = 0
  flag = 1
  while cnt < 1:
    # fill input file
    n = 9 * nextInt()
    print("n = %d"%(n))
    print("-------")
    elem1 = './graph_gen %d'%(n)
    os.system(elem1)
    # Now check output
    elemA = 'mpirun -n %d ./2dsparseFW'%(9)
    os.system(elemA)
    # Now compare
    out1 = []
    with open("output_sparse", mode='r') as o1:
      out1.append(o1.read().strip())
    
    elemB = './classicFW'
    os.system(elemB)
    # Now compare
    out2 = []
    with open("output_classic", mode='r') as o2:
      out2.append(o2.read().strip())

    if(out1 != out2):
      flag = 0
      break
    cnt += 1

  print(flag)