import sys
import numpy as np
import matplotlib.pylab as plt

SNAP = sys.argv[1]
level = int(sys.argv[2])
N = (2**(level) + 2)**3

N =  N
f = open(SNAP)

data = np.fromfile(f,dtype=float, count=(6*N ))

C1 = []
C2 = []
C3 = []
C4 = []
C5 = []

i=0
while (i < 5*N):
    C1.append(data[i])
    C2.append(data[i+1])
    C3.append(data[i+2])
    C4.append(data[i+3])
    C5.append(data[i+4])
    i=i+5

print len(C1)
plt.plot(C1)  
plt.show()
C1 = np.array(C1).reshape((2**level + 2, 2**level + 2, 2**level + 2)).T

#rho = np.resize(rho,(2**level, 2**level, 2**level))

#plt.hist(C1[1], histtype='step')
#plt.hist(C1[1], histtype='step')
#plt.hist(C1[2], histtype='step')

print C1

#plt.plot(C1)
#print C1[2]
plt.figure(figsize=(18,6))
plt.subplot(131)
plt.imshow(C1[1])

plt.subplot(132)
plt.imshow(C1[len(C1)/2])

plt.subplot(133)
plt.imshow(C1[len(C1) - 1])
plt.show()
#print len(data)
#print N
#print len(data)/(N*6)
print data
