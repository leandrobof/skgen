import numpy as np
import matplotlib.pyplot as plt
import sys


a=np.loadtxt("BN_band.dat")
b=np.loadtxt("BN_band_tot.dat")
r=range(0,len(b[:,0]))
Efermi=sys.argv[1]
efermi=float(b[240,-9])*np.ones(np.shape(b))


for i in range(2,len(a[0,:])):
    plt.plot(r,a[:,i],color="r")

for j in range(3,len(b[0,:-8])):
    plt.plot(r,(b[:,j]-efermi[:,1]),color="b")

plt.show()
