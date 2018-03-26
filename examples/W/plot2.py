import numpy as np
import matplotlib.pyplot as plt
import sys

pi=3.14159265359
a=np.loadtxt("w_band_tot.dat")
b=np.loadtxt("wbandqe.dat")
Efermi=sys.argv[1]
efermi=float(Efermi)*np.ones(np.shape(a))
b=np.sort(b)

for i in range(1,len(a[0,:])):
    plt.plot(a[:,0],(a[:,i]-efermi[:,1]),color="r")

for i in range(4,len(b[0,:])):
    plt.plot(a[:,0],b[:,i],color="b")

plt.show()
