import numpy as np
import matplotlib.pyplot as plt
import sys
import os


os.system("cp bulk dftb_in.hsd"  )
os.system("dftb+"  )
os.system("dp_bands band.out band"  )
b=np.loadtxt("wbandqe.dat")
for line in open("detailed.out","r"):
    if "Fermi energy" in line:
        Efermi=line.split()[4]


a=np.loadtxt("band_tot.dat")


efermi=float(Efermi)*np.ones(np.shape(a))
b=np.sort(b)

for i in range(1,len(a[0,:])):
    plt.plot(a[:,0],(a[:,i]-efermi[:,1]),color="r")

for i in range(0,len(b[0,:])):
    plt.plot(a[:,0],b[:,i],color="b")

plt.show()
