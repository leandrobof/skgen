XCDIR=/opt/etsf
INCLDIR=../include

all: main2.o Orbital.o Func.o xc_potential.o hartree.o Poisson.o Integrand.o gauss.o scf.o SKtable.o
	g++ main2.o scf.o Orbital.o Func.o xc_potential.o hartree.o Poisson.o Integrand.o gauss.o SKtable.o  -o   skgen -L$(XCDIR)/lib -lgsl -lopenblas -lm -lxc -std=c++11

main2.o: main2.cpp $(INCLDIR)/Orbital.h
	g++  main2.cpp -c -I$(INCLDIR) -I$(XCDIR)/include -std=c++11

SKtable.o: SKtable.cpp $(INCLDIR)/SKtable.h
	g++  SKtable.cpp -c -I$(INCLDIR) -std=c++11 

gauss.o:gauss.cpp $(INCLDIR)/gauss.h
	g++ gauss.cpp -c -I$(INCLDIR) -std=c++11

Integrand.o:Integrand.cpp $(INCLDIR)/Integrand.h
	g++ Integrand.cpp -c -I$(INCLDIR) -std=c++11

scf.o:scf.cpp $(INCLDIR)/scf.h
	g++ scf.cpp -c -I$(XCDIR)/include -I$(INCLDIR) -std=c++11 

Orbital.o: Orbital.cpp $(INCLDIR)/Orbital.h 
	g++  Orbital.cpp -c -O3 -I$(INCLDIR)  -std=c++11

Func.o:Func.cpp $(INCLDIR)/Func.h 
	g++ Func.cpp -c -O3 -I$(INCLDIR) -std=c++11

xc_potential.o: xc_potential.cpp $(INCLDIR)/xc_potential.h
	g++ xc_potential.cpp   -c -O3 -I$(XCDIR)/include -I$(INCLDIR)

hartree.o: hartree.cpp $(INCLDIR)/hartree.h Poisson.cpp $(INCLDIR)/Poisson.h
	g++ hartree.cpp -c -O3 -I$(INCLDIR) -std=c++11

Poisson.o: Poisson.cpp $(INCLDIR)/Poisson.h
	g++ Poisson.cpp -c -O3 -I$(INCLDIR)

clean:
	rm -rf *.o
	rm atomdft
