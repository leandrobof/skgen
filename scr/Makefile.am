AM_CXXFLAGS = -std=c++11
ACLOCAL_AMFLAGS = -I m4


bin_PROGRAMS= skgen
skgen_SOURCES=   main2.cpp Orbital.cpp Orbital.h Func.cpp Func.h xc_potential.cpp xc_potential.h hartree.cpp hartree.h Poisson.cpp Poisson.h Integrand.cpp Integrand.h  gauss.cpp gauss.h SKtable.cpp SKtable.h scf.cpp scf.h
skgen_CXXFLAGS = $(XC_CFLAGS) $(BOOST_CPPFLAGS) -std=c++11
skgen_LDADD =  $(GSL_LIBS) $(XC_LIBS)

