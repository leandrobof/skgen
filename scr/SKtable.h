/*
 * SKtable.h
 *
 *  Created on: Feb 27, 2018
 *      Author: leandro
 */
#include "globals.h"
#include "Orbital.h"
#include "gauss.h"


#ifndef SKTABLE_H_
#define SKTABLE_H_

class SKtable {
private:
	double step;
	double rmin;
	double rmax;
	int ngrid;
	vector<Integrand*> S;
	vector<Integrand*> V;
	bool homonuclear;
public:
	SKtable();
	virtual ~SKtable();
	void create( Orbital_spline **a,Orbital_spline **b,bool);
    void set_grid(int );
    void set_rmax(double);
    void set_rmin(double);
    void set_step(double);
    void run( Orbital_spline **A, Orbital_spline **B,Potential_spline &Veff,Potential_spline &Veff0b,Potential_spline &Veffb,Potential_spline &Vconf,vector<double> &e,vector<double> &U,double *ocupation,string archivo);
};

#endif /* SKTABLE_H_ */
