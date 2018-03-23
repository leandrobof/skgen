/*
 * Integrand.h
 *
 *  Created on: Apr 5, 2017
 *      Author: leandro
 */

#ifndef INTEGRAND_H_
#define INTEGRAND_H_
#include <gsl/gsl_spline.h>
#include <math.h>
#include <iostream>
#include "Orbital.h"

using namespace std;


class Integrand {
protected:


public:
	Integrand();
	virtual ~Integrand();
	virtual double calculate(double *R1,double *R2,double ,double,double,double &sin1,double &sin2,double &cos1,double &cos2) = 0;
    virtual double energy(){return 0;};
};
/*
class Gaussianas : public Integrand {
private:
    double z1;
    double z2;
    double N1;
    double N2;

public:
    Gaussianas(double ,double ,double );
    double gaussian1(double );
    double gaussian2(double );
    double real(double );
    virtual double operator()(double ,double );
};
#endif /* INTEGRAND_H_

class Density: public Integrand{
private:
	gsl_spline *rho1;
	gsl_spline *rho2;
	gsl_interp_accel *accrho1;
	gsl_interp_accel *accrho2;
    double r1max;
    double r2max;
public:
	Density(gsl_spline *,gsl_interp_accel *,gsl_spline *,gsl_interp_accel *,double,double,double);
	virtual double calculate(double &r1,double &r2,double &sin1,double &sin2,double &cos1,double &cos2);

};
*/
class S : public Integrand{
protected:
	int l1;
	int l2;

public:
    S();
    inline virtual double f(double sin1,double sin2,double cos1,double cos2){return 1;}
    virtual double calculate(double *R1,double *R2,double v1,double Db,double vconf,double &sin1,double &sin2,double &cos1,double &cos2);
};
class S_ss : public S{
public:
	S_ss( );
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return 0.5;}

};
class S_sp :public S{
public:
	S_sp( );
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return sqrt(3)/2.*cos2;}
};
class S_pp_sig: public S{
public:
	S_pp_sig( );
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return (3/2.)*cos1*cos2;}
};
class S_pp_pi : public S{
public:
	S_pp_pi( );
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return (3/4.)*sin1*sin2;}
};
class S_sd :public S{
public:
	S_sd( );
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return sqrt(5)/4.*(3*cos2*cos2-1);}
};

class S_pd_pi : public S{
public:
	S_pd_pi( );
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return (1*sqrt(45)/4.)*sin1*sin2*cos2;}
};
class S_pd_sig : public S{
public:
	S_pd_sig( );
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return 1*sqrt(15)/4.*cos1*(3*cos2*cos2-1);}
};
class S_dd_sig : public S{
public:
	S_dd_sig( );
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return 5/8.*(3*cos2*cos2-1)*(3*cos1*cos1-1);}
};
class S_dd_pi : public S{
public:
	S_dd_pi( );
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return (15/4.)*sin1*cos1*sin2*cos2;}
};
class S_dd_del : public S{
public:
	S_dd_del( );
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return (15/16.)*sin1*sin1*sin2*sin2;}
};
class V: public S{
private:
      double e;
public:
	V(double en)  ;
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return 1;}
	virtual double calculate(double *R1,double *R2,double v1,double Db,double vconf,double &sin1,double &sin2,double &cos1,double &cos2);
    virtual double energy(){return e;}
};

class V_ss : public V{
public:
	V_ss(double en);
    inline virtual double f(double sin1,double sin2,double cos1,double cos2){return 0.5;}

};
class V_sp :public V{
public:
	V_sp(double en);
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return sqrt(3)/2.*cos2;}
};
class V_pp_sig: public V{
public:
	V_pp_sig(double en);
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return (3/2.)*cos1*cos2;}
};
class V_pp_pi : public V{
public:
	V_pp_pi(double en);
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return (3/4.)*sin1*sin2;}
};
class V_sd :public V{
public:
	V_sd(double en);
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return sqrt(5)/4.*(3*cos2*cos2-1);}
};
class V_pd_pi : public V{
public:
	V_pd_pi(double en);
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return (1*sqrt(45)/4.)*sin1*sin2*cos2;}
};
class V_pd_sig :public V{
public:
	V_pd_sig(double en);
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return 1*sqrt(15)/4.*cos1*(3*cos2*cos2-1);}
};
class V_dd_del : public V{
public:
	V_dd_del(double en);
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return (15/16.)*sin1*sin1*sin2*sin2;}
};
class V_dd_pi : public V{
public:
	V_dd_pi(double en);
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return (15/4.)*sin1*cos1*sin2*cos2;}
};
class V_dd_sig :public V{
public:
	V_dd_sig(double en);
	inline virtual double f(double sin1,double sin2,double cos1,double cos2){return 5/8.*(3*cos2*cos2-1)*(3*cos1*cos1-1);}
};
#endif /* INTEGRAND_H_*/
