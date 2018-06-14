/*
 * main.cpp
 *
 *  Created on: Nov 15, 2016
 *      Author: leandro
 */

#include "Orbital.h"
#include "gauss.h"
#include "SKtable.h"
#include "scf.h"
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "toml.h"

using namespace std;
using namespace boost;


//double Hubbard(Scf &a,vector<Orbital*> &b,const int &c );

void read_table_opt(char* archivo,SKtable& sk);

struct atom{
	vector<vector<double>> param_dens;
	vector<vector<double>> param_basis;
	vector<string> symbol;
	vector<vector<double>> e;
	vector<vector<double>> U;
	};
atom read_table_atom(char* archivo);
int main(int argc,char *argv[]){


vector<string> simbolo;
vector<vector<double>> e;
vector<double*> ocupation;
vector<vector<double>> U;

vector<Potential_spline*> Veff0;
vector<Potential_spline*> Veff;
vector<Potential_spline*> vconf;

vector<Orbital_spline **> A;

vector<Scf*> Atoms;








//***Lee  Atomo el archivo, y los parametros de confinamiento
atom natom=read_table_atom(argv[1]);
simbolo=natom.symbol;

//********************************************
for(int i=0;i<simbolo.size();i++){
	A.push_back(new Orbital_spline* [3] {});
	e.push_back({0,0,0});
	ocupation.push_back(new double[3]);
	U.push_back({0.,0.,0.});
	Veff0.push_back(new Potential_spline);
	Veff.push_back(new Potential_spline);
	vconf.push_back(new Potential_spline);
	Atoms.push_back(new Scf);
// Inicializa y corre atomo sin confinamiento.
	Atoms[i]->initialize(simbolo[i]);
    Atoms[i]->run(0,1,1,0.2);
    Atoms[i]->energy(e[i],ocupation[i]);
// Corre  confinamiento densidad y obtienen pot.
    Atoms[i]->run(natom.param_dens[i][0],natom.param_dens[i][1],natom.param_dens[i][2],0.2);
    Atoms[i]->Veff_noconf(*(Veff0[i]));
// Corre confinamientos bases y obtiene bases y potenciales.
    Atoms[i]->run(natom.param_basis[i][0],natom.param_basis[i][1],natom.param_basis[i][2],0.2);
	Atoms[i]->orbital(A[i]);
	Atoms[i]->Vconf(*(vconf[i]));
	Atoms[i]->Veff_noconf(*(Veff[i]));
//Remplaza valores e[] y U[] si es necesario
	for (int j=0;j<3;j++){
		if (natom.U[i][j] >= 0){
			U[i][j]=natom.U[i][j];
		}
		if (natom.e[i][j] <= 2.){
			e[i][j]=natom.e[i][j];
		}
					
	}
    
}

cout<<"Calculating skf files"<<endl;
for (int i=0;i<Atoms.size();i++){
	for (int j=0;j<Atoms.size();j++){
		SKtable sk;
		read_table_opt(argv[1],sk);

		if (i !=j){
			sk.create(A[i],A[j],false);
		}
		else {
			sk.create(A[i],A[j],true);
		};

		sk.run(A[i],A[j],*(Veff0[i]),*(Veff0[j]),*(Veff[j]),*(vconf[j]),e[i],U[i],ocupation[i],simbolo[i]+"-"+simbolo[j]+".skf");

		cout<<simbolo[i]+"-"+simbolo[j]+".skf"<<endl;
	}
}


for(int k=0;k<Atoms.size();k++){
	for(int i=0;i<3;i++){
		delete A[k][i];
	}
}

for (int i=0;i<Atoms.size();i++){
	delete Atoms[i];
    delete Veff[i];
    delete vconf[i];
    delete Veff0[i];
    delete [] ocupation[i];


}


//*******************
return 0;
}




void read_table_opt(char* archivo,SKtable& sk){
	ifstream entrada(archivo);
		
	toml::ParseResult pr = toml::parse(entrada);
	toml::Value v  = pr.value;
    //Parseo de opciones
	if (!pr.valid()) {
	    cout << pr.errorReason << endl;
	}
	toml::Value* x = v.find("TableOptions.ngrid");
	if (x and x->is<int>()){
		sk.set_grid(x->as<int>());
	}
	
	x=v.find("TableOptions.rmax");	
	if (x and x->is<double>()){
		sk.set_rmax(x->as<double>());
	}
	
    x=v.find("TableOptions.rmin");	
	if (x and x->is<double>()){
		sk.set_rmin(x->as<double>());
	}
    
    x=v.find("TableOptions.step");	
	if (x and x->is<double>()){
		sk.set_step(x->as<double>());
	}
	entrada.close();
}	
atom read_table_atom(char* archivo)  {
    //Parseo de atomos
    ifstream entrada(archivo);
    toml::ParseResult pr = toml::parse(entrada);
	toml::Value v  = pr.value;
    //Parseo de opciones
	if (!pr.valid()) {
	    cout << pr.errorReason << endl;
	}
    atom atparam;
    const toml::Value* z = v.find("TableOptions.Atoms");
    const toml::Array& ar = z->as<toml::Array>();
	for (const toml::Value& v : ar){
		atparam.symbol.push_back(v.get<string>("element"));
				
		vector<double> bas={0,0,0};
		vector<double> dens={0,0,0};
		vector<double> Uparse={0.,0.,0.};
		vector<double> eparse={2.1,2.1,2.1};
		
		z=v.find("basis");
		if(z and z->is<vector<double>>()){
			bas=z->as<vector<double>>();
		}
		atparam.param_basis.push_back(bas);
		
		z=v.find("dens");
		if(z and z->is<vector<double>>()){
			 dens=z->as<vector<double>>();
		}
		else{dens=bas;} 
		atparam.param_dens.push_back(dens);
		
		z=v.find("Ud");
		if(z and z->is<double>()){
			 Uparse[2]=z->as<double>();
		}
		
		z=v.find("Up");
		if(z and z->is<double>()){
			 Uparse[1]=z->as<double>();
		}
		
		z=v.find("Us");
		if(z and z->is<double>()){
			 Uparse[0]=z->as<double>();
		}
		
		z=v.find("ed");
		if(z and z->is<double>()){
			 eparse[2]=z->as<double>();
		}
		
		z=v.find("ep");
		if(z and z->is<double>()){
			 eparse[1]=z->as<double>();
		}
		
		z=v.find("es");
		if(z and z->is<double>()){
			 eparse[0]=z->as<double>();
		}
		
		
		atparam.U.push_back(Uparse);
		atparam.e.push_back(eparse);
		
	}
			  

	entrada.close();
return atparam;
};
/*
double Hubbard(Scf &a,vector<Orbital*> &b,const int &c ){
	double h=0.001;
	double n=b[c]->ocup();
	double nf=n+0.5*h;
	double nb=n-0.5*h;
	double e=b[c]->energy();
	double ef;
	double eb;
	b[c]->set_ocup(nf);
	a.run(b,0,1,1,0.3);
	ef=b[c]->energy();
	b[c]->set_ocup(nb);
	a.run(b,0,1,1,0.3);
	eb=b[c]->energy();
	b[c]->set_ocup(n);
	a.run(b,0,1,1,0.3);
	return (ef-eb)/h;
}
*/

