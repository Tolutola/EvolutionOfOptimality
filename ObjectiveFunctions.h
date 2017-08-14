#include<vector>
#include<iostream>
#include<tuple>
#include<cmath>
#include <stdio.h>
#include<string.h>
#include<algorithm>

#include "mosek.h"

using namespace std;
#ifndef OBJ_H_
#define OBJ_H_

template<class T>
// maximization of biomass yield
vector<T> maxBiomass(const MSKint32t numvar, const MSKint32t numcon, vector<double> c, vector<MSKint32t> aptrb1,vector<MSKint32t> aptre1,
		vector<MSKint32t> asub1, vector<double> aval1, vector<MSKboundkeye> bkc1, vector<double> blc1, vector<double> buc1,
		vector<MSKboundkeye> bkx1, vector<double> blx1, vector<T> bux1,int Osense);

#include"maxBiomass.tpp"

template<class T>
// minimization of the euclidean norm of fluxes (used in maximization of biomass (and atp) yield  per unit flux
 vector<T> minEuclideanNorm(const MSKint32t NUMVAR, const MSKint32t NUMCON, vector<double> c1, vector<MSKint32t> aptrb1,vector<MSKint32t> aptre1,
		vector<MSKint32t> asub1, vector<double> aval1, vector<MSKboundkeye> bkc1, vector<double> blc1, vector<double> buc1,
		vector<MSKboundkeye> bkx1, vector<double> blx1, vector<T> bux2,
		int NUMANZ, int NUMQNZ,vector<MSKint32t> qsubi1, vector<MSKint32t> qsubj1,vector<double> qval1);

#include "minEuclideanNorm.tpp"


template<class T>
// maximization of biomass yield per unit flux
vector<T> maxBiomass_PF(const MSKint32t numvar, const MSKint32t numcon, vector<double> c, vector<MSKint32t> aptrb1,vector<MSKint32t> aptre1,
		vector<MSKint32t> asub1, vector<double> aval1, vector<MSKboundkeye> bkc1, vector<double> blc1, vector<double> buc1,
		vector<MSKboundkeye> bkx1, vector<double> blx1, vector<T> bux1);

#include "maxBiomass_PF.tpp"

template<class T>
// maximization of ATP yield
vector<T> maxATP(const MSKint32t numvar, const MSKint32t numcon, vector<double> c, vector<MSKint32t> aptrb1,vector<MSKint32t> aptre1,
		vector<MSKint32t> asub1, vector<double> aval1, vector<MSKboundkeye> bkc1, vector<double> blc1, vector<double> buc1,
		vector<MSKboundkeye> bkx1, vector<double> blx1, vector<T> bux1,T atpFlux);

#include "maxATP.tpp"

template<class T>
// maximization of ATP yield per unit flux
vector<T> maxATP_PF(const MSKint32t numvar, const MSKint32t numcon, vector<double> c, vector<MSKint32t> aptrb1,vector<MSKint32t> aptre1,
		vector<MSKint32t> asub1, vector<double> aval1, vector<MSKboundkeye> bkc1, vector<double> blc1, vector<double> buc1,
		vector<MSKboundkeye> bkx1, vector<double> blx1, vector<T> bux1, T atpFlux);

#include "maxATP_PF.tpp"


//template<class T>
//vector<T> RemoveLoops(const MSKint32t numvar, const MSKint32t numcon, vector<double> c, vector<MSKint32t> aptrb1,vector<MSKint32t> aptre1,
//		vector<MSKint32t> asub1, vector<double> aval1, vector<MSKboundkeye> bkc1, vector<double> blc1, vector<double> buc1,
//		vector<MSKboundkeye> bkx1, vector<double> blx1, vector<T> bux1);
//
//#include"RemoveLoops.tpp"

#endif
