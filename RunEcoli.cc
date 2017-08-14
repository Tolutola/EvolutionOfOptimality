#include <iostream>
#include<algorithm>

#include "Matrix.h"
#include "ObjectiveFunctions.h"
#include "auxillaryFunctions.h"

int main () {
	// load genome-scale metabolic reconstruction of E. coli
	// read in data from text files

	string filename ="exchangeRxns.txt";
	vector<int> ExcIndex=readfile(filename,1,1);
	
	// import lower and upper bounds (on flux values)
	filename="lb.txt";
	vector<double> Lb=readfile(filename,double(1),2);


	filename="ub.txt";
	vector<double> Ub=readfile(filename,double(1),1);
	
	// import bounds on constraints 

	filename="blc.txt";
	vector<double> blc=readfile(filename,double(1),2); 
	vector<double>buc=blc; // for linear equality constraints, upper and lower bounds are the same.
	
	// import objective function
	filename="objectiveVec.txt";
	vector<double> c=readfile(filename,double(1),1);
	

	// import Stoichiometric matrix as a sparse matrix in mosek format (Aval, aptrb, aptre, asub)
	filename="AVAL.txt"; //non zero values of S
	vector<double> aval=readfile(filename,double(1),1);
	

	filename="ASUB.txt"; // positions of the non zero values
	vector<MSKint32t> asub=readfile(filename,MSKint32t(2),1);
	

	filename="APTRB.txt"; // pointers to the first non zero index in a column
	vector<int> aptrb=readfile(filename,int(1),1);
	
	filename="APTRE.txt"; // pointers to the last non zero index in a column
	vector<MSKint32t> aptre=readfile(filename,MSKint32t(1),1);
	
	// the program simulates a batch reactor of ecoli cells growing on glucose (and acetate)
	// I use a forward euler type solver (instead of my own ODE solver implementation)
	// that way I can monitor what goes on easily

	//	// basic parameters of the simulation
	int noOfSteps=2000;
	int noOfRxns=2382;
	int noOfExcRxns=300;
	int noOfMets=1668;
	int bioIndex=1004;
	int atpIndex=374;
	int noObj=4;
	double timeStep=0.25; // hours
	double initBioConc=0.001; //starting concentration of ecoli cells of different types (gDCW/L)

	// change bounds on key reactions to reflect physiological conditions
	Lb[atpIndex]=8.39;
	Ub[atpIndex]=1000;

	Lb[848]=-9.4; // glucose mmol/gDCW/L
	Ub[848]=0;

	Lb[932]=-18; // oxygen mmol/gDCW/L

	
	//	// bound keys for the different bounds on constraints and variables (needed by mosek)
	vector<MSKboundkeye> bkx(noOfRxns);
	std::fill (bkx.begin(),bkx.end(),MSK_BK_RA); 
	
	vector<MSKboundkeye> bkc(noOfMets);
	std::fill (bkc.begin(),bkc.end(),MSK_BK_FX); 
	

	vector<int> substrateIndex={147}; // glucose uptake index in the external reactions vector
	vector<double> initSubstrateConc={20};//mM of glucose (other substrates could be co-fed hence a vector is used)


	const MSKint32t numvar=noOfRxns, numcon=noOfMets;

	//	//initialize biomassMatrix
	double tempBiomass=0;
	Matrix<double> biomassMatrix(1,noObj); 
	for (int i=0;i!=noObj;++i){
		biomassMatrix(0,i)=initBioConc;
		tempBiomass+=biomassMatrix(0,i);
	}
	//
	//	// initialize concentrations Matrix
	//	//OriginalBounds for exchange reactions
	vector<double> OriginalBounds(noOfExcRxns);
	for (int n=0;n!=noOfExcRxns;++n) OriginalBounds[n]=-Lb[ExcIndex[n]];
	//
	Matrix<double> concentrationMatrix(1,noOfExcRxns);
	//
	//	// take care of unknown concentrations (by first filling up all concentrations with a high number)
	for (int n=0;n!=noOfExcRxns;++n) {
		if (OriginalBounds[n]>0) concentrationMatrix(0,n)=1000; // mM
		else concentrationMatrix(0,n)=0;
	}
	//
	//	//fill up known starting reactor concentrations (in this case just glucose)
	for (int n=0;n!=substrateIndex.size();++n) concentrationMatrix(0,substrateIndex[n])=initSubstrateConc[n];
	//
	//	// initalize uptake bounds
	vector<double> uptakeBounds(noOfExcRxns);

	//
	for (int n=0;n!=noOfExcRxns;++n)  {
		uptakeBounds[n]=concentrationMatrix(0,n)/tempBiomass/timeStep;
		// make sure bounds are not higher than that specified by the metabolic reconstruction
		if (uptakeBounds[n]>OriginalBounds[n] & OriginalBounds[n]>0){
			uptakeBounds[n]=OriginalBounds[n];
		}
		// change the uptake bounds in the metabolic network
		Lb[ExcIndex[n]]=-uptakeBounds[n];

	}
	//
	//
	//	// initialize uptake fluxMatrix
	Matrix<double> uptakeFlux(noOfExcRxns,noObj);
	//
	//	// initialize solution vector
	vector<double> solution1(noOfRxns),solution2(noOfRxns),solution3(noOfRxns),solution4(noOfRxns);
	//
	//	// initialize temporary variables to be used in the loop 
	vector<double> tempGrate(noObj); //growth rate vector
	double tempConsumed; // substrate consumed
	//
	//	// initialize time vector
	vector<double> timeVector;
	timeVector.push_back(0);
	//
	//	// initialize flux matrices for each objective funtion
	//	// these would store the entire time profile of the flux ditributions for each objective function
	Matrix<double> FluxMatrix1(1,noOfRxns);
	Matrix<double> FluxMatrix2(1,noOfRxns);
	Matrix<double> FluxMatrix3(1,noOfRxns);
	Matrix<double> FluxMatrix4(1,noOfRxns);
	//
	//	// print
	cout<<"\n\nStepNo\t"<<"MaxBio\t"<<"MaxATP\t"<<"MaxBPF\t"<<"MaxAPF\t"<<endl;
	cout<<0<<"\t";
	for (int m=0;m!=noObj;++m) cout<<biomassMatrix(0,m)<<"\t";
	cout<<endl;

	// main loop 
	for (int i=0; i!=noOfSteps;++i){

		// solve the optimization problems for the different objective functions

		// maximize biomass production

		solution1= maxBiomass(numvar,numcon,c,aptrb,aptre,asub,aval, bkc, blc, buc,bkx, Lb,Ub,1);

		if (solution1[bioIndex]<=0.0000001) {
			cout<<"No biomass growth, nutrients exhausted"<<endl;
			break;
		}
		solution2=maxATP(numvar,numcon,c,aptrb,aptre,asub,aval, bkc, blc, buc,bkx, Lb,Ub, solution1[atpIndex]);
		solution3=maxBiomass_PF(numvar,numcon,c,aptrb,aptre,asub,aval, bkc, blc, buc,bkx, Lb,Ub,5);
		solution4=maxATP_PF(numvar,numcon,c,aptrb,aptre,asub,aval, bkc, blc, buc,bkx, Lb,Ub,solution1[atpIndex]);


		tempGrate={solution1[bioIndex],solution2[bioIndex],solution3[bioIndex],solution4[bioIndex]};

		// update the flux distribution matrices
		FluxMatrix1.addRow(0);
		FluxMatrix2.addRow(0);
		FluxMatrix3.addRow(0);
		FluxMatrix4.addRow(0);

		for (int m=0;m!=noOfRxns;++m) {
			FluxMatrix1(i+1,m)=solution1[m];
			FluxMatrix2(i+1,m)=solution2[m];
			FluxMatrix3(i+1,m)=solution3[m];
			FluxMatrix4(i+1,m)=solution4[m];
		}

		// update biomass Matrix
		biomassMatrix.addRow(0.0);
		tempBiomass=0;
		for (int m=0;m!=noObj;++m) {
			biomassMatrix(i+1,m)=biomassMatrix(i,m)*exp(tempGrate[m]*timeStep);
			tempBiomass+=biomassMatrix(i+1,m);
		}

		concentrationMatrix.addRow(0.0);
		for (int n=0;n!=noOfExcRxns;++n) {
			// fill up uptakeFlux Matrix
			uptakeFlux(n,0)=solution1[ExcIndex[n]];
			uptakeFlux(n,1)=solution2[ExcIndex[n]];
			uptakeFlux(n,2)=solution3[ExcIndex[n]];
			uptakeFlux(n,3)=solution4[ExcIndex[n]];

			// update concentrations
			tempConsumed=0;
			for (int m=0;m!=noObj;++m){
				tempConsumed+=(uptakeFlux(n,m)/tempGrate[m]*biomassMatrix(i+1,m)*(1-exp(tempGrate[m]*timeStep)));
			}
			concentrationMatrix(i+1,n)=concentrationMatrix(i,n)-tempConsumed;
			if (concentrationMatrix(i+1,n)<0) concentrationMatrix(i+1,n)=0; // no negative concetrations allowed

			// update bounds for uptake reactions
			uptakeBounds[n]= concentrationMatrix(i+1,n)/tempBiomass/timeStep;
			if (uptakeBounds[n]>1000) uptakeBounds[n]=1000; // to avoid numerical issues
			if (uptakeBounds[n]>OriginalBounds[n] & OriginalBounds[n]>0) uptakeBounds[n]=OriginalBounds[n]; // revert to the original
			if (abs(uptakeBounds[n])<0.0000000001) uptakeBounds[n]=0; // to avoid numerical issues

			// change the uptake bounds in the metabolic network
			Lb[ExcIndex[n]]=-uptakeBounds[n];

		}

		// uptake time vector
		timeVector.push_back(i*timeStep);

		// print out biomass concentrations in the reactor
		cout<<i+1<<"\t";
		for (int m=0;m!=noObj;++m) cout<<biomassMatrix(i+1,m)<<"\t";
		cout<<endl;


	}
	//store flux distribution matrices in text files
	writeFluxDistributions(FluxMatrix1,FluxMatrix2,FluxMatrix3,FluxMatrix4);

	cout<<"program completed\n"<<"flux distribution matrices stored as text files"<<endl;

}