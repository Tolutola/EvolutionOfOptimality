#include "ObjectiveFunctions.h"
#include "mosek.h"


template<class T>
vector<T> maxBiomass_PF(const MSKint32t numvar, const MSKint32t numcon, vector<double> c1, vector<MSKint32t> aptrb1,vector<MSKint32t> aptre1,
		vector<MSKint32t> asub1, vector<double> aval1, vector<MSKboundkeye> bkc1, vector<double> blc1, vector<double> buc1,
		vector<MSKboundkeye> bkx1, vector<double> blx1, vector<T> bux2,int noSteps) {

    vector<T> myResults(numvar);
    // first compute maximum and minimum biomass growth rates
    vector<T> Objsol= maxBiomass(numvar,numcon,c1,aptrb1,aptre1,asub1,aval1, bkc1, blc1, buc1,bkx1, blx1,bux2,1);   
    T bmMax=Objsol[1004]; 
    double tmpObj; // temporary variable for euclidean norm of fluxes    
    for(int k=0; k!=numvar;++k) tmpObj+=pow(Objsol[k],2);
    
    T bmMin =0 ; // set minimum biomass growth rate as zero
    
    int NUMANZ = 9236; // number of non zero elements in A
    int NUMQNZ=2382;  // number of non zero elements in Q
    
    vector<MSKint32t> qsubi1(numvar);
    vector<MSKint32t> qsubj1 (numvar);
    vector<double> qval1 (numvar);
    
    for (int i=0;i!=numvar;++i) {
    	qsubi1[i]=i;
    	qsubj1[i]=i;
    	qval1[i]=double(2);
    }
    
    c1[1004]=double(0); // no linear part of the quadratic programming problem (thus setting biomass coefficient to zero
    
    vector<double> objValues; // vector to store valid results of the optimization problem
    double tmp; // temporary variable to hold value of biomass growth
   	double Objmax=bmMax/tmpObj/2; //max biomass per unit flux
   	
   	// temporary vectors to hold results of optimization
    vector<double> tmpVec;
     
   
    for (int i=0;i!=noSteps;++i) {
    	 tmp=bmMin+i*(bmMax-bmMin)/noSteps;
    	 blx1[1004]=tmp;
    	 bux2[1004]=tmp;
    	 bkx1[1004]=MSK_BK_FX;
    	 
    	 vector<double> tmpVec;
    	 tmpVec= minEuclideanNorm(numvar,numcon,c1,aptrb1,aptre1,asub1,aval1, bkc1,  blc1, buc1,bkx1, blx1,bux2,NUMANZ, NUMQNZ,qsubi1,qsubj1,qval1);
    	 
    	 if (tmpVec[1004]>0) {
    	 	for(int k=0; k!=numvar;++k) tmpObj+=pow(tmpVec[k],2);
    	 	
    	 	if ((tmp/2/tmpObj)>Objmax){
    	 		Objmax=tmp/2/tmpObj;
    	 		Objsol=tmpVec;
    	 	}
  		
    	 }
    	 
   }
    
    for(int k=0;k!=numvar;++k) 	myResults[k]=T(Objsol[k]);
   	
	
	return myResults;
}