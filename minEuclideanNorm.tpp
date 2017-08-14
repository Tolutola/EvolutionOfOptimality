#include "ObjectiveFunctions.h"
#include "mosek.h"
#include <stdio.h>

	

template<class T>
vector<T> minEuclideanNorm(const MSKint32t NUMVAR, const MSKint32t NUMCON, vector<double> c1, vector<MSKint32t> aptrb1,vector<MSKint32t> aptre1,
		vector<MSKint32t> asub1, vector<double> aval1, vector<MSKboundkeye> bkc1, vector<double> blc1, vector<double> buc1,
		vector<MSKboundkeye> bkx1, vector<double> blx1, vector<T> bux2,
		int NUMANZ, int NUMQNZ,vector<MSKint32t> qsubi1, vector<MSKint32t> qsubj1,vector<double> qval1) {
	


		
  vector<T> myResults(NUMVAR);
  
 
  
  // initialize the variables
	vector<double> bux1(NUMVAR);
	for (int i=0;i!=NUMVAR;++i) bux1[i]=double(bux2[i]);
   
  MSKint32t     qsubi[NUMQNZ];
  MSKint32t     qsubj[NUMQNZ];
  double        qval[NUMQNZ];
  
  MSKint32t     i,j;
  double        xx[NUMVAR];

  MSKenv_t      env = NULL;
  MSKtask_t     task = NULL;
  MSKrescodee   r;
  
  
  
  // copy into arrays	
  MSKint32t aptrb[aptrb1.size()], aptre[aptre1.size()], asub[asub1.size()];
	double aval[aval1.size()],blc[NUMCON], buc[NUMCON],blx[NUMVAR],bux[NUMVAR],c[NUMVAR];
	MSKboundkeye bkc[NUMCON],bkx[NUMVAR];
  
	 for(int k=0;k!=NUMCON;++k) {
	 	bkc[k]=bkc1[k];
	 	blc[k]=blc1[k];
	 	buc[k]=buc1[k];
	 	
	 }
	
	 for(int k=0;k!=aptrb1.size();++k) aptrb[k]=aptrb1[k];
	 
	 for(int k=0;k!=aptre1.size();++k) aptre[k]=aptre1[k];
	 
	
	for(int k=0;k!=NUMVAR;++k) {
		c[k]=c1[k];
		bkx[k]=bkx1[k];
		blx[k]=blx1[k];
		bux[k]=bux1[k];
		
	}
	
	
	for(int k=0;k!=aval1.size();++k) aval[k]=aval1[k];
	for(int k=0;k!=asub1.size();++k) asub[k]=asub1[k];
	
		
	for(int k=0;k!=NUMQNZ;++k){
		 qsubi[k]=qsubi1[k];
		 qsubj[k]=qsubj1[k];
		 qval[k]=qval1[k];
	}
	
	
  
  /* Create the mosek environment. */
  r = MSK_makeenv(&env,NULL);

  if ( r==MSK_RES_OK )
  { 
    /* Create the optimization task. */
    r = MSK_maketask(env,NUMCON,NUMVAR,&task);

    if ( r==MSK_RES_OK )
    {
    
      /* Append 'NUMCON' empty constraints.
       The constraints will initially have no bounds. */
      if ( r == MSK_RES_OK )
        r = MSK_appendcons(task,NUMCON);
  
      /* Append 'NUMVAR' variables.
       The variables will initially be fixed at zero (x=0). */
      if ( r == MSK_RES_OK )
        r = MSK_appendvars(task,NUMVAR);
  
      /* Optionally add a constant term to the objective. */
      if ( r ==MSK_RES_OK )
        r = MSK_putcfix(task,0.0);
      for(j=0; j<NUMVAR && r == MSK_RES_OK; ++j)
      {
        /* Set the linear term c_j in the objective.*/  
        if(r == MSK_RES_OK)
          r = MSK_putcj(task,j,c[j]);
  
        /* Set the bounds on variable j.
         blx[j] <= x_j <= bux[j] */
        if(r == MSK_RES_OK)
          r = MSK_putvarbound(task,
                              j,           /* Index of variable.*/
                              bkx[j],      /* Bound key.*/
                              blx[j],      /* Numerical value of lower bound.*/
                              bux[j]);     /* Numerical value of upper bound.*/
  
        /* Input column j of A */   
        if(r == MSK_RES_OK)
          r = MSK_putacol(task,
                          j,                 /* Variable (column) index.*/
                          aptre[j]-aptrb[j], /* Number of non-zeros in column j.*/
                          asub+aptrb[j],     /* Pointer to row indexes of column j.*/
                          aval+aptrb[j]);    /* Pointer to Values of column j.*/
        
      }
  
      /* Set the bounds on constraints.
         for i=1, ...,NUMCON : blc[i] <= constraint i <= buc[i] */
      for(i=0; i<NUMCON && r==MSK_RES_OK; ++i)
        r = MSK_putconbound(task,
                            i,           /* Index of constraint.*/
                            bkc[i],      /* Bound key.*/
                            blc[i],      /* Numerical value of lower bound.*/
                            buc[i]);     /* Numerical value of upper bound.*/

      if ( r==MSK_RES_OK )
      {
       
        /* Input the Q for the objective. */

        r = MSK_putqobj(task,NUMQNZ,qsubi,qsubj,qval);
      }

      if ( r==MSK_RES_OK )
      {
        MSKrescodee trmcode;
        
        /* Run optimizer */
        r = MSK_optimizetrm(task,&trmcode);

       
        if ( r==MSK_RES_OK )
        {
          MSKsolstae solsta;
          int j;
          
          MSK_getsolsta (task,MSK_SOL_ITR,&solsta);
          
          switch(solsta)
          {
            case MSK_SOL_STA_OPTIMAL:   
            case MSK_SOL_STA_NEAR_OPTIMAL:
              MSK_getxx(task,
                       MSK_SOL_ITR,    /* Request the interior solution. */
                       xx);
              
           			for(int k=0;k!=NUMVAR;++k) {
						myResults[k]=T(xx[k]);
						
						}
              
              break;
            case MSK_SOL_STA_DUAL_INFEAS_CER:
            case MSK_SOL_STA_PRIM_INFEAS_CER:
            case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
            case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:  
              //printf("Primal or dual infeasibility certificate found.\n");
              break;
              
            case MSK_SOL_STA_UNKNOWN:
              //printf("The status of the solution could not be determined.\n");
              break;
            default:
              //printf("Other solution status.");
              break;
          }
        }
        else
        {
          printf("Error while optimizing.\n");
        }
      }
    
      if (r != MSK_RES_OK)
      {
        /* In case of an error print error code and description. */      
        char symname[MSK_MAX_STR_LEN];
        char desc[MSK_MAX_STR_LEN];
        
        printf("An error occurred while optimizing.\n");     
        MSK_getcodedesc (r,
                         symname,
                         desc);
        printf("Error %s - '%s'\n",symname,desc);
      }
    }
    MSK_deletetask(&task);
  }
  MSK_deleteenv(&env);

 return myResults;

}