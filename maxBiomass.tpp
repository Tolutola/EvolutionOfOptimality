#include "ObjectiveFunctions.h"
#include "mosek.h"

template<class T>
vector<T> maxBiomass(const MSKint32t numvar, const MSKint32t numcon, vector<double> c1, vector<MSKint32t> aptrb1,vector<MSKint32t> aptre1,
		vector<MSKint32t> asub1, vector<double> aval1, vector<MSKboundkeye> bkc1, vector<double> blc1, vector<double> buc1,
		vector<MSKboundkeye> bkx1, vector<double> blx1, vector<T> bux2,int Osense) {

    vector<T> myResults(numvar);
           	
	// convert the vectors to arrays used by mosek
	// initialize the variables
	vector<double> bux1(numvar);
	for (int i=0;i!=numvar;++i) bux1[i]=double(bux2[i]);
	
	MSKint32t aptrb[aptrb1.size()], aptre[aptre1.size()], asub[asub1.size()];
	double aval[aval1.size()],blc[numcon], buc[numcon],blx[numvar],bux[numvar],c[numvar];
	MSKboundkeye bkc[numcon],bkx[numvar];
	
		
	// copy into arrays	
	 for(int k=0;k!=numcon;++k) bkc[k]=bkc1[k];
	
	 for(int k=0;k!=bkx1.size();++k) bkx[k]=bkx1[k];

	 for(int k=0;k!=aptrb1.size();++k) aptrb[k]=aptrb1[k];
	 
	 for(int k=0;k!=aptre1.size();++k) aptre[k]=aptre1[k];
	 
	
	for(int k=0;k!=numvar;++k) c[k]=c1[k];
	
	
	for(int k=0;k!=aval1.size();++k) aval[k]=aval1[k];
	for(int k=0;k!=asub1.size();++k) asub[k]=asub1[k];
	
	for(int k=0;k!=numcon;++k) blc[k]=blc1[k];
	 for(int k=0;k!=numcon;++k) buc[k]=buc1[k];
	
	 for(int k=0;k!=numvar;++k) blx[k]=blx1[k];

	for(int k=0;k!=numvar;++k) bux[k]=bux1[k];
	 
	
	MSKenv_t     env  = NULL;
	MSKtask_t    task = NULL;
	MSKrescodee  r;
	MSKint32t    i,j;

	/* Create the mosek environment. */
	r = MSK_makeenv(&env,NULL);
	
	
	if ( r==MSK_RES_OK )

	{
		
		/* Create the optimization task. */
		r = MSK_maketask(env,numcon,numvar,&task);

		/* Directs the log task stream to the 'printstr' function. */
		if ( r==MSK_RES_OK )
			
			/* Append 'numcon' empty constraints.
     The constraints will initially have no bounds. */
			if ( r == MSK_RES_OK )
				r = MSK_appendcons(task,numcon);

		/* Append 'numvar' variables.
     The variables will initially be fixed at zero (x=0). */
		if ( r == MSK_RES_OK )
			r = MSK_appendvars(task,numvar);
       
		for(j=0; j<numvar && r == MSK_RES_OK; ++j)
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
       for i=1, ...,numcon : blc[i] <= constraint i <= buc[i] */
		for(i=0; i<numcon && r==MSK_RES_OK; ++i)
			r = MSK_putconbound(task,
					i,           /* Index of constraint.*/
					bkc[i],      /* Bound key.*/
					blc[i],      /* Numerical value of lower bound.*/
					buc[i]);     /* Numerical value of upper bound.*/

		/* Maximize objective function. */
		if (r == MSK_RES_OK)
			if (Osense==1)
				r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);
			else
				r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);

		if ( r==MSK_RES_OK )
		{
			MSKrescodee trmcode;

			/* Run optimizer */
			r = MSK_optimizetrm(task,&trmcode);

			/* Print a summary containing information
         about the solution for debugging purposes. */
			//MSK_solutionsummary (task,MSK_STREAM_LOG);

			if ( r==MSK_RES_OK )
			{
				MSKsolstae solsta;

				if ( r==MSK_RES_OK )
					r = MSK_getsolsta (task,
							MSK_SOL_BAS,
							&solsta);
				switch(solsta)
				{
				case MSK_SOL_STA_OPTIMAL:   
				case MSK_SOL_STA_NEAR_OPTIMAL:
				{
					double *xx = (double*) calloc(numvar,sizeof(double));
					if ( xx )

					{
						MSK_getxx(task,
								MSK_SOL_BAS,    /* Request the basic solution. */
								xx);
								
									
					for(int k=0;k!=numvar;++k) {
						myResults[k]=T(xx[k]);
						
						}
					//cout<<myResults[848]<<endl;

					free(xx);
						
					}
					else 
						r = MSK_RES_ERR_SPACE;

					break;
				}
				case MSK_SOL_STA_DUAL_INFEAS_CER:
				case MSK_SOL_STA_PRIM_INFEAS_CER:
				case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
				case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:  
					//printf("Primal or dual infeasibility certificate found.\n");
					break;
				case MSK_SOL_STA_UNKNOWN:
				{
					char symname[MSK_MAX_STR_LEN];
					char desc[MSK_MAX_STR_LEN];

					/* If the solutions status is unknown, print the termination code
               indicating why the optimizer terminated prematurely. */

					MSK_getcodedesc(trmcode,
							symname,
							desc);

					printf("The solution status is unknown.\n");
					printf("The optimizer terminitated with code: %s\n",symname);
					break;
				}
				default:
					printf("Other solution status.\n");
					break;
				}
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

		/* Delete the task and the associated data. */
		MSK_deletetask(&task);
	}

	return myResults;
	/* Delete the environment and the associated data. */
	MSK_deleteenv(&env);
}