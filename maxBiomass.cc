#include "ObjectiveFunctions.h"
#include "mosek.h"


vector<double> maxBiomass(const MSKint32t numvar, const MSKint32t numcon, vector<double> c, vector<MSKint32t> aptrb1,vector<MSKint32t> aptre1,
		vector<MSKint32t> asub1, vector<double> aval1, vector<MSKboundkeye> bkc1, vector<double> blc1, vector<double> buc1,
		vector<MSKboundkeye> bkx1, vector<double> blx1, vector<double> bux1) {

//	const MSKint32t numvar = 4,
//			numcon = 3;
	
//	double       c[]     = {3.0, 1.0, 5.0, 1.0};
	/* Below is the sparse representation of the A
     matrix stored by column. */
//	MSKint32t    aptrb[] = {0, 2, 5, 7},
//			aptre[] = {2, 5, 7, 9},
//			asub[]  = { 0, 1,
//					0, 1, 2,
//					0, 1,
//					1, 2};
//	double       aval[]  = { 3.0, 2.0,
//			1.0, 1.0, 2.0,
//			2.0, 3.0,
//			1.0, 3.0};
	
	// convert the vectors to arrays used by mosek
	// initialize the variables
	cout<<"inside"<<numvar<<endl;
	
	MSKint32t aptrb[numvar], aptre[numvar], asub[numvar];
	double aval[aval1.size()],blc[numcon], buc[numcon],blx[numvar],bux[numvar];
	MSKboundkeye bkc[numcon],bkx[numvar];
	
	
	
	// copy into arrays	
	 memcpy( bkc, &bkc1[0], sizeof( MSKboundkeye ) * bkc1.size() );
	 memcpy( bkx, &bkx1[0], sizeof( MSKboundkeye ) * bkx1.size() );
	 memcpy( aptrb, &aptrb1[0], sizeof( MSKint32t ) * aptrb1.size() );
	 memcpy( aptre, &aptre1[0], sizeof( MSKint32t ) * aptre1.size() );
	 memcpy( asub, &asub1[0], sizeof( MSKint32t ) * asub1.size() );
	 memcpy( aval, &aval1[0], sizeof( double ) * aval1.size() );
	 memcpy( blc, &blc1[0], sizeof( double ) * blc1.size() );
	 memcpy( buc, &buc1[0], sizeof( double ) * buc1.size() );
	 memcpy( blx, &blx1[0], sizeof( double ) * blx1.size() );
	 memcpy( bux, &bux1[0], sizeof( double ) * bux1.size() );

	
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
			// r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);

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
			r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MAXIMIZE);

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
						
					}
					else 
						r = MSK_RES_ERR_SPACE;

					break;
				}
				case MSK_SOL_STA_DUAL_INFEAS_CER:
				case MSK_SOL_STA_PRIM_INFEAS_CER:
				case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
				case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:  
					printf("Primal or dual infeasibility certificate found.\n");
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

	/* Delete the environment and the associated data. */
	MSK_deleteenv(&env);

	return xx;
}