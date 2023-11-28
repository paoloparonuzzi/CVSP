

#include "CLIQUE.h"

#define number_of_CPU 1


///////////////////////////////////////////////////////////////////////////////
//#define CPLEX_OUTPUT

//#define print_clique_model

//#define print_sol_clique
///////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////
CPXENVptr env_clique;
CPXLPptr  lp_clique;

int stat_mipopt;
int status_mipopt;
int status_getbestobjval;
int status_getobjval;
//int status;

double cplex_ub,cplex_lb,cplex_lp;
int cplex_nodes;
double* cplex_sol;

//int ccnt;
//double* obj;
//double* rhs;
//char* sense;
//int* rmatbeg, *rmatind;
//int* cmatbeg, *cmatind;
//double* rmatval,* cmatval, *lb, *ub;
//int nzcnt;
//int rcnt;
char* _ctype;
/////////////////////////////////////////////////////////////////////////////


/***************************************************************************/
void clique_free_cplex (graphFF G)
/***************************************************************************/
{

	CPXfreeprob(env_clique,&lp_clique);

	status = CPXcloseCPLEX (&env_clique);
}


/***************************************************************************/
void clique_load_cplex (graphFF G_bar)
/***************************************************************************/
{

	env_clique=CPXopenCPLEX(&status);
	if (status!=0)
	{
		printf("cannot open CPLEX environment \n ");
		exit(-1);

	}
	lp_clique=CPXcreateprob(env_clique,&status,"clique");
	if (status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}
	CPXchgobjsen(env_clique,lp_clique,CPX_MAX);

	///////////////////////////////////////////////////
	CPXsetintparam (env_clique, CPX_PARAM_THREADS, number_of_CPU);
	///////////////////////////////////////////////////

	ccnt=G_bar->n;
	obj=(double*) calloc(ccnt,sizeof(double));
	lb=(double*) calloc(ccnt,sizeof(double));
	ub=(double*) calloc(ccnt,sizeof(double));
	_ctype=(char*) calloc(ccnt,sizeof(char));


	char **colname=(char**) calloc(ccnt,sizeof(char*));
	for(int i=0;i<ccnt;i++){colname[i]=(char*) calloc(1000,sizeof(char));}


	for(int i=0; i<ccnt; i++)
	{
		obj[i]=G_bar->W[i];
		lb[i]=0.0;
		ub[i]=1.0;
		sprintf(colname[i], "x%d",i);

		_ctype[i]='B';
	}

	status=CPXnewcols(env_clique,lp_clique,ccnt,obj,lb,ub,_ctype,colname);
	if(status!=0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	for(int i=0;i<ccnt;i++){free(colname[i]);}
	free(colname);




	for(int i=0; i<G_bar->m; i++)
	{


		rcnt=1;
		nzcnt=2;
		rhs=(double*) calloc(rcnt,sizeof(double));
		sense=(char*) calloc(rcnt,sizeof(double));

		rhs[0]=1.0;
		sense[0]='L';

		rmatbeg=(int*) calloc(rcnt,sizeof(int));
		rmatind=(int*) calloc(nzcnt,sizeof(int));
		rmatval=(double*) calloc(nzcnt,sizeof(double));


		rmatval[0]=1.0;
		rmatind[0]=G_bar->T[i];
		rmatval[1]=1.0;
		rmatind[1]=G_bar->H[i];

		rmatbeg[0]=0;

		status=CPXaddrows(env_clique,lp_clique,0,rcnt,nzcnt,rhs,sense,rmatbeg,rmatind,rmatval,NULL,NULL);
		if (status!=0)
		{
			printf("error in CPXaddrows\n");
			exit(-1);
		}
		free(rmatbeg);
		free(rmatval);
		free(rmatind);
		free(rhs);
		free(sense);

	}

#ifdef print_clique_model
	cout << "printing LP\n";
	//cin.get();
	status=CPXwriteprob(env_clique,lp_clique,"clique.lp",NULL);
	if(status!=0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
#endif

	free(obj);
	free(lb);
	free(ub);
	free(_ctype);
}


/***************************************************************************/
void clique_update_cplex (int *clique,graphFF G,graphFF G_bar)
/***************************************************************************/
{

//	status=CPXwriteprob(env_clique,lp_clique,"clique_round_before.lp",NULL);
//	if(status!=0)
//	{
//		printf("error in CPXwriteprob\n");
//		exit(-1);
//	}
//	cout << "printing LP clique\n";
//	cin.get();



	for(int i=0; i<G->n; i++)
	{
		for(int j=i+1; j<G->n; j++)
		{

			if(clique[i]>0.5 && clique[j]>0.5){
				rcnt=1;
				nzcnt=2;
				rhs=(double*) calloc(rcnt,sizeof(double));
				sense=(char*) calloc(rcnt,sizeof(double));

				rhs[0]=1.0;
				sense[0]='L';

				rmatbeg=(int*) calloc(rcnt,sizeof(int));
				rmatind=(int*) calloc(nzcnt,sizeof(int));
				rmatval=(double*) calloc(nzcnt,sizeof(double));


				rmatval[0]=1.0;
				rmatind[0]=i;
				rmatval[1]=1.0;
				rmatind[1]=j;

				rmatbeg[0]=0;

				status=CPXaddrows(env_clique,lp_clique,0,rcnt,nzcnt,rhs,sense,rmatbeg,rmatind,rmatval,NULL,NULL);
				if (status!=0)
				{
					printf("error in CPXaddrows\n");
					exit(-1);
				}
				free(rmatbeg);
				free(rmatval);
				free(rmatind);
				free(rhs);
				free(sense);
			}

		}
	}

//	status=CPXwriteprob(env_clique,lp_clique,"clique_round_after.lp",NULL);
//	if(status!=0)
//	{
//		printf("error in CPXwriteprob\n");
//		exit(-1);
//	}
//	cout << "printing LP clique\n";
//	cin.get();


}

void clique_fix_node (int v)
/***************************************************************************/
{
	int cnt=1;
	int indices = v;
	char lu='B';
	double bd=1;
	
	status = CPXchgbds(env_clique,lp_clique, cnt, &indices, &lu, &bd);

	
	//status=CPXwriteprob(env_clique,lp_clique,"clique_round_before.lp",NULL);
	if(status!=0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
//	cout << "printing LP clique\n";
//	cin.get();


}

void clique_remove_node (int v)
/***************************************************************************/
{
	int cnt=1;
	int indices = v;
	char lu='B';
	double bd=0;
	
	status = CPXchgbds(env_clique,lp_clique, cnt, &indices, &lu, &bd);

	
	//status=CPXwriteprob(env_clique,lp_clique,"clique_round_before.lp",NULL);
	if(status!=0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
//	cout << "printing LP clique\n";
//	cin.get();

}

void clique_free_node (int v)
/***************************************************************************/
{
	int cnt=1;
	int indices = v;
	char lu='U';
	double bd=1;
	
	status = CPXchgbds(env_clique,lp_clique, cnt, &indices, &lu, &bd);

	lu='L';
	bd=0;
	
	status = CPXchgbds(env_clique,lp_clique, cnt, &indices, &lu, &bd);

	
	//status=CPXwriteprob(env_clique,lp_clique,"clique_round_before.lp",NULL);
	if(status!=0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
//	cout << "printing LP clique\n";
//	cin.get();


}


/***************************************************************************/
int clique_solve_cplex(int *solution,int n)
/***************************************************************************/
{


#ifdef CPLEX_OUTPUT
	CPXsetintparam (env_clique, CPX_PARAM_SCRIND, CPX_ON);
#endif

	status_mipopt=CPXmipopt(env_clique,lp_clique);
	if(status_mipopt!=0)
	{
		printf("error in CPXmipopt\n");;
		exit(-1);
	}


	/////////////////////////////////////////////////////////////////////////
	double *solution_CPLEX=new double[n];

	status=CPXgetmipx(env_clique,lp_clique,solution_CPLEX,0,n-1);
	if(status!=0)
	{
		exit(-1);
	}


	for(int i=0; i<n; i++)
	{
		if(solution_CPLEX[i]>0.5){
			solution[i]=1;
		}
		else{
			solution[i]=0;
		}
	}

	delete []solution_CPLEX;
	/////////////////////////////////////////////////////////////////////////


	status=CPXgetmipobjval(env_clique,lp_clique,&cplex_lb);
	if(status!=0)
	{
		cout << "problem in CPXgetmipobjval";
		exit(-1);
	}


	//cout << "OBJ_CLIQUE ->\t " << cplex_lb << endl;

//	status=CPXwriteprob(env_clique,lp_clique,"clique_fixing.lp",NULL);
//	if(status!=0) {printf("error in CPXwriteprob\n");exit(-1);}
//	cout << "PRINTED AFTER OPTIMIZATION\n";
//	cin.get();

	return (int)(cplex_lb+0.5);
}



