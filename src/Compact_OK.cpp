#define _CRT_SECURE_NO_WARNINGS

#include "Compact.h"

#define number_of_CPU 1
#define print_sol_details_lp


//to use the callbacks
//#define BC
// int vertex_number; //this must correspond to G->n


typedef struct{
	graphFF G;
	int k;}
parametresCplex;

/***************************************************************************/
int position_x_v_i(int v, int i)
/***************************************************************************/
{
	return v * k + i;
}

/***************************************************************************/
int CPXPUBLIC CPXsetusercutcallbackfunc(CPXCENVptr env,void *cbdata,int wherefrom,void *cbhandle,int *useraction_p)
/***************************************************************************/
{

	(*useraction_p)=CPX_CALLBACK_SET;

	int status;
	int numvar=k*vertex_number;
	int i;

	double *xstar=new double[numvar];
	status=CPXgetcallbacknodex(env,cbdata,wherefrom,xstar,0,numvar-1);
	if(status!=0) {printf("cannot get the x\n");cin.get();exit(-1);}
	cout << "\nCurrent Solution\n";for (i=0;i<numvar;i++) {cout << xstar[i] << " ";}cout << endl;

	double obj_value;
	status=CPXgetcallbacknodeobjval(env,cbdata,wherefrom,&obj_value);
	if(status!=0) {printf("cannot get the x\n");cin.get();exit(-1);}
	cout << "Current Value: ->\t" << obj_value << endl;

	cin.get();

	//add your favourite cat
	//	bool cut_found=false;
	//
	//	if(cut_found){
	//
	//		int nzcnt=??;
	//
	//		int *cut_matind=(int*) calloc(nzcnt,sizeof(int));
	//		double *cut_rmatval=(double*) calloc(nzcnt,sizeof(double));
	//
	//		for (int f=0 ; f<nzcnt ; f++){
	//			cut_rmatval[f]=??; //value
	//			cut_matind[f]=??; //index of the variable
	//		}
	//		status=CPXcutcallbackadd (env,cbdata,wherefrom,nzcnt,(double)1.0,'G',cut_matind,cut_rmatval,0);
	//		if(status!=0) {printf("cannot add the cut\n");exit(-1);}
	//
	//		free(cut_rmatval);free(cut_matind);
	//	}

	delete[] xstar;

	return 0;

}//mycutcallback()


/***************************************************************************/
int CPXPUBLIC CPXsetlazycutcallbackfunc(CPXCENVptr env,void *cbdata,int wherefrom,void *temp,int *useraction_p)
/***************************************************************************/
{
	parametresCplex* argtemp=(parametresCplex*)temp;
	graphFF G=argtemp->G;
	int k=argtemp->k;
	(*useraction_p)=CPX_CALLBACK_SET;
	cout<<" instance "<<G->n<<" "<<k<<endl;
	int status;
	int numvar=k*G->n;
	int i;

	double *xstar=new double[numvar];
	status=CPXgetcallbacknodex(env,cbdata,wherefrom,xstar,0,numvar-1);
	if(status!=0) {printf("cannot get the x\n");cin.get();exit(-1);}
	cout << "\nCurrent Solution\n";
	for (i=0;i<numvar;i++) {cout << xstar[i] << " - ";}cout << endl;

	for (int v = 0; v < G->n; v++) {
		for (int u = 0; u < G->n; u++) {
			if(G->AMatrix[u][v] == 1 || G->AMatrix[v][u] == 1) {
				int pos1=-1;
				for (int i = 0; i < k; i++) {
					if(xstar[v*k+i]==1) {
						pos1=i;
						break;
					}
				}
				if(pos1==-1)continue;
				int pos2=-1;
				for (int i = 0; i < k; i++) {
					if(xstar[u*k+i]==1) {
						pos2=i;
						break;
					}
				}
				if(pos2==-1)continue;

			}
		}
	}
	double obj_value;
	status=CPXgetcallbacknodeobjval(env,cbdata,wherefrom,&obj_value);
	if(status!=0) {printf("cannot get the x\n");cin.get();exit(-1);}
	cout << "Current Value: ->\t" << obj_value << endl;

	cin.get();

	//add your favourite cat
	//	bool cut_found=false;
	//
	//	if(cut_found){
	//
	//		int nzcnt=??;
	//
	//		int *cut_matind=(int*) calloc(nzcnt,sizeof(int));
	//		double *cut_rmatval=(double*) calloc(nzcnt,sizeof(double));
	//
	//		for (int f=0 ; f<nzcnt ; f++){
	//			cut_rmatval[f]=??; //value
	//			cut_matind[f]=??; //index of the variable
	//		}
	//		status=CPXcutcallbackadd (env,cbdata,wherefrom,nzcnt,(double)1.0,'G',cut_matind,cut_rmatval,0);
	//		if(status!=0) {printf("cannot add the cut\n");exit(-1);}
	//
	//		free(cut_rmatval);free(cut_matind);
	//	}

	delete[] xstar;

	return 0;

}//mycutcallback()


/***************************************************************************/
void compact_model_free()
/***************************************************************************/
{

	//////////////////////////////////////////////////////////////
	status = CPXfreeprob(env_COMPACT, &lp_COMPACT);
	if (status != 0) {
		cout << "error in CPXfreeprob\n";
		exit(-1);
	}
	//////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////
	status = CPXcloseCPLEX(&env_COMPACT);
	if (status != 0) {
		cout << "cannot close CPLEX environment\n";
		exit(-1);
	}

	////////////////////////////////////////////////////////////////////

}

/***************************************************************************/
int * separate_clique_inequalities(graphFF G, int k, int q, int b,
		double ** vector_xiv)
/***************************************************************************/
{

	int * solution = (int*) malloc(sizeof(int) * k);
	int * vertices_in_clique = (int*) malloc(sizeof(int) * G->n);
	int * set_taken = (int*) malloc(sizeof(int) * k);

	for (int ik = 0; ik < k; ik++) {
		set_taken[ik] = 0;
	}
	for (int v = 0; v < G->n; v++) {
		vertices_in_clique[v] = 0;
	}
	for (int ik = 0; ik < k; ik++) {
		int max_k = -1;
		int max_v = -1;
		int max_value = -1;
		for (int ikmax = 0; ikmax < k; ikmax++) {
			if (set_taken[ikmax] == 0) {
				for (int v = 0; v < G->n; v++) {
					int test = 1;
					for (int u = 0; u < G->n; u++) {
						if (vertices_in_clique[u] == 1 && G->AMatrix[u][v] == 0
								&& G->AMatrix[v][u] == 0) {
							test = 0;
							break;
						}
					}
					if (test && vector_xiv[ikmax][v] > max_value) {
						max_k = ikmax;
						max_v = v;
						max_value = vector_xiv[ikmax][v];
					}
				}
			}
		}
		solution[max_k] = max_v;
	}
	free(set_taken);
	free(vertices_in_clique);

	return solution;
}

/***************************************************************************/
int **clique_maximal_heuristic(graphFF G, int k, int q, int b)
/***************************************************************************/
{
	//sort by degree
	int *order = (int*) malloc(sizeof(int) * G->n);

	int *flag_order = (int*) malloc(sizeof(int) * G->n);
	for (int u = 0; u < G->n; u++) {
		flag_order[u] = 0;
	}
	for (int u = 0; u < G->n; u++) {
		int max = -1;
		int posmax = -1;
		for (int v = 0; v < G->n; v++) {
			if (flag_order == 0 && G->DT[v] > max) {
				max = G->DT[v];
				posmax = v;
			}
		}
		flag_order[posmax] = 1;
		order[u] = posmax;
	}
	free(flag_order);
	int *node_in_clique = (int*) malloc(sizeof(int) * G->n);

	for (int i = 0; i < G->n; i++) {
		node_in_clique[i] = 0;
	}
	int **cliques = (int**) malloc(sizeof(int*) * k);
	for (int i = 0; i < k; i++) {
		cliques[i] = (int*) malloc(sizeof(int) * G->n);
	}
	for (int ik = 0; ik < k; ik++) {
		for (int u = 0; u < G->n; u++) {
			if (node_in_clique[u] == 0) {
				int test = 1;
				for (int v = 0; v < G->n; v++) {
					if (cliques[ik][v] == 1 && G->AMatrix[u][v] == 0
							&& G->AMatrix[v][u] == 0) {
						test = 0;
						break;
					}
				}
				if (test) {
					cliques[ik][u] = 1;
					node_in_clique[u] = 1;
				}
			}
		}
	}
	free(node_in_clique);
	return cliques;
}

/***************************************************************************/
int ** clique_maximum_exact(graphFF G, int k, int q, int b)
/***************************************************************************/
{

	int *node_in_clique = (int*) malloc(sizeof(int) * G->n);

	for (int i = 0; i < G->n; i++) {
		node_in_clique[i] = 0;
	}

	int **cliques = (int**) malloc(sizeof(int*) * k);
	for (int i = 0; i < k; i++) {
		cliques[i] = (int*) malloc(sizeof(int) * G->n);
	}

	for (int ik = 0; ik < k; ik++) {
		CPXENVptr env_CLIQUE = CPXopenCPLEX(&status);
		if (status != 0) {
			printf("cannot open CPLEX environment\n");
			exit(-1);
		}
		CPXLPptr lp_CLIQUE = CPXcreateprob(env_CLIQUE, &status, "K-SEP");
		if (status != 0) {
			printf("cannot create problem\n");
			exit(-1);
		}

		CPXchgobjsen(env_CLIQUE, lp_CLIQUE, CPX_MAX);

		ccnt = G->n;
		obj = (double*) calloc(ccnt, sizeof(double));
		lb = (double*) calloc(ccnt, sizeof(double));
		ub = (double*) calloc(ccnt, sizeof(double));
		c_type = (char*) calloc(ccnt, sizeof(char));

		char **colname = new char*[ccnt];
		for (int i = 0; i < ccnt; i++) {
			colname[i] = new char[100];
		}

		//creating variables
		int dummy = 0;
		for (int v = 0; v < G->n; v++) {
			obj[dummy] = 1.0;
			lb[dummy] = 0.0;
			if (node_in_clique[dummy] == 1)
				ub[dummy] = 0.0;
			else
				ub[dummy] = 1.0;
			c_type[dummy] = 'B';
			sprintf(colname[dummy], "x(%d).%d", v, dummy);
			dummy++;
		}
		//routine cplex per creare le colonne
		status = CPXnewcols(env_CLIQUE, lp_CLIQUE, ccnt, obj, lb, ub, c_type,
				colname);
		if (status != 0) {
			printf("error in CPXnewcols\n");
			exit(-1);
		}

		for (int i = 0; i < ccnt; i++) {
			delete[] colname[i];
		}
		delete[] colname;

		free(obj);
		free(lb);
		free(ub);
		free(c_type);

		//For each no edge in G ...


		for (int u = 0; u < G->n - 1; u++) {
			for (int v = u + 1; v < G->n; v++) {
				if (G->AMatrix[u][v] == 0 && G->AMatrix[v][u] == 0) {
					rcnt = 1;
					nzcnt = k;

					rhs = (double*) calloc(rcnt, sizeof(double));
					sense = (char*) calloc(rcnt, sizeof(double));

					rhs[0] = 1.0;
					sense[0] = 'L';

					rmatbeg = (int*) calloc(rcnt, sizeof(int));
					rmatind = (int*) calloc(nzcnt, sizeof(int));
					rmatval = (double*) calloc(nzcnt, sizeof(double));

					rmatval[0] = 1.0;
					rmatind[0] = u;
					rmatval[1] = 1.0;
					rmatind[1] = v;

					//for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
					//for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

					rmatbeg[0] = 0;
					rmatbeg[1] = 0;

					status = CPXaddrows(env_CLIQUE, lp_CLIQUE, 0, rcnt, nzcnt,
							rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
					if (status != 0) {
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
		status = CPXlpopt(env_CLIQUE, lp_CLIQUE);
		if (status != 0) {
			cout << "err_FILEor in CPXlpopt slave solve\n";
			exit(-1);
		}

		double objval_lp;

		status = CPXgetobjval(env_CLIQUE, lp_CLIQUE, &objval_lp);
		if (status != 0) {
			printf("error in CPXgetmipobjval\n");
			exit(-1);
		}

		cout << "\n\n Clique size ->\t " << objval_lp << endl;
		double* x_exact;
		x_exact = (double*) calloc(G->n, sizeof(double));

		//getting the solution
		status = CPXgetmipx(env_CLIQUE, lp_CLIQUE, x_exact, 0, G->n - 1);
		if (status != 0) {
			cout << endl << "\nNO SOLUTION!!" << endl;
			exit(-1);
			//err_FILE <<"error in CPXgetmipx model_solve\n";
		}

		cout << "\n\nCURRENT VARIABLES:\n";
		for (int v = 0; v < G->n; v++) {
			if (x_exact[v] > 0.9) {
				cliques[ik][v] = 1;
				node_in_clique[v] = 1;
			} else
				cliques[ik][v] = 0;

		}

	}
	free(node_in_clique);
	return cliques;
}

/***************************************************************************/
int * neighboors(graphFF G, int k, int q, int b, int u)
/***************************************************************************/
{
	int* neigh = (int*) calloc(G->n, sizeof(double));
	for (int v = 0; v < G->n; v++) {
		if (G->AMatrix[u][v] == 1 || G->AMatrix[v][u] == 1) {
			neigh[v] = 1;
		} else {
			neigh[v] = 0;
		}

	}
	neigh[u] = 1;

	return neigh;

}

/***************************************************************************/
int stable_maximum_exact(graphFF G, int k, int q, int b, int *subgraph)
/***************************************************************************/
{

	CPXENVptr env_STABLE = CPXopenCPLEX(&status);
	if (status != 0) {
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}
	CPXLPptr ILP_STABLE = CPXcreateprob(env_STABLE, &status, "K-SEP");
	if (status != 0) {
		printf("cannot create problem\n");
		exit(-1);
	}

	CPXchgobjsen(env_STABLE, ILP_STABLE, CPX_MAX);

	ccnt = G->n;
	obj = (double*) calloc(ccnt, sizeof(double));
	lb = (double*) calloc(ccnt, sizeof(double));
	ub = (double*) calloc(ccnt, sizeof(double));
	c_type = (char*) calloc(ccnt, sizeof(char));

	char **colname = new char*[ccnt];
	for (int i = 0; i < ccnt; i++) {
		colname[i] = new char[100];
	}

	//creating variables
	int dummy = 0;
	for (int v = 0; v < G->n; v++) {
		obj[dummy] = 1.0;
		lb[dummy] = 0.0;
		if (subgraph[dummy] == 1)
			ub[dummy] = 1.0;
		else
			ub[dummy] = 0.0;
		c_type[dummy] = 'B';
		sprintf(colname[dummy], "x(%d).%d", v, dummy);
		dummy++;
	}
	//routine cplex per creare le colonne
	status = CPXnewcols(env_STABLE, ILP_STABLE, ccnt, obj, lb, ub, c_type,
			colname);
	if (status != 0) {
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	for (int i = 0; i < ccnt; i++) {
		delete[] colname[i];
	}
	delete[] colname;

	free(obj);
	free(lb);
	free(ub);
	free(c_type);

	//	//For each no edge in G ...


	for (int u = 0; u < G->n - 1; u++) {
		if (subgraph[u] == 1) {
			for (int v = u + 1; v < G->n; v++) {
				if (subgraph[v] == 1) {
					if (G->AMatrix[u][v] == 1 || G->AMatrix[v][u] == 1) {
						rcnt = 1;
						nzcnt = 2;

						rhs = (double*) calloc(rcnt, sizeof(double));
						sense = (char*) calloc(rcnt, sizeof(double));

						rhs[0] = 1.0;
						sense[0] = 'L';

						rmatbeg = (int*) calloc(rcnt, sizeof(int));
						rmatind = (int*) calloc(nzcnt, sizeof(int));
						rmatval = (double*) calloc(nzcnt, sizeof(double));

						rmatval[0] = 1.0;
						rmatind[0] = u;
						rmatval[1] = 1.0;
						rmatind[1] = v;

						//for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
						//for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

						rmatbeg[0] = 0;
						rmatbeg[1] = 0;

						status = CPXaddrows(env_STABLE, ILP_STABLE, 0, rcnt,
								nzcnt, rhs, sense, rmatbeg, rmatind, rmatval,
								NULL, NULL);
						if (status != 0) {
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
		}
	}

	status = CPXmipopt(env_STABLE, ILP_STABLE);
	if (status != 0) {
		cout << "err_FILEor in CPXlpopt slave solve\n";
		exit(-1);
	}
	/*	cout << "\n\nwrite compact model\n";
	 char *_name=new char[100];
	 sprintf(_name,"model_base.lp");

	 status=CPXwriteprob(env_STABLE, ILP_STABLE,_name,NULL);





	 ////////////////////////////////////
	 double objval_lp;

	 status=CPXgetobjval(env_STABLE,ILP_STABLE,&objval_lp);
	 if(status!=0)
	 {
	 printf("error in CPXgetmipobjval\n");
	 exit(-1);
	 }

	 cout << "\n\n stable size ->\t " << objval_lp << endl;
	 double* x_exact;
	 x_exact=(double*) calloc(G->n,sizeof(double));

	 //getting the solution
	 status=CPXgetmipx(env_STABLE,ILP_STABLE,x_exact,0,G->n-1);
	 if(status!=0)
	 {
	 cout << endl << "\nNO SOLUTION!!" << endl;
	 exit(-1);
	 //err_FILE <<"error in CPXgetmipx model_solve\n";
	 }



	 cout << "\n\nCURRENT VARIABLES:\n";
	 for(int v=0; v< G->n; v++)
	 {
	 if(x_exact[v]>0.9){
	 cout<<v<<" ";
	 }

	 }
	 cout<<endl;
	 ///////////////////////////////////

	 //*/
	double objval_lp;

	status = CPXgetobjval(env_STABLE, ILP_STABLE, &objval_lp);
	if (status != 0) {
		printf("error in CPXgetmipobjval\n");
		exit(-1);
	}
	free(subgraph);
	return (int) (objval_lp + 0.1);
}

/***************************************************************************/
void maximize_clique(graphFF G, int* sol)
/***************************************************************************/
{
	for (int u = 0; u < G->n; u++) {
		if (sol[u] == 0) {
			int add = 1;
			for (int v = 0; v < G->n; v++) {
				if (sol[v] == 1) {
					if (G->AMatrix[u][v] == 0 && G->AMatrix[v][u] == 0) {
						add = 0;
						break;
					}

				}
			}
			if (add) {
				sol[u] = 1;
			}
		}
	}

}

/***************************************************************************/
int clique_maximum_exact(graphFF G, int k, int q, int b, int *subgraph, int* sol)
/***************************************************************************/
{

	CPXENVptr env_CLIQUE = CPXopenCPLEX(&status);
	if (status != 0) {
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}
	CPXLPptr ILP_CLIQUE = CPXcreateprob(env_CLIQUE, &status, "K-SEP");
	if (status != 0) {
		printf("cannot create problem\n");
		exit(-1);
	}

	CPXchgobjsen(env_CLIQUE, ILP_CLIQUE, CPX_MAX);

	ccnt = G->n;
	obj = (double*) calloc(ccnt, sizeof(double));
	lb = (double*) calloc(ccnt, sizeof(double));
	ub = (double*) calloc(ccnt, sizeof(double));
	c_type = (char*) calloc(ccnt, sizeof(char));

	char **colname = new char*[ccnt];
	for (int i = 0; i < ccnt; i++) {
		colname[i] = new char[100];
	}

	//creating variables
	int dummy = 0;
	for (int v = 0; v < G->n; v++) {
		obj[dummy] = 1.0;
		lb[dummy] = 0.0;
		if (subgraph[dummy] == 1)
			ub[dummy] = 1.0;
		else
			ub[dummy] = 0.0;
		c_type[dummy] = 'B';
		sprintf(colname[dummy], "x(%d).%d", v, dummy);
		dummy++;
	}
	//routine cplex per creare le colonne
	status = CPXnewcols(env_CLIQUE, ILP_CLIQUE, ccnt, obj, lb, ub, c_type,
			colname);
	if (status != 0) {
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	for (int i = 0; i < ccnt; i++) {
		delete[] colname[i];
	}
	delete[] colname;

	free(obj);
	free(lb);
	free(ub);
	free(c_type);

	//	//For each no edge in G ...


	for (int u = 0; u < G->n - 1; u++) {
		if (subgraph[u] == 1) {
			for (int v = u + 1; v < G->n; v++) {
				if (subgraph[v] == 1) {
					if (G->AMatrix[u][v] == 0 && G->AMatrix[v][u] == 0) {
						rcnt = 1;
						nzcnt = 2;

						rhs = (double*) calloc(rcnt, sizeof(double));
						sense = (char*) calloc(rcnt, sizeof(double));

						rhs[0] = 1.0;
						sense[0] = 'L';

						rmatbeg = (int*) calloc(rcnt, sizeof(int));
						rmatind = (int*) calloc(nzcnt, sizeof(int));
						rmatval = (double*) calloc(nzcnt, sizeof(double));

						rmatval[0] = 1.0;
						rmatind[0] = u;
						rmatval[1] = 1.0;
						rmatind[1] = v;

						//for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
						//for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

						rmatbeg[0] = 0;
						rmatbeg[1] = 0;

						status = CPXaddrows(env_CLIQUE, ILP_CLIQUE, 0, rcnt,
								nzcnt, rhs, sense, rmatbeg, rmatind, rmatval,
								NULL, NULL);
						if (status != 0) {
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
		}
	}

	status = CPXmipopt(env_CLIQUE, ILP_CLIQUE);
	if (status != 0) {
		cout << "err_FILEor in CPXlpopt slave solve\n";
		exit(-1);
	}
	/*	cout << "\n\nwrite compact model\n";
	 char *_name=new char[100];
	 sprintf(_name,"model_base.lp");

	 status=CPXwriteprob(env_STABLE, ILP_STABLE,_name,NULL);





	 ////////////////////////////////////
	 double objval_lp;

	 status=CPXgetobjval(env_STABLE,ILP_STABLE,&objval_lp);
	 if(status!=0)
	 {
	 printf("error in CPXgetmipobjval\n");
	 exit(-1);
	 }

	 cout << "\n\n stable size ->\t " << objval_lp << endl;
	 */
	double* x_exact;
	x_exact = (double*) calloc(G->n, sizeof(double));

	//getting the solution
	status = CPXgetmipx(env_CLIQUE, ILP_CLIQUE, x_exact, 0, G->n - 1);
	if (status != 0) {
		free(x_exact);
		cout << endl << "\nNO SOLUTION!!" << endl;
		exit(-1);
		//err_FILE <<"error in CPXgetmipx model_solve\n";
	}

	//	cout << "\n\nCURRENT VARIABLES:\n";
	for (int v = 0; v < G->n; v++) {
		if (x_exact[v] > 0.9) {
			//			cout<<v<<" ";
			sol[v] = 1;
		} else {
			sol[v] = 0;
		}

	}
	free(x_exact);
	//cout<<endl;
	///////////////////////////////////

	//*/
	double objval_lp;

	status = CPXgetobjval(env_CLIQUE, ILP_CLIQUE, &objval_lp);
	if (status != 0) {
		printf("error in CPXgetmipobjval\n");
		exit(-1);
	}
	//	free(subgraph);
	return (int) (objval_lp + 0.1);
}

/***************************************************************************/
void compact_model_load(graphFF G, int k, int q, int b)
/***************************************************************************/
{
	int option1 = 1;//ct8
	int option2 = 0;//ct9 and ct10
	int option3 = 0;//ct11
	int option4 = 0;//ct12
	int optionClique = 0;//search k-1 cliques, the first one must be assign to V_1 or V_0, the second one must be assign to V_2, v_1 or V_0 ...
	printf("Compact\n");
		//cin.get();
	int i;

	//////////////////////////////////////////////////////////////////////////////////////////////////
	env_COMPACT = CPXopenCPLEX(&status);
	if (status != 0) {
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	lp_COMPACT = CPXcreateprob(env_COMPACT, &status, "K-SEP");
	if (status != 0) {
		printf("cannot create problem\n");
		exit(-1);
	}

	int _dimension = G->n * k;

	CPXchgobjsen(env_COMPACT, lp_COMPACT, CPX_MAX);

	////////////////////////////////////////////////////////////////////////////
	CPXsetintparam(env_COMPACT, CPX_PARAM_THREADS, number_of_CPU);
	CPXsetdblparam(env_COMPACT, CPX_PARAM_TILIM, timeLimit);
	//CPXsetdblparam (env_COMPACT, CPX_PARAM_EPGAP ,1e-7);//tolleranza nella soluzione
	//CPXsetdblparam (env_COMPACT, CPX_PARAM_EPAGAP,1e-15);//tolleranza nell'albero
	//CPXsetintparam(env_COMPACT,CPX_PARAM_CLIQUES,-1);//toglie le clique
	//CPXsetintparam(env_COMPACT,CPX_PARAM_CUTPASS,-1);//toglie tutti i tagli
	//CPXsetintparam(env_COMPACT,CPX_PARAM_PREIND,0);//toglie il preprocessamento
	//CPXfreepresolve (env_COMPACT, slave_mip[clique_number]);
	////////////////////////////////////////////////////////////////////////////


#ifdef use_traditional_branch_and_cut
	status = CPXsetintparam (env_COMPACT, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
#endif
	CPXsetintparam(env_COMPACT, CPX_PARAM_SCRIND, CPX_ON);
#ifdef print_CPLEX_ON

	//CPXsetintparam (env_COMPACT, CPX_PARAM_MIPINTERVAL, 1);
#endif

#ifdef set_tolerance
	// * Set relative tolerance *
	status = CPXsetdblparam (env_COMPACT, CPX_PARAM_EPAGAP, 0.9);
	if (status)
	{
		fprintf (stderr, "error for CPX_PARAM_EPAGAP\n");
	}

	// * Set relative tolerance *
	status = CPXsetdblparam (env_COMPACT, CPX_PARAM_EPGAP, 0.000001);
	if (status)
	{
		fprintf (stderr, "error for CPX_PARAM_EPGAP\n");
	}
#endif

#ifdef no_cplex_heur
	status = CPXsetintparam (env_COMPACT, CPX_PARAM_HEURFREQ , -1);
#endif

#ifdef INPUT_CPX_policy
	// * Set MIP emphasis *
	//0 [CPX_MIPEMPHASIS_BALANCED] Balance optimality and feasibility
	//1 [CPX_MIPEMPHASIS_FEASIBILITY] Emphasize feasibility over optimality
	//2 [CPX_MIPEMPHASIS_OPTIMALITY] Emphasize optimality over feasibility
	//3 [CPX_MIPEMPHASIS_BESTBOUND] Emphasize moving best bound
	//4 [CPX_MIPEMPHASIS_HIDDENFEAS] Emphasize hidden feasibility
	status = CPXsetintparam (env_COMPACT, CPX_PARAM_MIPEMPHASIS , INPUT_CPX_policy);
#endif

#ifdef set_ub
	// * Set up upper bound *  ( CPX_PARAM_CUTLO the other one)
	status= CPXsetdblparam(env_COMPACT, CPX_PARAM_CUTUP, (double )heur_val+_tollerance);
	if ( status )
	{
		fprintf (stderr, "error for CPX_PARAM_CUTUP\n");
	}
#endif

#ifdef no_cuts
	CPXsetintparam(env_COMPACT,CPX_PARAM_CUTPASS,-1);//no cuts
#endif

#ifdef no_preprocessing
	CPXsetintparam (env_COMPACT,CPX_PARAM_REDUCE, CPX_OFF);
	CPXsetintparam(env_COMPACT,CPX_PARAM_PRELINEAR,CPX_OFF);
	CPXsetintparam(env_COMPACT,CPX_PARAM_PREIND,CPX_OFF);
	CPXsetintparam(env_COMPACT,CPX_PARAM_PROBE,CPX_OFF);
#endif

#ifdef strong_cuts
	CPXsetintparam (env_COMPACT, CPX_PARAM_CLIQUES,2);
	CPXsetintparam (env_COMPACT, CPX_PARAM_COVERS,2);
	CPXsetintparam (env_COMPACT, CPX_PARAM_FLOWCOVERS,2);
	CPXsetintparam (env_COMPACT, CPX_PARAM_IMPLBD,2);
	CPXsetintparam (env_COMPACT, CPX_PARAM_GUBCOVERS,2);
	CPXsetintparam (env_COMPACT, CPX_PARAM_FRACCUTS,2);
	CPXsetintparam (env_COMPACT, CPX_PARAM_FLOWPATHS,2);
	CPXsetintparam (env_COMPACT, CPX_PARAM_MIRCUTS,2);
	CPXsetintparam (env_COMPACT, CPX_PARAM_DISJCUTS,2);
	CPXsetintparam (env_COMPACT, CPX_PARAM_ZEROHALFCUTS,2);
	CPXsetintparam (env_COMPACT, CPX_PARAM_MCFCUTS,2);
#endif

	//CPXsetintparam(env_COMPACT,CPX_PARAM_NODELIM,0);//solo il nodo radice

#ifdef strong_branch
	CPXsetintparam (env_COMPACT, CPX_PARAM_VARSEL,3); //fa lo strong branching
#endif

	//CPXsetintparam(env_COMPACT,CPX_PARAM_PREIND,0);//toglie il preprocessamento // PROPRIO NO ( Q NON SEMIDEFINITA POSITIVA)
	//CPXfreepresolve (env_COMPACT, lp_COMPACT);
	//CPXsetintparam (env_COMPACT,CPX_PARAM_AGGFILL , CPX_OFF); //SOLO NO
	//CPXsetintparam(env_COMPACT,CPX_PARAM_PRELINEAR,CPX_OFF); //SOLO NO
	//CPXsetintparam (env_COMPACT,CPX_PARAM_AGGIND , CPX_OFF);
	//CPXsetintparam (env_COMPACT,CPX_PARAM_REDUCE, CPX_OFF);


	if (optionClique) {
		int ** cliquesT = (int**) malloc(sizeof(int*) * k);
		int * subgraphT = (int*) malloc(sizeof(int) * G->n);
		for (int v = 0; v < G->n; v++) {
			subgraphT[v] = 1;
		}
		for (int i = 0; i < k; i++) {
			cliquesT[i] = (int*) malloc(sizeof(int) * G->n);
			clique_maximum_exact(G, k, q, b, subgraphT, cliquesT[i]);
			/*for(int v=0;v<G->n;v++){
			 cliquesT[i][v]=0;
			 }
			 cliquesT[i][i]=1;
			 */
			for (int v = 0; v < G->n; v++) {
				if (cliquesT[i][v] == 1)
					subgraphT[v] = 0;
			}
		}
		int ** variablesT = (int**) malloc(sizeof(int*) * k);
		for (int i = 0; i < k; i++) {
			variablesT[i] = (int*) malloc(sizeof(int) * G->n);
			for (int v = 0; v < G->n; v++) {
				variablesT[i][v] = 1;
			}
		}
		for (int i = 0; i < k - 1; i++) {
			for (int v = 0; v < G->n; v++) {
				if (cliquesT[i][v] == 1) {
					for (int j = i + 1; j < k; j++) {
						variablesT[j][v] = 0;
					}
				}
			}
		}

		//creo le colonne
		ccnt = _dimension;
		obj = (double*) calloc(ccnt, sizeof(double));
		lb = (double*) calloc(ccnt, sizeof(double));
		ub = (double*) calloc(ccnt, sizeof(double));
		c_type = (char*) calloc(ccnt, sizeof(char));

		char **colname = new char*[ccnt];
		for (int i = 0; i < ccnt; i++) {
			colname[i] = new char[100];
		}

		int dummy = 0;
		for (int v = 0; v < G->n; v++) {
			for (int i = 0; i < k; i++) {
				//    		if(i==0){
				//    			obj[dummy]=0.0;
				//    		}
				//    		else{
				obj[dummy] = 1.0;
				//    		}

				lb[dummy] = 0.0;
				if (variablesT[i][v] == 1) {
					ub[dummy] = 1.0;
				} else {
					ub[dummy] = 0.0;
				}

				c_type[dummy] = 'B';
				sprintf(colname[dummy], "x(%d.%d).%d", v, i, dummy);
				dummy++;
			}
		}
		status = CPXnewcols(env_COMPACT, lp_COMPACT, ccnt, obj, lb, ub, c_type,
				colname);
		if (status != 0) {
			printf("error in CPXnewcols\n");
			exit(-1);
		}

		for (int i = 0; i < ccnt; i++) {
			delete[] colname[i];
		}
		delete[] colname;

		free(obj);
		free(lb);
		free(ub);
		free(c_type);

		//		cin.get();
		//		cout<<"begin free"<<endl;
		free(subgraphT);
		for (int i = 0; i < k; i++) {
			free(variablesT[i]);
		}
		free(variablesT);
		for (int i = 0; i < k; i++) {
			free(cliquesT[i]);
		}
		free(cliquesT);
		//		cout<<"free"<<endl;
		//		cin.get();

	} else {
		//creo le colonne
		ccnt = _dimension;
		obj = (double*) calloc(ccnt, sizeof(double));
		lb = (double*) calloc(ccnt, sizeof(double));
		ub = (double*) calloc(ccnt, sizeof(double));
		c_type = (char*) calloc(ccnt, sizeof(char));

		char **colname = new char*[ccnt];
		for (i = 0; i < ccnt; i++) {
			colname[i] = new char[100];
		}

		//creating variables
		int dummy = 0;
		for (int v = 0; v < G->n; v++) {
			for (int i = 0; i < k; i++) {
				//    		if(i==0){
				//    			obj[dummy]=0.0;
				//    		}
				//    		else{
				obj[dummy] = 1.0;
				//    		}

				lb[dummy] = 0.0;
				ub[dummy] = 1.0;
				c_type[dummy] = 'B';
				sprintf(colname[dummy], "x(%d.%d).%d", v, i, dummy);
				dummy++;
			}
		}
		status = CPXnewcols(env_COMPACT, lp_COMPACT, ccnt, obj, lb, ub, c_type,
				colname);
		if (status != 0) {
			printf("error in CPXnewcols\n");
			exit(-1);
		}

		for (i = 0; i < ccnt; i++) {
			delete[] colname[i];
		}
		delete[] colname;

		free(obj);
		free(lb);
		free(ub);
		free(c_type);

	}
	//routine cplex per creare le colonne

	//each vertex must go in at most one partition
	for (int v = 0; v < G->n; v++) {

		rcnt = 1;
		nzcnt = k;

		rhs = (double*) calloc(rcnt, sizeof(double));
		sense = (char*) calloc(rcnt, sizeof(double));

		rhs[0] = 1.0;
		sense[0] = 'L';

		rmatbeg = (int*) calloc(rcnt, sizeof(int));
		rmatind = (int*) calloc(nzcnt, sizeof(int));
		rmatval = (double*) calloc(nzcnt, sizeof(double));

		for (int i = 0; i < k; i++) {
			rmatval[i] = 1.0;
			rmatind[i] = position_x_v_i(v, i);
		}

		//for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
		//for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

		rmatbeg[0] = 0;

		status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt, nzcnt, rhs,
				sense, rmatbeg, rmatind, rmatval, NULL, NULL);
		if (status != 0) {
			printf("error in CPXaddrows\n");
			exit(-1);
		}

		free(rmatbeg);
		free(rmatval);
		free(rmatind);
		free(rhs);
		free(sense);
	}

	if (option4) {
		int *no_done = (int*) malloc(sizeof(int) * G->n);
		for (int i = 0; i < G->n; i++) {
			no_done[i] = 1;
		}
		int nb_done = 0;
		int find = 1;
		for (int v = 0; v < G->n; v++) {

			//		while (nb_done<G->n-1) {
			//			cin.get();
			//			cout<<" "<<nb_done<<" - "<<G->n<<endl;
			int *clique = (int*) malloc(sizeof(int) * G->n);
			//clique_maximum_exact(G, k, q, b, no_done, clique);
			for (int u = 0; u < G->n; u++) {
				clique[u] = 0;
			}
			clique[v] = 1;
			maximize_clique(G, clique);
			int size_clique = 0;
			for (int u = 0; u < G->n; u++) {
				if (clique[u] == 1)
					size_clique++;
			}
			if (size_clique > 2) {

				/*cout<<"clique "<<endl;

				 for(int deb=0;deb<G->n;deb++)
				 cout<<" "<<clique[deb];
				 cout<<endl;
				 cin.get();

				 cout<<"nogood "<<endl;

				 for(int deb=0;deb<G->n;deb++)
				 cout<<" "<<no_done[deb];
				 cout<<endl;
				 */

				rcnt = 1;
				nzcnt = G->n * k;

				rhs = (double*) calloc(rcnt, sizeof(double));
				sense = (char*) calloc(rcnt, sizeof(double));

				rhs[0] = k - 1;
				sense[0] = 'G';

				rmatbeg = (int*) calloc(rcnt, sizeof(int));
				rmatind = (int*) calloc(nzcnt, sizeof(int));
				rmatval = (double*) calloc(nzcnt, sizeof(double));

				int pos_current = 0;
				for (int vp = 0; vp < G->n; vp++) {
					if (clique[vp] == 0) {
						no_done[vp] = 0;
						nb_done++;
						for (int ip = 0; ip < k; ip++) {
							//attention i
							rmatval[pos_current] = 1.0;
							rmatind[pos_current] = position_x_v_i(vp, ip);
							pos_current++;
						}
					} else {
						for (int ip = 0; ip < k; ip++) {
							rmatval[pos_current] = 0;
							rmatind[pos_current] = position_x_v_i(vp, ip);
							pos_current++;
						}
					}
				}

				rmatbeg[0] = 0;

				status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt, nzcnt,
						rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status != 0) {
					printf("error in CPXaddrows\n");
					exit(-1);
				}

				free(rmatbeg);
				free(rmatval);
				free(rmatind);
				free(rhs);
				free(sense);
			}
			free(clique);
		}
		free(no_done);
	}

	if (option3) {
		for (int v = 0; v < G->n; v++) {
			for (int u = 0; u < G->n; u++) {
				if (u == v) {
					continue;
				}
				int comonNeighbour = 0;
				for (int w = 0; w < G->n; w++) {
					if ((G->AMatrix[w][v] == 1 || G->AMatrix[v][w] == 1)
							&& (G->AMatrix[w][u] == 1 || G->AMatrix[u][w] == 1)) {
						comonNeighbour = 1;
					}
				}
				if (!comonNeighbour) {
					int alpha_u = stable_maximum_exact(G, k, q, b,
							neighboors(G, k, q, b, u));
					int alpha_v = stable_maximum_exact(G, k, q, b,
							neighboors(G, k, q, b, v));
					if (alpha_u > 1 && alpha_v > 1) {
						int* subgraph = (int*) malloc(G->n * sizeof(int));
						for (int vp = 0; vp < G->n; vp++) {
							subgraph[vp] = 0;
						}
						for (int e = 0; e < G->m; e++) {
							if (G->H[e] == v || G->H[e] == u) {
								subgraph[G->T[e]] = 1;
							}
							if (G->T[e] == v || G->T[e] == u) {
								subgraph[G->H[e]] = 1;
							}
						}
						subgraph[v] = 1;
						subgraph[u] = 1;
						///
						for (int i1 = 0; i1 < k; i1++) {
							for (int i2 = 0; i2 < k; i2++) {
								if (i1 == i2) {
									continue;
								}
								rcnt = 1;
								nzcnt = G->n * k;

								rhs = (double*) calloc(rcnt, sizeof(double));
								sense = (char*) calloc(rcnt, sizeof(double));

								rhs[0] = k - alpha_v - 1 - alpha_u - 1;
								sense[0] = 'G';

								rmatbeg = (int*) calloc(rcnt, sizeof(int));
								rmatind = (int*) calloc(nzcnt, sizeof(int));
								rmatval = (double*) calloc(nzcnt,
										sizeof(double));

								int pos_current = 0;
								for (int vp = 0; vp < G->n; vp++) {
									for (int ip = 0; ip < k; ip++) {
										//attention i
										if (subgraph[vp] == 0 && ip != i1 && ip
												!= i2) {
											rmatval[pos_current] = 1.0;
											rmatind[pos_current]
													= position_x_v_i(vp, ip);
											pos_current++;
										} else {
											if (vp == v && ip == i1) {
												rmatval[pos_current] = -alpha_v;
												rmatind[pos_current]
														= position_x_v_i(vp, ip);
												pos_current++;
											} else if (vp == u && ip == i2) {
												rmatval[pos_current] = -alpha_u;
												rmatind[pos_current]
														= position_x_v_i(vp, ip);
												pos_current++;
											} else if (subgraph[vp] == 1 && ip
													== i) {
												rmatval[pos_current] = 0;
												rmatind[pos_current]
														= position_x_v_i(vp, ip);
												pos_current++;
											} else {
												rmatval[pos_current] = 0.0;
												rmatind[pos_current]
														= position_x_v_i(vp, ip);
												pos_current++;

											}
										}
									}
								}

								rmatbeg[0] = 0;

								status = CPXaddrows(env_COMPACT, lp_COMPACT, 0,
										rcnt, nzcnt, rhs, sense, rmatbeg,
										rmatind, rmatval, NULL, NULL);
								if (status != 0) {
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

						///
						free(subgraph);
					}
				}
			}
		}

	}

	//no neighbours inequalities
	//int option2 = 1;
	if (option2) {
		for (int v = 0; v < G->n; v++) {
			int alpha_v = stable_maximum_exact(G, k, q, b,
					neighboors(G, k, q, b, v));//G->n;//replace by the stable set size in $v \cup N(v)$
			int N_v = 0;
			int* subgraph = (int*) malloc(G->n * sizeof(int));
			for (int vp = 0; vp < G->n; vp++)
				subgraph[vp] = 0;
			for (int e = 0; e < G->m; e++) {
				if (G->H[e] == v) {
					subgraph[G->T[e]] = 1;
					N_v++;
				}
				if (G->T[e] == v) {
					subgraph[G->H[e]] = 1;
					N_v++;
				}
			}
			subgraph[v] = 1;

			//			alpha_v=N_v;
			//ct10
			if (alpha_v > 1) {
				for (i = 0; i < k; i++) {

					rcnt = 1;
					nzcnt = G->n * k;

					rhs = (double*) calloc(rcnt, sizeof(double));
					sense = (char*) calloc(rcnt, sizeof(double));

					rhs[0] = k - alpha_v - 1;
					sense[0] = 'G';

					rmatbeg = (int*) calloc(rcnt, sizeof(int));
					rmatind = (int*) calloc(nzcnt, sizeof(int));
					rmatval = (double*) calloc(nzcnt, sizeof(double));

					int pos_current = 0;
					for (int vp = 0; vp < G->n; vp++) {
						for (int ip = 0; ip < k; ip++) {
							//attention i
							if (subgraph[vp] == 0 && ip != i) {
								rmatval[pos_current] = 1.0;
								rmatind[pos_current] = position_x_v_i(vp, ip);
								pos_current++;
							} else {
								if (vp == v && ip == i) {
									rmatval[pos_current] = -alpha_v;
									rmatind[pos_current] = position_x_v_i(vp,
											ip);
									pos_current++;
								} else if (subgraph[vp] == 1 && ip == i) {
									rmatval[pos_current] = 0;
									rmatind[pos_current] = position_x_v_i(vp,
											ip);
									pos_current++;
								} else {
									rmatval[pos_current] = 0.0;
									rmatind[pos_current] = position_x_v_i(vp,
											ip);
									pos_current++;

								}
							}
						}
					}

					rmatbeg[0] = 0;

					status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt,
							nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL,
							NULL);
					if (status != 0) {
						printf("error in CPXaddrows\n");
						exit(-1);
					}

					free(rmatbeg);
					free(rmatval);
					free(rmatind);
					free(rhs);
					free(sense);

				}

				//ct9

				rcnt = 1;
				nzcnt = G->n * k;

				rhs = (double*) calloc(rcnt, sizeof(double));
				sense = (char*) calloc(rcnt, sizeof(double));

				rhs[0] = k - alpha_v;
				sense[0] = 'G';

				rmatbeg = (int*) calloc(rcnt, sizeof(int));
				rmatind = (int*) calloc(nzcnt, sizeof(int));
				rmatval = (double*) calloc(nzcnt, sizeof(double));

				int pos_current = 0;
				for (int vp = 0; vp < G->n; vp++) {
					for (int ip = 0; ip < k; ip++) {
						if (subgraph[vp] == 0 && ip != i) {
							rmatval[pos_current] = 1.0;
							rmatind[pos_current] = position_x_v_i(vp, ip);
							pos_current++;
						} else {
							if (vp == v) {
								rmatval[pos_current] = 1 - alpha_v;
								rmatind[pos_current] = position_x_v_i(vp, ip);
								pos_current++;
							} else {
								rmatval[pos_current] = 0.0;
								rmatind[pos_current] = position_x_v_i(vp, ip);
								pos_current++;

							}
						}
					}
				}

				rmatbeg[0] = 0;

				status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt, nzcnt,
						rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status != 0) {
					printf("error in CPXaddrows\n");
					exit(-1);
				}
				free(subgraph);
				free(rmatbeg);
				free(rmatval);
				free(rmatind);
				free(rhs);
				free(sense);
			}
		}
	}

	//if there is an edge the vertices must go in different partitions (or in the separator)
	//int option1=1;
	if (option1) {
		for (int i = 0; i < k; i++) {
			//			for (int j = i + 1; j < k; j++) {

			//    		if(i==j){continue;}

			for (int e = 0; e < G->m; e++) {

				//arc v -> u

				rcnt = 1;
				nzcnt = k;

				rhs = (double*) calloc(rcnt, sizeof(double));
				sense = (char*) calloc(rcnt, sizeof(double));

				rhs[0] = 1.0;
				sense[0] = 'L';

				rmatbeg = (int*) calloc(rcnt, sizeof(int));
				rmatind = (int*) calloc(nzcnt, sizeof(int));
				rmatval = (double*) calloc(nzcnt, sizeof(double));

				rmatval[0] = 1.0;
				rmatind[0] = position_x_v_i(G->H[e], i);
				int pos_current = 1;
				for (int j = 0; j < k; j++) {
					if (j != i) {
						rmatval[pos_current] = 1.0;
						rmatind[pos_current] = position_x_v_i(G->T[e], j);
						pos_current++;
					}
				}
				//for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
				//for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

				rmatbeg[0] = 0;

				status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt, nzcnt,
						rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status != 0) {
					printf("error in CPXaddrows\n");
					exit(-1);
				}

				free(rmatbeg);
				free(rmatval);
				free(rmatind);
				free(rhs);
				free(sense);

				//arc u -> v
				rcnt = 1;
				nzcnt = k;

				rhs = (double*) calloc(rcnt, sizeof(double));
				sense = (char*) calloc(rcnt, sizeof(double));

				rhs[0] = 1.0;
				sense[0] = 'L';

				rmatbeg = (int*) calloc(rcnt, sizeof(int));
				rmatind = (int*) calloc(nzcnt, sizeof(int));
				rmatval = (double*) calloc(nzcnt, sizeof(double));

				rmatval[0] = 1.0;
				rmatind[0] = position_x_v_i(G->T[e], i);

				pos_current = 1;
				for (int j = 0; j < k; j++) {
					if (j != i) {
						rmatval[pos_current] = 1.0;
						rmatind[pos_current] = position_x_v_i(G->H[e], j);
						pos_current++;
					}
				}

				//for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
				//for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

				rmatbeg[0] = 0;

				status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt, nzcnt,
						rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status != 0) {
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

	} else {
		/*
		 for (int i = 0; i < k ; i++) {
		 for (int j = 0; j < k; j++) {

		 if(i==j){continue;}

		 for (int e = 0; e < G->m; e++) {

		 //arc v -> u

		 rcnt = 1;
		 nzcnt = 2;

		 rhs = (double*) calloc(rcnt, sizeof(double));
		 sense = (char*) calloc(rcnt, sizeof(double));

		 rhs[0] = 1.0;
		 sense[0] = 'L';

		 rmatbeg = (int*) calloc(rcnt, sizeof(int));
		 rmatind = (int*) calloc(nzcnt, sizeof(int));
		 rmatval = (double*) calloc(nzcnt, sizeof(double));

		 rmatval[0] = 1.0;
		 rmatind[0] = position_x_v_i(G->H[e], i);

		 rmatval[1] = 1.0;
		 rmatind[1] = position_x_v_i(G->T[e], j);

		 //for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
		 //for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

		 rmatbeg[0] = 0;

		 status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt,
		 nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL,
		 NULL);
		 if (status != 0) {
		 printf("error in CPXaddrows\n");
		 exit(-1);
		 }

		 free(rmatbeg);
		 free(rmatval);
		 free(rmatind);
		 free(rhs);
		 free(sense);

		 //arc u -> v
		 rcnt = 1;
		 nzcnt = 2;

		 rhs = (double*) calloc(rcnt, sizeof(double));
		 sense = (char*) calloc(rcnt, sizeof(double));

		 rhs[0] = 1.0;
		 sense[0] = 'L';

		 rmatbeg = (int*) calloc(rcnt, sizeof(int));
		 rmatind = (int*) calloc(nzcnt, sizeof(int));
		 rmatval = (double*) calloc(nzcnt, sizeof(double));

		 rmatval[0] = 1.0;
		 rmatind[0] = position_x_v_i(G->T[e], i);

		 rmatval[1] = 1.0;
		 rmatind[1] = position_x_v_i(G->H[e], j);

		 //for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
		 //for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

		 rmatbeg[0] = 0;

		 status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt,
		 nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, NULL,
		 NULL);
		 if (status != 0) {
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
		 }*/
		for (int i = 0; i < k; i++) {
			for (int j = 0; j < k; j++) {

				if (i != j) {

					for (int e = 0; e < G->m; e++) {

						//arc v -> u

						rcnt = 1;
						nzcnt = 2;

						rhs = (double*) calloc(rcnt, sizeof(double));
						sense = (char*) calloc(rcnt, sizeof(double));

						rhs[0] = 1.0;
						sense[0] = 'L';

						rmatbeg = (int*) calloc(rcnt, sizeof(int));
						rmatind = (int*) calloc(nzcnt, sizeof(int));
						rmatval = (double*) calloc(nzcnt, sizeof(double));

						rmatval[0] = 1.0;
						rmatind[0] = position_x_v_i(G->H[e], i);

						rmatval[1] = 1.0;
						rmatind[1] = position_x_v_i(G->T[e], j);

						//for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
						//for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

						rmatbeg[0] = 0;

						status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt,
								nzcnt, rhs, sense, rmatbeg, rmatind, rmatval,
								NULL, NULL);
						if (status != 0) {
							printf("error in CPXaddrows\n");
							exit(-1);
						}

						free(rmatbeg);
						free(rmatval);
						free(rmatind);
						free(rhs);
						free(sense);

						//arc u -> v
						rcnt = 1;
						nzcnt = 2;

						rhs = (double*) calloc(rcnt, sizeof(double));
						sense = (char*) calloc(rcnt, sizeof(double));

						rhs[0] = 1.0;
						sense[0] = 'L';

						rmatbeg = (int*) calloc(rcnt, sizeof(int));
						rmatind = (int*) calloc(nzcnt, sizeof(int));
						rmatval = (double*) calloc(nzcnt, sizeof(double));

						rmatval[0] = 1.0;
						rmatind[0] = position_x_v_i(G->T[e], i);

						rmatval[1] = 1.0;
						rmatind[1] = position_x_v_i(G->H[e], j);

						//for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
						//for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

						rmatbeg[0] = 0;

						status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt,
								nzcnt, rhs, sense, rmatbeg, rmatind, rmatval,
								NULL, NULL);
						if (status != 0) {
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
		}
	}
	// no empty partitions
	for (int i = 0; i < k; i++) {

		rcnt = 1;
		nzcnt = G->n;

		rhs = (double*) calloc(rcnt, sizeof(double));
		sense = (char*) calloc(rcnt, sizeof(double));

		rhs[0] = 1.0;
		sense[0] = 'G';

		rmatbeg = (int*) calloc(rcnt, sizeof(int));
		rmatind = (int*) calloc(nzcnt, sizeof(int));
		rmatval = (double*) calloc(nzcnt, sizeof(double));

		for (int v = 0; v < G->n; v++) {
			rmatval[v] = 1.0;
			rmatind[v] = position_x_v_i(v, i);
		}

		//for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
		//for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

		rmatbeg[0] = 0;

		status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt, nzcnt, rhs,
				sense, rmatbeg, rmatind, rmatval, NULL, NULL);
		if (status != 0) {
			printf("error in CPXaddrows\n");
			exit(-1);
		}

		free(rmatbeg);
		free(rmatval);
		free(rmatind);
		free(rhs);
		free(sense);
	}

	//balancing the number of nodes between the different partitions
	if (q != -1) {
		for (int i = 0; i < k - 1; i++) {
			for (int j = i + 1; j < k; j++) {

				// couple i -> j

				rcnt = 1;
				nzcnt = 2 * G->n;

				rhs = (double*) calloc(rcnt, sizeof(double));
				sense = (char*) calloc(rcnt, sizeof(double));

				rhs[0] = q;
				sense[0] = 'L';

				rmatbeg = (int*) calloc(rcnt, sizeof(int));
				rmatind = (int*) calloc(nzcnt, sizeof(int));
				rmatval = (double*) calloc(nzcnt, sizeof(double));

				int counter_local = 0;
				for (int v = 0; v < G->n; v++) {
					rmatval[counter_local] = 1.0;
					rmatind[counter_local] = position_x_v_i(v, i);
					counter_local++;
				}
				for (int v = 0; v < G->n; v++) {
					rmatval[counter_local] = -1.0;
					rmatind[counter_local] = position_x_v_i(v, j);
					counter_local++;
				}

				//for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
				//for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

				rmatbeg[0] = 0;

				status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt, nzcnt,
						rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status != 0) {
					printf("error in CPXaddrows\n");
					exit(-1);
				}

				free(rmatbeg);
				free(rmatval);
				free(rmatind);
				free(rhs);
				free(sense);

				// couple j -> i


				rcnt = 1;
				nzcnt = 2 * G->n;

				rhs = (double*) calloc(rcnt, sizeof(double));
				sense = (char*) calloc(rcnt, sizeof(double));

				rhs[0] = q;
				sense[0] = 'L';

				rmatbeg = (int*) calloc(rcnt, sizeof(int));
				rmatind = (int*) calloc(nzcnt, sizeof(int));
				rmatval = (double*) calloc(nzcnt, sizeof(double));

				counter_local = 0;
				for (int v = 0; v < G->n; v++) {
					rmatval[counter_local] = 1.0;
					rmatind[counter_local] = position_x_v_i(v, j);
					counter_local++;
				}
				for (int v = 0; v < G->n; v++) {
					rmatval[counter_local] = -1.0;
					rmatind[counter_local] = position_x_v_i(v, i);
					counter_local++;
				}

				//for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
				//for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

				rmatbeg[0] = 0;

				status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt, nzcnt,
						rhs, sense, rmatbeg, rmatind, rmatval, NULL, NULL);
				if (status != 0) {
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

	// max vertices per partition
	if (b != -1) {
		for (int i = 0; i < k; i++) {

			rcnt = 1;
			nzcnt = G->n;

			rhs = (double*) calloc(rcnt, sizeof(double));
			sense = (char*) calloc(rcnt, sizeof(double));

			rhs[0] = b;
			sense[0] = 'L';

			rmatbeg = (int*) calloc(rcnt, sizeof(int));
			rmatind = (int*) calloc(nzcnt, sizeof(int));
			rmatval = (double*) calloc(nzcnt, sizeof(double));

			for (int v = 0; v < G->n; v++) {
				rmatval[v] = 1.0;
				rmatind[v] = position_x_v_i(v, i);
			}

			//for(int jj=0; jj<n_position_total; jj++){cout << rmatval[jj] << "\t";}cout << endl;
			//for(int jj=0; jj<n_position_total; jj++){cout << rmatind[jj] << "\t";}cout << endl;

			rmatbeg[0] = 0;

			status = CPXaddrows(env_COMPACT, lp_COMPACT, 0, rcnt, nzcnt, rhs,
					sense, rmatbeg, rmatind, rmatval, NULL, NULL);
			if (status != 0) {
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

#ifdef write_problem

	cout << "\n\nwrite compact model\n";
	char *_name=new char[100];
	sprintf(_name,"model_base.lp");

	status=CPXwriteprob(env_COMPACT,lp_COMPACT,_name,NULL);
	if(status!=0) {printf("error in CPXwriteprob\n"); while(1);exit(-1);}

	delete [] _name;

	exit(-1);
#endif

}

/***************************************************************************/
double compact_model_solve(char *istname, graphFF G, int k, int q, int b)
/***************************************************************************/
{

	//////////////////////////////////////////////////////////////////////////////////////////////////
	int _dimension = G->n * k;
	double* x_exact;

	bool _TimeLimit = false;

	x_exact = (double*) calloc(_dimension, sizeof(double));

	for (int i = 0; i < _dimension; i++) {
		x_exact[i] = 0;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////
	status = CPXwriteprob(env_COMPACT, lp_COMPACT, "lp.lp", NULL);

	if (option == 2) {


		clock_t time_start = clock();



		status = CPXchgprobtype(env_COMPACT, lp_COMPACT, CPXPROB_LP);

		status = CPXwriteprob(env_COMPACT, lp_COMPACT, "lp.lp", NULL);
		if (status != 0) {
			printf("error in CPXwriteprob\n");
			exit(-1);
		}

		status = CPXlpopt(env_COMPACT, lp_COMPACT);
		if (status != 0) {
			cout << "err_FILEor in CPXlpopt slave solve\n";
			exit(-1);
		}

		clock_t time_end = clock();
		double solution_time = (double) (time_end - time_start)
								/ (double) CLOCKS_PER_SEC;


		int lp_stat = CPXgetstat(env_COMPACT, lp_COMPACT);

		double objval_lp=-1;

		status = CPXgetobjval(env_COMPACT, lp_COMPACT, &objval_lp);
		if (status != 0) {
			printf("error in CPXgetmipobjval\n");
		}

		cout << "\n\nCR ->\t " << objval_lp  << "\t time \t " << solution_time << "\t stat \t " << lp_stat << "\t" << istname << "\t" << k << endl;

#ifdef print_sol_details_lp
		//getting the solution
		status= CPXgetmipx(env_COMPACT, lp_COMPACT, x_exact, 0,
				_dimension - 1);
		if (status != 0) {
			cout << endl << "\nNO SOLUTION!!" << endl;
			exit(-1);
			//err_FILE <<"error in CPXgetmipx model_solve\n";
		}


		cout << "\n\nCURRENT VARIABLES:\n";
		for (int v = 0; v < G->n; v++) {
			for (int i = 0; i < k; i++) {
				printf("%.3f\t", x_exact[position_x_v_i(v, i)]);
			}
			cout << endl;
		}
#endif

		exit(-1);

	}

	//only the root node
	if (option == 3) {
		CPXsetintparam(env_COMPACT, CPX_PARAM_NODELIM, 0);
	}

	//////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////

#ifdef BC
	CPXsetintparam(env_COMPACT,CPX_PARAM_PRELINEAR,CPX_OFF);
	CPXsetintparam(env_COMPACT,CPX_PARAM_MIPCBREDLP,CPX_OFF);
	vertex_number=G->n;
	///////////////////////////////////////////////
	CPXsetusercutcallbackfunc(env_COMPACT,CPXsetusercutcallbackfunc ,NULL);
	///////////////////////////////////////////////
	//for the B&C with Lazy
	parametresCplex paramC;
	paramC.G=G;
	paramC.k=k;
	CPXsetlazyconstraintcallbackfunc(env_COMPACT, CPXsetlazycutcallbackfunc, &paramC);

#endif



	clock_t time_start = clock();

	//solving the problem
	status = CPXmipopt(env_COMPACT, lp_COMPACT);
	if (status != 0) {
		printf("error in CPXmipopt\n");
		exit(-1);
	}

	clock_t time_end = clock();
	double solution_time = (double) (time_end - time_start)
					/ (double) CLOCKS_PER_SEC;

#ifdef BC
	///////////////////////////////////////////////
	CPXsetusercutcallbackfunc(env_COMPACT,NULL,NULL);
	///////////////////////////////////////////////
#endif

	int status_mipopt = CPXgetstat(env_COMPACT, lp_COMPACT);
	if (status_mipopt == 107 || status_mipopt == 108) {
		_TimeLimit = true;
		cout << endl << "\nTIME LIMIT!!" << endl;
	}

	bool sol_ok = true;
	//getting the solution
	status = CPXgetmipx(env_COMPACT, lp_COMPACT, x_exact, 0, _dimension - 1);
	if (status != 0) {
		sol_ok = false;
		cout << endl << "\nNO SOLUTION!!" << endl;
		//err_FILE <<"error in CPXgetmipx model_solve\n";
	}

	double optValue;

	status = CPXgetmipobjval(env_COMPACT, lp_COMPACT, &optValue);
	cout << "hello " << endl;

	cout << "\n\nwrite compact model\n";

	status = CPXwriteprob(env_COMPACT, lp_COMPACT, "test.lp", NULL);
	//	if(status!=0) {printf("error in CPXwriteprob\n");	while(1);exit(-1);}

	//	exit(-1);


	if (status != 0) {
		cout << endl << "\nNO SOLUTION!!" << endl;
		sol_ok = false;
		optValue = 0;
	}

	double best_ub;

	status = CPXgetbestobjval(env_COMPACT, lp_COMPACT, &best_ub);
	if (status != 0) {
		cout << endl << "\nNO SOLUTION!!" << endl;
		best_ub = 0;
	}

	int node_number = CPXgetnodecnt(env_COMPACT, lp_COMPACT);
	cout << "---->" << best_ub << endl;
	cout << fixed << "\nOBJ VAL\t\t\t" << optValue << endl;
	cout << fixed << "NODES\t\t\t" << node_number << endl;

	int *vertex_partition = new int[k];
	for (int z = 0; z < k; z++) {
		vertex_partition[z] = 0;
	}

	if (sol_ok) {
		for (int v = 0; v < G->n; v++) {
			for (int i = 0; i < k; i++) {
				if (x_exact[position_x_v_i(v, i)] > 0.5) {
					vertex_partition[i] = vertex_partition[i] + 1;
				}
			}
		}

		cout << "Vertices per partition:\n";
		for (int i = 0; i < k; i++) {
			cout << vertex_partition[i] << "\t";
		}
		cout << endl;

	}

#ifdef print_sol_details
	cout << "\n\nCURRENT VARIABLES:\n";
	for (int v = 0; v < G->n; v++) {
		for (int i = 0; i < k; i++) {
			cout << (int) (x_exact[position_x_v_i(v, i)] + 0.5) << "\t";
		}
		cout << endl;
	}
#endif

	int nodes_compact = CPXgetnodecnt(env_COMPACT, lp_COMPACT);
	int lp_stat = CPXgetstat(env_COMPACT, lp_COMPACT);

	cout << "lp_stat " << lp_stat << endl;

	int numcols_C = CPXgetnumcols(env_COMPACT, lp_COMPACT);
	int numrows_C = CPXgetnumrows(env_COMPACT, lp_COMPACT);
	/*
	ofstream info_SUMMARY("info_compact.txt", ios::app);
	info_SUMMARY << fixed << istname << "\t" << G->n << "\t" << G->m << "\t"
			<< algorithm << "\t" << option << "\t" << k << "\t" << q << "\t"
			<< b << "\t" << optValue << "\t" << best_ub << "\t"
			<< nodes_compact << "\t" << numcols_C << "\t" << numrows_C << "\t"
			<< solution_time << "\t" << lp_stat << "\t\t";

	for (int i = 0; i < k; i++) {
		info_SUMMARY << vertex_partition[i] << "\t";
	}
	info_SUMMARY << endl;

	info_SUMMARY.close();
	*/
	if(output.is_open())
	{
			output.precision(10);
			output 	<< (G->n - optValue) << "\t" << (G->n-best_ub) << "\t"  // << rootObjective << "\t" 
				<< solution_time << "\t" << nodes_compact << "\t" 
				<< numcols_C << "\t" << numrows_C << "\t" << lp_stat;
	}
	/////////////////////////
	free(x_exact);
	////////////////////////


	return optValue;

}

