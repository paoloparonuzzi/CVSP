#ifndef VARIABLE_local_HEADER
#define VARIABLE_local_HEADER


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <string>
#include <iomanip> // for precision
#include "ilcplex/ilocplex.h"
#include "Graph_v4.h"

//#include </Program Files\IBM\ILOG\CPLEX_Studio126\cplex\include\ilcplex/cplex.h>
//#include </opt/ibm/ILOG/CPLEX_Studio126/cplex/include/ilcplex/cplex.h>
#include </opt/ibm/ILOG/CPLEX_Studio1271/cplex/include/ilcplex/cplex.h>

using namespace std;

//////////////////////////////////////////////////
extern graphFF G;

extern int counter_clique;
extern int **clique;
extern int *sizeClique;

extern int disj_counter_clique;
extern int **disj_clique;
extern int *disj_sizeClique;

extern int *parent;
extern double *wweight;
extern int *iindex;
extern int *vindex;
extern double *vweight;
extern int *bondersCoeff;

extern int *statusDFS;
//made global to be used inside the user cut, if you don't need, remake local!
extern int **symmAdjMatrix;

extern bool *visited;
extern bool *removed;
extern int *goodVertex;
extern int *cliqueCut;
extern int **rGraph;

//Used for Connecticity cut
extern int* depth;
extern int* low;
extern int* localDim;

extern int userCuts;
extern int lazyCuts;

extern int bppCuts;
extern int lazyCall;
extern int compCuts;
extern int kConnCuts;

extern IloExprArray capConstrArray;
//////////////////////////////////////////////////


/////////////////////////
extern int algorithm;
extern int option;
extern int k;
extern int q;
extern int b;
extern double cap;
extern double timeLimit;
extern double tolerance;
extern double treeLimit;
extern int userFreq;
extern int seed;
extern char* istname;
/////////////////////////

/////////////////////////////////////CPLEX/////////////////////////////////////

extern CPXENVptr env_COMPACT;
extern CPXLPptr lp_COMPACT;

extern CPXENVptr env_kp; 
extern CPXLPptr lp_kp; 

extern int status,ccnt,rcnt,nzcnt;
extern int* rmatbeg, *rmatind,*cmatbeg, *cmatind;
extern double* rmatval,*cmatval,*x,*pi,*obj, *lb, *ub,*rhs;
extern char *c_type,* sense;

extern IloNumArray userSolx;
extern IloNumArray lazySolx;

extern IloIntArray goodCliqueU;
extern IloIntArray degreeTree;

/////////////////////////////////////////
extern int edge_number;
extern int vertex_number;
extern int *to;

extern ofstream output;
//////////////////////////////////////////

//NB: they are necessary for the Bin Packing Problem
//This vector collects the dimensions of each connected component in the graph
extern vector<int> weights;
//This 2D vector collects the list of vertices for each connected component in the graph
extern vector<vector<int> > compVertices;

#endif
