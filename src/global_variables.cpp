

#include "global_variables.h"

//////////////////////////////////////////////////
graphFF G;

int counter_clique;
int **clique;
int *sizeClique;

int disj_counter_clique;
int **disj_clique;
int *disj_sizeClique;

int *parent;
double *wweight;
int *iindex;
int *vindex;
double *vweight;
int *bondersCoeff;

int *statusDFS;
bool *visited;
bool *removed;
int *goodVertex;
int *cliqueCut;
int **rGraph;

//Used for Connecticity cut
int* depth;
int* low;
int* localDim;

//made global to be used inside the user cut, if you don't need, remake local!
int **symmAdjMatrix;

int userCuts=0;
int lazyCuts=0;

int bppCuts=0;
int lazyCall=0;
int compCuts=0;
int kConnCuts=0;

IloExprArray capConstrArray;
//////////////////////////////////////////////////


/////////////////////////
int algorithm;
int option;
int k;
int q;
int b;
double cap;
double timeLimit;
double tolerance;
double treeLimit;
int userFreq;
int seed;
char* istname;
/////////////////////////

/////////////////////////////////////CPLEX/////////////////////////////////////

CPXENVptr env_COMPACT;
CPXLPptr lp_COMPACT;

CPXENVptr env_kp; 
CPXLPptr lp_kp; 

int status,ccnt,rcnt,nzcnt; 
int* rmatbeg, *rmatind,*cmatbeg, *cmatind; 
double* rmatval,*cmatval,*x,*pi,*obj, *lb, *ub,*rhs;
char *c_type,* sense;

IloNumArray userSolx;
IloNumArray lazySolx;

IloIntArray goodCliqueU;
IloIntArray degreeTree;

/////////////////////////////////////////
int edge_number;
int vertex_number;
int *fr;
int *to;

ofstream output;
//////////////////////////////////////////

//NB: they are necessary for the Bin Packing Problem
//This vector collects the dimensions of each connected component in the graph
vector<int> weights;
//This 2D vector collects the list of vertices for each connected component in the graph
vector<vector<int> > compVertices;
