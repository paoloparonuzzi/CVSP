#ifndef FUNCION_local_HEADER
#define FUNCION_local_HEADER


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
#include <iomanip> // for precision
#include <queue>
#include <list>

#include "global_variables.h"
#include "Graph_v4.h"
#include <random>


using namespace std;

/***************************************************************************/
int    goldberg(int, int, int, int, int*, int*, double*, double*, double*, int*, short int *);
/***************************************************************************/

/***************************************************************************/
void SORT_NON_INCR_INT(int *item,int *score,int n);
/***************************************************************************/

/***************************************************************************/
void kp_load_cplex(int n_item, double C, double *weights);
/***************************************************************************/

/***************************************************************************/
double kp_solve_cplex(int n_item,double *profits, double *solution);
/***************************************************************************/

// Clique algorithm for separation of clique path constraint
/***************************************************************************/
double findGoodClique(int uIndex, int vIndex, graphFF G, IloNumArray zStar, IloIntArray goodClique);
/***************************************************************************/

/***************************************************************************/
void lazyArrayDFS(IloExpr *expr, IloIntVarArray x, int v, int nComp, int *dim);
/***************************************************************************/

/***************************************************************************/
void lazyArrayBFS(IloExpr *expr, IloIntVarArray x, int v, int nComp, int *dim);
/***************************************************************************/

// Clique algorithm for separation of clique in the capacity version
/***************************************************************************/
double findCapClique(int uIndex, graphFF G, IloNumArray xStar);
/***************************************************************************/

/***************************************************************************/
bool bfs(int s, int t);
/***************************************************************************/

/***************************************************************************/
int fordFulkerson(int s, int t, int dim, int nComp, int minDegree);
/***************************************************************************/

/***************************************************************************/
int kVertexConnect(int dim, int nComp, int minDegree);
/***************************************************************************/

/***************************************************************************/
void simpleDFS(int v, int nComp, int *dim);
/***************************************************************************/

/***************************************************************************/
void compTreeDFS(int v, int nComp, int *dim);
/***************************************************************************/

/***************************************************************************/
void componentDFS(int v, int nComp, int *dim);
/***************************************************************************/

/***************************************************************************/
IloInt findCompCoeff(int dim, int vIndexComp, int nComp);
/***************************************************************************/

/***************************************************************************/
void connectivityDFS(int v, int nComp, int *dim, int *connValue, int pred, int d);
/***************************************************************************/

//DFS that (try to) find a biconnected component and, as usual, define connected component
/***************************************************************************/
void mixedDFS(int v, int nComp, int *dim, int *connValue, int pred, int d);
/***************************************************************************/

/***************************************************************************/
int bendersDFS(int v, int nComp, int *dim);
/***************************************************************************/

/***************************************************************************/
void bendersBFS(int v, int nComp, int *dim);
/***************************************************************************/

//Component BFS that enumerate the component and find and define a tree inside a component
/***************************************************************************/
void componentBFS(int v, int nComp, int *dim);
/***************************************************************************/

//Component BFS that enumerate the component and find and define a tree inside a component
/***************************************************************************/
void mixedBFS(int v, int nComp);
/***************************************************************************/

#endif
