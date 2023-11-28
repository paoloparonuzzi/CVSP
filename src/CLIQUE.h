#ifndef CLIQUE_HEADER
#define CLIQUE_HEADER


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

#include "global_variables.h"

#include "Graph_v4.h"

using namespace std;

/********************************************/
void   clique_load_cplex (graphFF G_bar);
/********************************************/

/***************************************************************************/
void clique_update_cplex (int *clique,graphFF G,graphFF G_bar);
/***************************************************************************/

/***************************************************************************/
int clique_solve_cplex(int *solution,int n);
/***************************************************************************/

/***************************************************************************/
double clique_solve_cplex_edge_plus(double *node_weights_EICP, double *solution,int n);
/***************************************************************************/

/********************************************/
double clique_solve_cplex_edge(int **fixing, double *solution,int n);
/********************************************/

/********************************************/
double clique_solve_cplex_node(int *fixing, double *solution,int n);
/********************************************/

/********************************************/
void clique_free_cplex (graphFF G);
/********************************************/

void clique_fix_node (int v);

void clique_remove_node (int v);

void clique_free_node (int v);




#endif
