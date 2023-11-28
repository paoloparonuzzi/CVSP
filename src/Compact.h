#ifndef COMPACT_local_HEADER
#define COMPACT_local_HEADER


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
#include "global_functions.h"
#include "Graph_v4.h"


using namespace std;

/***************************************************************************/
double compact_model_solve(char *istname,graphFF G,int k,int q,int b);
/***************************************************************************/

/***************************************************************************/
void compact_model_load(graphFF G,int k,int q,int b);
/***************************************************************************/

/***************************************************************************/
void compact_model_free();
/***************************************************************************/


/***************************************************************************/
int ** clique_maximum_exact(graphFF G, int k, int q, int b);
/***************************************************************************/

/***************************************************************************/
int **clique_maximal_heuristic(graphFF G, int k, int q, int b);
/***************************************************************************/

/***************************************************************************/
int * separate_clique_inequalities(graphFF G, int k, int q, int b, double ** vector_xiv);
/***************************************************************************/

#endif
