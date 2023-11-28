//
// Created by paolo on 15/11/18.
//

#ifndef KSEP_CAPA_CONNECTIVITYCUT_H
#define KSEP_CAPA_CONNECTIVITYCUT_H

#include "ilcplex/ilocplex.h"
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

//DFS that (try to) find a biconnected component and, as usual, define connected component
/***************************************************************************/
void userConnectivityDFS(int v, int nComp, int *dim, int *connValue, int pred, int d);
/***************************************************************************/

void ConnectivityCutModel();

#endif //KSEP_CAPA_CONNECTIVITYCUT_H
