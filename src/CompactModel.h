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

void LP_CompactModel();
void compactEdgeModel();
void CompactModelToCompare();
void CompactModelWithClique();
void minMaxC_CompactClique();
void minMaxC_CompactEdge();
void minMaxC_CompactFollower();
