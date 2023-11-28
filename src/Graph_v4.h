#ifndef graphFFV4_HEADER
#define graphFFV4_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <unordered_set>
#include <ilconcert/iloexpression.h>


using namespace std;

/*---------------------------------------------*/
// Struct to contain a graphFF.
// GraphFF ---> is the struct kind declared
// graphFF ---> is the pointer kind to GraphFF
/*---------------------------------------------*/
typedef struct {
	int n;        // node number
	int m;        // arcs number
	double *W;    // W[i] = weight of the node
	int *H;       // H[i] = head  of the arc
	int *T;       // T[i] = tail of the arc
	double *P;    // P[i] = weight of the arc
	int *NFS;     // NFS[i]: begin of FS(i) in AFS
	int *AFS;     // arcs index ordered for FS
	int *NBS;     // NBS[i]: begin of BS(i) in ABS
	int *ABS;     // arcs index ordered for BS
	int LABuilt;  // structure built check
	int Deleted;  // structures deleted check
	int **AMatrix;
	int AMBuilt;  // Adjacency Matrix structure built check
	int *DT;  // number of neighbours total
	int *DP;  // number of neighbours delta+
	int *DM;  // number of neighbours delta-
	vector<int> *neighb; //vector with all neighbour
	vector<vector<int>> edgeClique;
} GraphFF, *graphFF;

/********************************************************************************************************************************/
int buildCliqueEdgeCover(graphFF& G);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
int CheckBS(graphFF G);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
int CheckFS(graphFF G);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void printAM(graphFF G);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void printFS(graphFF G);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void printBS(graphFF G);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void printGRAPH(graphFF G);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void deleteGraphFF(graphFF G);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
graphFF buildGraphFF(int node_number,int arc_number,int *heads, int *tails,double *weights_node,double *weights_arcs,int AMBuilt);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void SORT_NON_DECR(int *item,double *score,int n);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void SORT_NON_INCR(int *item,double *score,int n);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void shuffle(int *array, size_t n, int seed);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
graphFF  order_vertices_non_increasing_degree(graphFF G,int AMBuilt);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
graphFF  order_vertices_non_decreasing_degree(graphFF G,int AMBuilt);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
graphFF  order_vertices_order(graphFF G,int AMBuilt,int *order);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
int  random_directed_graphFF_data_load_density(int nodes,double density,int *tails,int *heads,double *weights_nodes,double *weights_edges,int min_weight,int max_weight);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void random_directed_graphFF_data_load(int nodes,int arcs,int *tails,int *heads,double *weights_nodes,double *weights_arcs,int min_weight,int max_weight);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
int  random_undirected_graphFF_data_load_density(int nodes,double density,int *tails,int *heads,double *weights_nodes,double *weights_edges,int min_weight,int max_weight);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void random_undirected_graphFF_data_load(int nodes,int edges,int *tails,int *heads,double *weights_nodes,double *weights_edges,int min_weight,int max_weight);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void ReadDIMACSFile(char * inFile,int *nodes,int *edges,int *tails,int *heads,bool checks);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void ReadDIMACSFile_C(char* fileName,int *nodes,int *edges,int *tails,int *heads);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void Floyd_Warshall(int vertex_number,double **weight,double **dist,int **pred);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
graphFF buildComplementaryGraphFF_undirected(graphFF G,int AMBuilt);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
graphFF buildComplementaryGraphFF_directed(graphFF G,int AMBuilt);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void read_weights(graphFF G,graphFF G_bar,char * inFile);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void checks_info(graphFF G,char* istname,int vertex_number,int edge_number,int *tails,int *heads,bool _basic);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void funzione_crea_figura_grafo_dot(char *nome_file,int n_vertex,int n_layer,int *n_out,int *nxt_out,int *head);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void funzione_crea_figura_grafo_gml(char *nome_file,double *vet_x,double *vet_y,int numero_punti,double coef_scala,double **matrice_Dist,int flag_archi);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void funzione_crea_figura_grafo_gml_2(char *nome_file,double *vet_x,double *vet_y,int numero_punti,double coef_scala,double **matrice_Dist,int flag_archi);
/********************************************************************************************************************************/

/********************************************************************************************************************************/
void funzione_crea_figura_grafo_gml_3(char *nome_file,double *vet_x,double *vet_y,int numero_punti,double coef_scala,double **matrice_Dist,int flag_archi, vector<int> xSol);
/********************************************************************************************************************************/


#endif