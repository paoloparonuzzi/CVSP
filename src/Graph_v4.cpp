

#include "Graph_v4.h"
#include "global_variables.h"

/* Sorting function:

Sorting of n items (integer) according to non-increasing values of their scores (double)
 */
/********************************************************************************************************************************/
void SORT_NON_INCR(int *item,double *score,int n)
/********************************************************************************************************************************/
{
	int salto, i, j, tempItem;
	double tempScore;

	for(salto=n/2; salto>0; salto/=2){
		for(i=salto; i<n; i++){
			for(j=i - salto; j>=0; j-=salto){
				if(score[j]>=score[j + salto]) break;
				tempScore=score[j];
				score[j]=score[j + salto];
				score[j + salto]=tempScore;
				tempItem=item[j];
				item[j]=item[j + salto];
				item[j + salto]=tempItem;
			}
		}
	}
}

/* Sorting function:

Sorting of n items (integer) according to non-increasing values of their scores (double)
 */
/********************************************************************************************************************************/
void SORT_NON_INCR(int *item,double *score, int* node, int n)
/********************************************************************************************************************************/
{
	int salto,i,j,tempItem;
	double tempScore;

	for (salto = n/2; salto > 0; salto /=2)
		for (i = salto; i < n; i++)
			for (j = i-salto; j >= 0; j-=salto) {
				if (score[j] >= score[j+salto]) break;
				tempScore = score[j]; score[j]=score[j+salto]; score[j+salto]=tempScore;
				tempItem = item[j]; item[j]=item[j+salto]; item[j+salto]=tempItem;
			}
}


/* Sorting function:

Sorting of n items (integer) according to non-decreasing values of their scores (double)
 */
/********************************************************************************************************************************/
void SORT_NON_DECR(int *item,double *score,int n)
/********************************************************************************************************************************/
{
	int salto,i,j,tempItem;
	double tempScore;

	for (salto = n/2; salto > 0; salto /=2)
		for (i = salto; i < n; i++)
			for (j = i-salto; j >= 0; j-=salto) {
				if (score[j] <= score[j+salto]) break;
				tempScore = score[j]; score[j]=score[j+salto]; score[j+salto]=tempScore;
				tempItem = item[j]; item[j]=item[j+salto]; item[j+salto]=tempItem;
			}
}

/********************************************************************************************************************************/
void shuffle(int *array, size_t n, int seed)
/********************************************************************************************************************************/
{
	srand(seed);

    if (n > 1)
    {
        size_t i;
        for (i = 0; i < n - 1; i++)
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}

/********************************************************************************************************************************/
graphFF order_vertices_non_increasing_degree(graphFF G,int AMBuilt)
/********************************************************************************************************************************/
{


	int *item=(int*)calloc(G->n,sizeof(int));
	double *score=(double*)calloc(G->n,sizeof(double));
	int *map=(int*)calloc(G->n,sizeof(int));

	for(int i=0;i<G->n;i++){
		item[i]=i;
		score[i]=G->DT[i];
	}

	//	cout << "BEFORE\n";
	//	for(int i=0;i<G->n;i++){ cout << item[i] << "\t";} cout << endl;

	SORT_NON_INCR(item,score,G->n);

	for(int i=0;i<G->n;i++){
		for(int j=0;j<G->n;j++){
			if(i==item[j]){map[i]=j;break;}
		}
	}

	//	cout << "AFTER\n";
	//	for(int i=0;i<G->n;i++){ cout << item[i] << "\t";} cout << endl;
	//	for(int i=0;i<G->n;i++){ cout << map[i] << "\t";} cout << endl;

	int *heads=(int*)calloc(G->m,sizeof(int));
	int *tails=(int*)calloc(G->m,sizeof(int));
	double *weights_arcs=(double*)calloc(G->m,sizeof(double));
	double *weights_nodes=(double*)calloc(G->n,sizeof(double));


	for(int i=0;i<G->m;i++){

		//cout << "old\t\t" << G->T[i] << "\t" << G->H[i] << endl;

		int n_1=map[G->T[i]];
		int n_2=map[G->H[i]];

		//cout << "new\t\t"<< min(n_1,n_2) << "\t" << max(n_1,n_2) << endl;

		tails[i]=min(n_1,n_2);
		heads[i]=max(n_1,n_2);
	}

	for(int i=0;i<G->m;i++){
		weights_arcs[i]=G->P[i];
	}

	for(int i=0;i<G->n;i++){
		weights_nodes[i]=G->W[item[i]];
	}

	int n=G->n;
	int m=G->m;

	deleteGraphFF(G);

	G = buildGraphFF(n,m,heads,tails,weights_nodes,weights_arcs,1);

	//	printGRAPH(GG);
	//	printFS(GG);
	//	printBS(GG);
	//	printAM(GG);

	free(map);
	free(heads);
	free(tails);
	free(weights_arcs);
	free(weights_nodes);
	free(item);
	free(score);

	return(G);

}

/////////////////////////////////////////////////////////////////////////////////
graphFF  order_vertices_non_decreasing_degree(graphFF G,int AMBuilt)
/////////////////////////////////////////////////////////////////////////////////
{

	int *item=(int*)calloc(G->n,sizeof(int));
	double *score=(double*)calloc(G->n,sizeof(double));
	int *map=(int*)calloc(G->n,sizeof(int));


	for(int i=0;i<G->n;i++){
		item[i]=i;
		score[i]=G->DT[i];
	}

	//	cout << "BEFORE\n";
	//	for(int i=0;i<G->n;i++){ cout << item[i] << "\t";} cout << endl;

	SORT_NON_DECR(item,score,G->n);

	for(int i=0;i<G->n;i++){
		for(int j=0;j<G->n;j++){
			if(i==item[j]){map[i]=j;break;}
		}
	}

	//	cout << "AFTER\n";
	//	for(int i=0;i<G->n;i++){ cout << item[i] << "\t";} cout << endl;
	//	for(int i=0;i<G->n;i++){ cout << map[i] << "\t";} cout << endl;


	int *heads=(int*)calloc(G->m,sizeof(int));
	int *tails=(int*)calloc(G->m,sizeof(int));
	double *weights_arcs=(double*)calloc(G->m,sizeof(double));
	double *weights_nodes=(double*)calloc(G->n,sizeof(double));



	for(int i=0;i<G->m;i++){

		//cout << "old\t\t" << G->T[i] << "\t" << G->H[i] << endl;

		int n_1=map[G->T[i]];
		int n_2=map[G->H[i]];

		//cout << "new\t\t"<< min(n_1,n_2) << "\t" << max(n_1,n_2) << endl;

		tails[i]=min(n_1,n_2);
		heads[i]=max(n_1,n_2);
	}

	for(int i=0;i<G->m;i++){
		weights_arcs[i]=G->P[i];
	}

	for(int i=0;i<G->n;i++){
		weights_nodes[i]=G->W[item[i]];
	}

	int n=G->n;
	int m=G->m;

	deleteGraphFF(G);

	G = buildGraphFF(n,m,heads,tails,weights_nodes,weights_arcs,1);

	//	printGRAPH(GG);
	//	printFS(GG);
	//	printBS(GG);
	//	printAM(GG);

	free(map);
	free(heads);
	free(tails);
	free(weights_arcs);
	free(weights_nodes);
	free(item);
	free(score);

	return(G);


}


/********************************************************************************************************************************/
graphFF  order_vertices_order(graphFF G,int AMBuilt,int *order)
/********************************************************************************************************************************/
{

	int *map=(int*)calloc(G->n,sizeof(int));


	for(int i=0;i<G->n;i++){
		for(int j=0;j<G->n;j++){
			if(i==order[j]){map[i]=j;break;}
		}
	}

	//	cout << "AFTER\n";
	//	for(int i=0;i<G->n;i++){ cout << item[i] << "\t";} cout << endl;
	//	for(int i=0;i<G->n;i++){ cout << map[i] << "\t";} cout << endl;

	int *heads=(int*)calloc(G->m,sizeof(int));
	int *tails=(int*)calloc(G->m,sizeof(int));
	double *weights_arcs=(double*)calloc(G->m,sizeof(double));
	double *weights_nodes=(double*)calloc(G->n,sizeof(double));



	for(int i=0;i<G->m;i++){

		//cout << "old\t\t" << G->T[i] << "\t" << G->H[i] << endl;

		int n_1=map[G->T[i]];
		int n_2=map[G->H[i]];

		//cout << "new\t\t"<< min(n_1,n_2) << "\t" << max(n_1,n_2) << endl;

		tails[i]=min(n_1,n_2);
		heads[i]=max(n_1,n_2);
	}

	for(int i=0;i<G->m;i++){
		weights_arcs[i]=G->P[i];
	}

	for(int i=0;i<G->n;i++){
		weights_nodes[i]=G->W[order[i]];
	}

	int n=G->n;
	int m=G->m;

	deleteGraphFF(G);

	G = buildGraphFF(n,m,heads,tails,weights_nodes,weights_arcs,1);

	//	printGRAPH(GG);
	//	printFS(GG);
	//	printBS(GG);
	//	printAM(GG);

	free(map);
	free(heads);
	free(tails);
	free(weights_arcs);
	free(weights_nodes);

	return(G);


}

/********************************************************************************************************************************/
void checks_info(graphFF G,char* istname,int vertex_number,int edge_number,int *tails,int *heads,bool _basic)
/********************************************************************************************************************************/
{



	bool OK_undirected=true;

	if(_basic==false){
		for(int i=0;i<edge_number && OK_undirected;i++){
			for(int j=0;j<edge_number;j++){
				if ((heads[i] == tails[j]) && (tails[i] == heads[j])){
					//						cout << "DIRECTED GRAPH\n";
					//						cout << tails[j] << "\t" << heads[j] << endl;
					//						cout << tails[i] << "\t" << heads[i] << endl;
					OK_undirected=false;
				}
			}
		}
	}

	bool OK_connected=true;

	if(_basic==false){

		double no_path=9999;
		double **weight=new double*[vertex_number];
		double **dist=new double*[vertex_number];
		int **pred=new int*[vertex_number];

		for(int i=0;i<vertex_number;i++){
			weight[i]=new double[vertex_number];
			dist[i]=new double[vertex_number];
			pred[i]=new int[vertex_number];

		}

		for(int i=0;i<vertex_number;i++){
			for(int j=0;j<vertex_number;j++){
				weight[i][j]=no_path;
			}
		}

		for(int i=0;i<edge_number;i++){
			weight[tails[i]][heads[i]]=1;
			weight[heads[i]][tails[i]]=1;
		}

		Floyd_Warshall(vertex_number,weight,dist,pred);


		for(int i=0;i<vertex_number;i++){
			for(int j=i+i;j<vertex_number;j++){
				if(dist[i][j]==no_path){
					OK_connected=false;
				}
			}
		}


		for(int i=0;i<vertex_number;i++){
			delete [] weight[i];
			delete [] dist[i];
			delete [] pred[i];
		}
		delete [] weight;
		delete [] dist;
		delete [] pred;
	}

	double density= ((double) edge_number) / (( vertex_number * (vertex_number-1) ) / 2);

	int min_degree=G->n;
	int max_degree=0;

	for(int i=0;i<G->n;i++){

		if(min_degree>G->DT[i]){min_degree=G->DT[i];}
		if(max_degree<G->DT[i]){max_degree=G->DT[i];}

	}

	ofstream info_SUMMARY("info_graphFFs.txt", ios::app);
	info_SUMMARY << fixed
			<< istname << "\t"
			<< OK_undirected << "\t"
			<< OK_connected << "\t"
			<< vertex_number << "\t"
			<< edge_number << "\t"
			<< density << "\t"
			<< _basic << "\t"
			<< min_degree << "\t"
			<< max_degree << "\t"			
			<< endl;

	info_SUMMARY.close();


	exit(-1);
}



/*-------------------------------------------------------*/
/*                Shortest Path Floyd_Warshall           */
/*-------------------------------------------------------*/

//only positive weights
/********************************************************************************************************************************/
void Floyd_Warshall(int vertex_number,double **weight,double **dist,int **pred)
/********************************************************************************************************************************/
{

	//initialization
	for(int i=0;i<vertex_number;i++){
		for(int j=0;j<vertex_number;j++){
			dist[i][j] = weight[i][j];
			pred[i][j] = i;
		}
	}

	//floyd-warshall
	for(int h=0;h<vertex_number;h++){
		for(int i=0;i<vertex_number;i++){
			for(int j=0;j<vertex_number;j++){
				if (dist[i][j] > dist[i][h] + dist[h][j] /*+0.000001*/) {
					dist[i][j] = dist[i][h] + dist[h][j];
					pred[i][j] = pred[h][j];
				}
			}
		}
	}

}

/*-------------------------------------------------------*/
/*                       Print GraphFF                     */
/*-------------------------------------------------------*/

/********************************************************************************************************************************/
void printGRAPH(graphFF G)
/********************************************************************************************************************************/
{

	cout << "Number of nodes\t"<< G->n <<    endl;
	cout << "Number of arcs\t"<< G->m <<    endl;



	for(int i=0;i<G->m;i++){
		cout << "tail\t" << G->T[i] << "\thead\t" << G->H[i] << "\tweights\t" << G->P[i] << endl;
	}
	cout << endl;

	for(int i=0;i<G->n;i++){
		cout << "Node\t" << i << "\t weight\t" << G->W[i] << "\tneighours total\t" << G->DT[i] << "\t neighours +\t" << G->DP[i] << "\t neighours -\t" << G->DM[i] << endl;

	}

}

/*-------------------------------------------------------*/
/*                Print Adjacency Matrix                 */
/*-------------------------------------------------------*/

/********************************************************************************************************************************/
void printAM(graphFF G)
/********************************************************************************************************************************/
{

	cout << "Adjacency Matrix"<<   endl;

	for(int i=0;i<G->n;i++){
		for(int ii=0;ii<G->n;ii++){
			cout << G->AMatrix[i][ii]+G->AMatrix[ii][i] << " " ;
		}
		cout << endl;
	}
	cout << endl;

}


/*-------------------------------------------------------*/
/*                Print FS and BS                        */
/*-------------------------------------------------------*/

/********************************************************************************************************************************/
void printFS(graphFF G)
/********************************************************************************************************************************/
{

	for(int i=0;i<G->n;i++){
		cout << "Forward star of\t" << i << endl;
		for (int  k = G->NFS[i]; k < G->NFS[i+1]; k++ )
		{
			cout << "Arc\t" << G->AFS[k] << "\ttail\t"<< G->T[G->AFS[k]] << "\thead\t" << G->H[G->AFS[k]] << endl;
		}
	}
}

/********************************************************************************************************************************/
void printBS(graphFF G)
/********************************************************************************************************************************/
{

	for(int i=0;i<G->n;i++){
		cout << "Backward star of\t" << i << endl;
		for (int  k = G->NBS[i]; k < G->NBS[i+1]; k++ )
		{
			cout << "Arc\t" << G->ABS[k] << "\ttail\t"<< G->T[G->ABS[k]] << "\thead\t" << G->H[G->ABS[k]] << endl;
		}
	}

}


/*-------------------------------------------------------*/
/* Build the undirected graphFF from Standard Input        */
/*-------------------------------------------------------*/


/********************************************************************************************************************************/
void deleteGraphFF(graphFF G)
/********************************************************************************************************************************/
{

	if (G->Deleted != 0) return;

	delete[] G->NFS;
	delete[] G->AFS;
	delete[] G->NBS;
	delete[] G->ABS;

	delete[] G->P;
	delete[] G->H;
	delete[] G->T;
	delete[] G->W;

	delete[] G->DM;
	delete[] G->DP;
	delete[] G->DT;

	if(G->AMBuilt==1){
		for(int i=0;i<G->n;i++){delete[] G->AMatrix[i];}
		delete[] G->AMatrix;
		G->AMBuilt=0;
	}

	for(int i=0; i<G->n; i++)
		G->neighb[i].end();
	delete [] G->neighb;

	G->Deleted = 1;

	delete G;
}



/*-------------------------------------------------------------------------------------------------------*/
/*
 * --> Build complementary graphFF ****consider the original graphFF directed***
 * --> it requires AMBuilt built for G
 * --> all weights = 1
 * --> no autoarcs
 */
/*-------------------------------------------------------------------------------------------------------*/

/********************************************************************************************************************************/
graphFF buildComplementaryGraphFF_directed(graphFF G,int AMBuilt)
/********************************************************************************************************************************/
{


	if(G->AMBuilt==0){cout << "AMBuilt not built for G\n";exit(-1);}

	graphFF G_bar;

	int arcs= (  ( G->n * G->n  ) - G->n ) - G->m;
	int nodes= G->n;
	int *heads=new int[arcs];
	int *tails=new int[arcs];
	double *weights_nodes=new double[nodes];
	double *weights_arcs=new double[arcs];

	int counter=0;
	for(int i=0; i<G->n; i++)
	{
		for(int j=0; j<G->n; j++)
		{
			if(G->AMatrix[i][j]==1 || i==j){continue;}
			tails[counter]=i;
			heads[counter]=j;
			counter++;
		}
	}
	if(counter!=(arcs)){cout << "ERROR DIMENSION IN buildComplementaryGraphFF\t" << counter << "\t" << arcs << endl ; exit(-1);}

	for(int a=0;a<arcs;a++){weights_arcs[a]=1;}
	for(int v=0;v<nodes;v++){weights_nodes[v]=1;}

	G_bar=buildGraphFF(nodes,arcs,heads,tails,weights_nodes,weights_arcs,1);

	delete []heads;
	delete []tails;
	delete []weights_nodes;
	delete []weights_arcs;


	return(G_bar);
}




/*-------------------------------------------------------------------------------------------------------*/
/*
 * --> Build complementary graphFF ****consider the original graphFF undirected***
 * --> it requires AMBuilt built for G
 * --> all weights = 1
 * --> no autoarcs
 */
/*-------------------------------------------------------------------------------------------------------*/

/********************************************************************************************************************************/
graphFF buildComplementaryGraphFF_undirected(graphFF G,int AMBuilt)
/********************************************************************************************************************************/
{


	if(G->AMBuilt==0){cout << "AMBuilt not built for G\n";exit(-1);}

	graphFF G_bar;

	int edges= ( ( ( G->n * ( G->n - 1 ) ) ) / 2 ) - G->m;
	int nodes= G->n;
	int *heads=new int[edges];
	int *tails=new int[edges];
	double *weights_nodes=new double[nodes];
	double *weights_arcs=new double[edges];

	int counter=0;
	for(int i=0; i<G->n; i++)
	{
		for(int j=i+1; j<G->n; j++)
		{
			if(G->AMatrix[i][j]==1 || G->AMatrix[j][i]==1){continue;}
			tails[counter]=i;
			heads[counter]=j;
			counter++;
		}
	}
	if(counter!=(edges)){cout << "ERROR DIMENSION IN buildComplementaryGraphFF\t" << counter << "\t" << edges << endl ; exit(-1);}

	for(int e=0;e<edges;e++){weights_arcs[e]=1;}
	for(int v=0;v<nodes;v++){weights_nodes[v]=1;}

	G_bar=buildGraphFF(nodes,edges,heads,tails,weights_nodes,weights_arcs,1);

	delete []heads;
	delete []tails;
	delete []weights_nodes;
	delete []weights_arcs;


	return(G_bar);
}



/*-------------------------------------------------------*/
/*   Build the directed graphFF from Standard Input        */
/*-------------------------------------------------------*/

/********************************************************************************************************************************/
graphFF buildGraphFF(int node_number,int arc_number,int *heads, int *tails,double *weights_node,double *weights_arcs,int AMBuilt)
/********************************************************************************************************************************/
{

	int i, k, j;
	int n;        // nodes
	int m;        // arcs
	int *H;       // H[i] = head of the edge
	int *T;       // T[i] = tail of the edge
	double *P;    // P[i] = weight of the edge
	double *W;    // W[i] = weight of the node

	graphFF G = new GraphFF;

	n=node_number;  // number of nodes
	m=arc_number;  // number of arcs


	/**** Memory Allocation ****/
	H = new int[m];
	T = new int[m];
	P = new double[m];
	W = new double[n];


	/* store the data in the graphFF structure */
	G->H = H;
	G->T = T;
	G->P = P;
	G->n = n;
	G->m = m;
	G->W = W;



	for(int e=0;e<arc_number;e++){
		H[e]=heads[e];
		T[e]=tails[e];
		P[e]=weights_arcs[e];
	}

	for(int v=0;v<node_number;v++){
		W[v]=weights_node[v];
	}


	// BS and FS still to build
	G->LABuilt = 0;
	// momory to be deallocated
	G->Deleted = 0;


	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////


	if(G->LABuilt != 0)
	{ printf("Forward/Backward already built.\n");
	return(G); //
	}
	if(G->Deleted != 0)
	{ printf("structures already deallocated.\n");
	return(G); //
	}

	/* Allocation of NFS, AFS*/
	int * NFS = new int[n+1];
	int * AFS = new int[m];
	int * NBS = new int[n+1];
	int * ABS = new int[m];

	//{ Phase 1: dimensioning FS and BS}
	for(int i=0;i<n;i++){
		NFS[i] = 0;
		NBS[i] = 0;
	}

	for(int e=0;e<m;e++){
		i = T[e];
		NFS[i] = NFS[i] + 1;
		j = H[e];
		NBS[j] = NBS[j] + 1;

	}

	//{ Phase 2: filling NFS and NBS}
	NFS[n ] = m ;
	NBS[n ] = m ;

	for(int i=n-1;i>=0;i--){
		NFS[i] = NFS[i + 1] - NFS[i];
		NBS[i] = NBS[i + 1] - NBS[i];
	}

	//{ Phase 3: filling AFS }
	for(int e=0;e<m;e++){

		i = T[e];
		k = NFS[i];
		AFS[k] = e;
		NFS[i] = NFS[i] + 1;

		j = H[e];
		k = NBS[j];
		ABS[k] = e;
		NBS[j] = NBS[j] + 1;
	}


	//{ Phase 4: restore NFS and NBS}
	for(int i=n-1;i>=1;i--){
		NFS[i] = NFS[i-1];
		NBS[i] = NBS[i-1];
	}
	NFS[0] = 0;
	NBS[0] = 0;

	///////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////

	//	cout << "NFS\n";
	//	for(int i=0;i<n+1;i++){
	//		cout << NFS[i] << endl;
	//	}
	//
	//	cout << "AFS\n";
	//	for(int i=0;i<m;i++){
	//		cout << AFS[i] << endl;
	//	}
	//
	//	cout << "NBS\n";
	//	for(int i=0;i<n+1;i++){
	//		cout << NBS[i] << endl;
	//	}
	//
	//	cout << "ABS\n";
	//	for(int i=0;i<m;i++){
	//		cout << ABS[i] << endl;
	//	}

	/* store in G */
	G->NFS = NFS;
	G->AFS = AFS;

	G->NBS = NBS;
	G->ABS = ABS;

	// structures built
	G->LABuilt=1;

	G->DT=new int[n];  // number of neighbours total
	G->DP=new int[n];  // number of neighbours delta+
	G->DM=new int[n];  // number of neighbours delta-

	for(int i=0;i<G->n;i++){
		G->DT[i]=0;
		G->DP[i]=0;
		G->DM[i]=0;
		for (int  k = G->NFS[i]; k < G->NFS[i+1]; k++ )
		{
			G->DT[i]=G->DT[i]+1;
			G->DP[i]=G->DP[i]+1;
		}
		for (int  k = G->NBS[i]; k < G->NBS[i+1]; k++ )
		{
			G->DT[i]=G->DT[i]+1;
			G->DM[i]=G->DM[i]+1;
		}
	}

	G->neighb = new vector<int>[n];
	for(int i=0; i<n; i++){
		G->neighb[i].resize(static_cast<unsigned long>(G->DT[i]));
		int j=0;
		for(int kk=G->NFS[i]; kk<G->NFS[i+1]; kk++)
		{
			G->neighb[i][j]=G->H[G->AFS[kk]];
			j++;
		}
		for(int kk=G->NBS[i]; kk<G->NBS[i+1]; kk++)
		{
			G->neighb[i][j]=G->T[G->ABS[kk]];
			j++;
		}
		if(j!=G->DT[i]){
			cout << "Error in building new neighborhood" << endl;
		}
	}

	if(AMBuilt==1){

		G->AMBuilt=1;

		G->AMatrix = new int*[n];
		for(int i=0;i<n;i++){
			G->AMatrix[i] = new int[n];
		}

		for(int i=0;i<n;i++){
			for(int ii=0;ii<n;ii++){
				G->AMatrix[i][ii]=0;
			}
		}

		for(int e=0;e<m;e++){
			G->AMatrix[G->T[e]][G->H[e]]=1;
		}

		//		for(int i=0;i<n;i++){
		//			for(int ii=0;ii<n;ii++){
		//				cout << G->AMatrix[i][ii] ;
		//			}
		//			cout << endl;
		//		}
		//		cout << endl;

	}
	else{

		G->AMBuilt=0;
	}

	return(G);
}

int buildCliqueEdgeCover(graphFF& G){

	//STEP 1
	int nClique=0;
	unordered_set<int> W;
	unordered_set<int> V;

	for(int i=0; i<vertex_number; i++){
		//STEP 2
		W.clear();
		for(int j=0; j<i; j++){
			if(G->AMatrix[i][j] == 1 || G->AMatrix[j][i] == 1)
				W.insert(j);
		}
		//STEP 3
		if(W.empty()){
			nClique++;
			G->edgeClique.resize(static_cast<unsigned long>(nClique));
			G->edgeClique[nClique-1].push_back(i);
			continue;
		}
		//STEP 4
		int m=0;
		V.clear();
		while(m<nClique && V!=W){
			bool contained=true;
			int count=0;
			while(contained && count<G->edgeClique[m].size()){
				unordered_set<int>::const_iterator got = W.find(G->edgeClique[m][count]);
				contained = (got != W.end());
				count++;
			}
			if(contained){
				G->edgeClique[m].push_back(i);
				V.insert(G->edgeClique[m].begin(), G->edgeClique[m].end());
			}
			m++;
		}
		//STEP 5
		for(int v : V)
			W.erase(v);

		//STEP 6
		while(!W.empty()){
			int maxSize=-1;
			int indexMaxSize=-1;
			for(int mm=0; mm<nClique; mm++){
				int checkSize=0;
				for(int s : G->edgeClique[mm]){
					unordered_set<int>::const_iterator got = W.find(s);
					if(got != V.end())
						checkSize++;
				}
				if(checkSize>maxSize){
					maxSize=checkSize;
					indexMaxSize=mm;
				}
			}
			nClique++;
			G->edgeClique.resize(static_cast<unsigned long>(nClique));
			for(int s : G->edgeClique[indexMaxSize]){
				unordered_set<int>::const_iterator got = W.find(s);
				if(got != W.end())
					G->edgeClique[nClique-1].push_back(s);
				W.erase(s);
			}
			G->edgeClique[nClique-1].push_back(i);
		}
	}

	//output << nClique << "\t";
	
	//STEP 7
	for(int m=0; m<nClique; m++){
		assert(nClique == G->edgeClique.size());

		bool contained = true;
		for(int i=0; i<G->edgeClique[m].size()-1; i++){
			for(int j=i+1; j<G->edgeClique[m].size(); j++){
				int u = G->edgeClique[m][i];
				int v = G->edgeClique[m][j];
				bool edgeFound = false;
				for(int mm=0; mm<nClique; mm++){
					if(mm != m){
						unordered_set<int> copyClique;
						copyClique.insert(G->edgeClique[mm].begin(), G->edgeClique[mm].end());
						unordered_set<int>::const_iterator got1 = copyClique.find(u);
						unordered_set<int>::const_iterator got2 = copyClique.find(v);
						edgeFound = edgeFound || (got1 != copyClique.end() && got2 != copyClique.end());
						copyClique.clear();
					}
				}
				contained = contained && edgeFound;
			}
		}
		
		if(contained){
			G->edgeClique.erase(G->edgeClique.begin()+m);
			nClique--;
			m--;
		}
	}
	return nClique;
}

/* Check Backward Star */
/********************************************************************************************************************************/
int CheckBS(graphFF G)
/********************************************************************************************************************************/
{

	int i, j, a, errors;
	int n, m;
	int *NBS;
	int *ABS;
	int *T;
	int *H;
	int *A;

	n = G->n;
	m = G->m;
	T = G->T;
	H = G->H;
	NBS = G->NBS;
	ABS = G->ABS;
	A = new int[m];
	for (i=0; i<m; i++) A[i]=1;

	errors = 0;

	for (i=0; i<n; i++)
		for (j=NBS[i]; j < NBS[i+1]; j++)
		{ a = ABS[j];
		A[a] -= 1;
		if ( H[a] != i )
		{ printf("ERROR: Arc %1d in BS(%1d), head=%1d\n",
				a,        i,        H[a] );
		errors += 1;
		}
		}

	for (i=0; i<m; i++)
	{ if ( A[i] != 0 )
	{ printf("ERROR: Arc %1d=(%1d,%1d) counted %1d times\n",
			i, H[i], T[i], (1-A[i]) );
	errors += 1;
	}
	}
	printf("End check Backward Star\n");
	delete[] A;
	return(errors);
}

/* Check Forward Star */
/********************************************************************************************************************************/
int CheckFS(graphFF G)
/********************************************************************************************************************************/
{

	int i, j, a, errors;
	int n, m;
	int *NFS;
	int *AFS;
	int *T;
	int *H;
	int *A;

	n = G->n;
	m = G->m;
	T = G->T;
	H = G->H;
	NFS = G->NFS;
	AFS = G->AFS;
	A = new int[m];
	for (i=0; i<m; i++) A[i]=1;

	errors = 0;

	for (i=0; i<n; i++)
		for (j=NFS[i]; j < NFS[i+1]; j++)
		{ a = AFS[j];
		A[a] -= 1;
		if ( T[a] != i )
		{ printf("ERROR: Arc %1d in FS(%1d), tail=%1d\n",
				a,        i,        T[a] );
		errors += 1;
		}
		}

	for (i=0; i<m; i++)
	{ if ( A[i] != 0 )
	{ printf("ERROR: Arc %1d=(%1d,%1d) counted %1d times\n",
			i, T[i], H[i], (1-A[i]) );
	errors += 1;
	}
	}
	printf("End check Forward Star\n");
	delete[] A;
	return(errors);
}



/********************************************************************************************************************************/
int  random_undirected_graphFF_data_load_density(int nodes,double density,int *tails,int *heads,double *weights_nodes,double *weights_edges,int min_weight,int max_weight)
/********************************************************************************************************************************/
{


	int edges=0;

	if(density > 1.0 ){cout << "Too dense\n\n"; exit(-1);}

	for(int i=0;i<nodes;i++){
		for(int j=i+1;j<nodes;j++){

			float random=(float)rand() / (float)RAND_MAX ;

//			cout << "random\t" << random << endl;
//			cin.get();

			if(random>density){continue;}

			tails[edges]=i;
			heads[edges]=j;
			edges++;

		}
	}

//	cout << "EGDES\n";
//	for(int i=0;i<edges;i++){
//		cout << tails[i]  << "\t" << heads[i] << endl;
//	}


	cout << "density\t"  << density << "\t real density \t" << (double)(edges)/((nodes*(nodes-1))/2.0) << "\t edges \t" << edges << endl;
//	cin.get();

	int random_1;
	if(max_weight>min_weight){
		for(int i=0;i<edges;i++){
			random_1=(int)rand();
			random_1=(random_1%(max_weight-1-min_weight+1));
			weights_edges[i]=(int)random_1;
		}
		for(int i=0;i<nodes;i++){
			random_1=(int)rand();
			random_1=(random_1%(max_weight-1-min_weight+1));
			weights_nodes[i]=(int)random_1;
		}
	}
	else{
		for(int i=0;i<edges;i++){
			weights_edges[i]=max_weight;
		}
		for(int i=0;i<nodes;i++){
			weights_nodes[i]=max_weight;
		}
	}

	return edges;

}


/********************************************************************************************************************************/
void random_undirected_graphFF_data_load(int nodes,int edges,int *tails,int *heads,double *weights_nodes,double *weights_edges,int min_weight,int max_weight)
/********************************************************************************************************************************/
{


	if(edges > ((nodes*(nodes-1))/2) ){cout << "Too many edges\n\n"; exit(-1);}

	int random_1,random_2;
	for(int i=0;i<edges;i++){
		bool _OK= false;
		do{

			random_1=(int)rand();
			random_1=(random_1%(nodes-1-0+1));

			random_2=(int)rand();
			random_2=(random_2%(nodes-1-0+1));

			if(random_1!=random_2){
				_OK=true;
				for(int j=0; j <i; j++){
					if (
							heads[j] == (int)min(random_1,random_2)
							&&
							tails[j] == (int)max(random_1,random_2)
					)
					{
						_OK=false;
					}

				}
			}
		}while(_OK==false);

		heads[i]=(int)min(random_1,random_2);
		tails[i]=(int)max(random_1,random_2);

	}

	if(max_weight>min_weight){
		for(int i=0;i<edges;i++){
			random_1=(int)rand();
			random_1=(random_1%(max_weight-1-min_weight+1));
			weights_edges[i]=(int)random_1;
		}
		for(int i=0;i<nodes;i++){
			random_1=(int)rand();
			random_1=(random_1%(max_weight-1-min_weight+1));
			weights_nodes[i]=(int)random_1;
		}
	}
	else{
		for(int i=0;i<edges;i++){
			weights_edges[i]=max_weight;
		}
		for(int i=0;i<nodes;i++){
			weights_nodes[i]=max_weight;
		}
	}


}


/********************************************************************************************************************************/
int  random_directed_graphFF_data_load_density(int nodes,double density,int *tails,int *heads,double *weights_nodes,double *weights_edges,int min_weight,int max_weight)
/********************************************************************************************************************************/
{


	int edges=0;

	if(density > 1.0 ){cout << "Too dense\n\n"; exit(-1);}

	for(int i=0;i<nodes;i++){
		for(int j=0;j<nodes;j++){

			if(i==j){continue;}

			float random=(float)rand() / (float)RAND_MAX ;

//			cout << "random\t" << random << endl;
//			cin.get();

			if(random>density){continue;}

			tails[edges]=i;
			heads[edges]=j;
			edges++;

		}
	}

//	cout << "EGDES\n";
//	for(int i=0;i<edges;i++){
//		cout << tails[i]  << "\t" << heads[i] << endl;
//	}


	cout << "density\t"  << density << "\t real density \t" << (double)(edges)/((nodes*(nodes-1))/2.0) << "\t edges \t" << edges << endl;
//	cin.get();

	int random_1;
	if(max_weight>min_weight){
		for(int i=0;i<edges;i++){
			random_1=(int)rand();
			random_1=(random_1%(max_weight-1-min_weight+1));
			weights_edges[i]=(int)random_1;
		}
		for(int i=0;i<nodes;i++){
			random_1=(int)rand();
			random_1=(random_1%(max_weight-1-min_weight+1));
			weights_nodes[i]=(int)random_1;
		}
	}
	else{
		for(int i=0;i<edges;i++){
			weights_edges[i]=max_weight;
		}
		for(int i=0;i<nodes;i++){
			weights_nodes[i]=max_weight;
		}
	}

	return edges;

}


/********************************************************************************************************************************/
void random_directed_graphFF_data_load(int nodes,int arcs,int *tails,int *heads,double *weights_nodes,double *weights_arcs,int min_weight,int max_weight)
/********************************************************************************************************************************/
{

	if(arcs > ((nodes*(nodes-1))) ){cout << "Too many arcs\n\n"; exit(-1);}


	int random_1,random_2;
	for(int i=0;i<arcs;i++){
		bool _OK= false;
		do{

			random_1=(int)rand();
			random_1=(random_1%(nodes-1-0+1));

			random_2=(int)rand();
			random_2=(random_2%(nodes-1-0+1));

			if(random_1!=random_2){
				_OK=true;
				for(int j=0; j <i; j++){
					if (
							heads[j] == (int)random_1
							&&
							tails[j] == (int)random_2
					)
					{
						_OK=false;
					}

				}
			}
		}while(_OK==false);

		heads[i]=(int)random_1;
		tails[i]=(int)random_2;

	}

	if(max_weight>min_weight){
		for(int i=0;i<arcs;i++){
			random_1=(int)rand();
			random_1=(random_1%(max_weight-1-min_weight+1));
			weights_arcs[i]=(int)random_1;
		}
		for(int i=0;i<nodes;i++){
			random_1=(int)rand();
			random_1=(random_1%(max_weight-1-min_weight+1));
			weights_nodes[i]=(int)random_1;
		}
	}
	else{
		for(int i=0;i<arcs;i++){
			weights_arcs[i]=max_weight;
		}
		for(int i=0;i<nodes;i++){
			weights_nodes[i]=max_weight;
		}
	}


}

/********************************************************************************************************************************/
void ReadDIMACSFile_C(char* fileName,int *nodes,int *edges,int *tails,int *heads)
/********************************************************************************************************************************/
{


	/* file di input */
	FILE* inputFile;

	/* matrice di adiacenza implementata come un unico vettore */
	//	int* adjacencyMatrix;

	/* dimensione della matrice di adiacenza */
	//int matrixDimension;

	/* indici generici */
	int i, /*j,*/ v1, v2;

	/* indice per la costruzione del vettore degli adiacenti */
	//int indexV;

	/* controllo che ogni vertice abbia un adiacente */
	//int check;

	char c;
	char str[255];


	/* apertura del file dell'istanza */
	if( (inputFile=fopen(fileName, "r")) == NULL ){
		printf("File non esistente: %s\n", fileName);
		exit(2);
	}

	/* lettura intestazione */
	fscanf(inputFile, "%s", &c);
	while (c=='c'){
		printf("\n read comment line");
		fgets(str,255,inputFile);
		fscanf(inputFile, "%s", &c);
	}

	int vertexNumber,edgeNumber;

	/* lettura delle dimensioni del grafo */
	fscanf(inputFile, "%*s %d %d", &vertexNumber, &edgeNumber);
	//ignora "p edge" oppure "e edge"

	*nodes=vertexNumber;
	*edges=edgeNumber;

	if( (vertexNumber<2) || (edgeNumber<1) ){
		puts("readInputGraphFF: Errore nelle dimensioni del problema.");
		exit(3);
	}

	printf("\nIstanza: %s\n\n", fileName);
	printf("Numero vertici = %d\nNumero lati = %d\n\n", vertexNumber, edgeNumber);


	/* lettura dei lati */
	for(i=0; i<edgeNumber; i++){
		if(fscanf(inputFile, "%*s %d %d", &v1, &v2) !=2){
			//ignora "e" all'inizio di ogni riga
			puts("readInputGraphFF: Errore nei dati in input durante la lettura dei lati.");
			exit(4);
		}

		if( v1==v2 || v1<1 || v1>vertexNumber || v2<1 || v2>vertexNumber ){
			puts("readInputGraphFF: Errore nei dati in input: incoerenza nei vertici.");
			exit(5);
		}

		tails[i]=v1-1;
		heads[i]=v2-1;


	}//for

	/* controllo che il file sia terminato, consumando i caratteri bianchi */
	if(fscanf(inputFile, "%*s") != EOF){
		puts("readInputGraphFF: Errore al termine della lettura dell'input: il file contiene ulteriori dati.");
		exit(6);
	}
	fclose(inputFile);


}//readInputGraphFF


/********************************************************************************************************************************/
void ReadDIMACSFile(char * inFile,int *nodes,int *edges,int *tails,int *heads,bool checks)
/********************************************************************************************************************************/
{


	cout << "DIMACS INSTANCE\n";
	cout << "\n-> original graphFF\n";

	ifstream in(inFile);
	if(!in)
	{
		cout << "File could not be opened. " << endl;
		exit(1);
	}


	char *buf_1=new char[5000];

	string str1;
	string str2;
	string str3;

	str2='e';
	str3='p';

	for( string line; getline( in, line ); )
	{

		istringstream buffer(line);

		buffer >> str1;

		//cout << "str1\t"<< str1 << endl;

		if(!str3.compare(str1)){

			buffer >> buf_1;
			//			cout << "buf_1\t" << buf_1 << endl;

			buffer >> *nodes;

			buffer >> *edges;

			break;
		}

	}

	/*
	cout << "\n\nPrint Info and exit....\n";

		ofstream info_SUMMARY("info_graphFFs.txt", ios::app);
		info_SUMMARY << fixed
				<< inFile << "\t"
				<< *nodes << "\t"
				<< *edges << "\n";
		info_SUMMARY.close();
	exit(-1);
	 */

	bool OK_undirected=true;
	bool OK_connected=true;

	int i=0;
	for( string line; getline( in, line ); )
	{

		istringstream buffer(line);
		//		cout << line << endl;

		buffer >> str1;
		//				cout << "str1\t"<< str1 << endl;

		if(!str2.compare(str1)){

			//			cout << "edge number\t" << i << endl;
			int a,b;
			buffer >> a;
			buffer >> b;
			//			cout << a << "\t" << b << endl;
			//			cin.get();

			tails[i]=a-1;
			heads[i]=b-1;

			if(checks)
			{
				for(int j=0;j<i;j++){

					if ((heads[i] == tails[j]) && (tails[i] == heads[j])){
						//						cout << "DIRECTED GRAPH\n";
						//						cout << tails[j] << "\t" << heads[j] << endl;
						//						cout << tails[i] << "\t" << heads[i] << endl;
						OK_undirected=false;
					}

				}
			}
			//			cout << heads[i] << "\t" << tails[i] << endl;
			//			cin.get();

			if(i==(*edges-1)){
				break;
			}

			i++;
		}

	}

	if(checks){

		double no_path=999999;
		double **weight=new double*[*nodes];
		double **dist=new double*[*nodes];
		int **pred=new int*[*nodes];

		for(int i=0;i<*nodes;i++){
			weight[i]=new double[*nodes];
			dist[i]=new double[*nodes];
			pred[i]=new int[*nodes];

		}

		for(int i=0;i<*nodes;i++){
			for(int j=0;j<*nodes;j++){
				weight[i][j]=no_path;
			}
		}

		for(int i=0;i<*edges;i++){
			weight[tails[i]][heads[i]]=1;
			weight[heads[i]][tails[i]]=1;
		}

		Floyd_Warshall(*nodes,weight,dist,pred);


		for(int i=0;i<*nodes;i++){
			for(int j=i+i;j<*nodes;j++){
				if(dist[i][j]==no_path){
					OK_connected=false;
				}
			}
		}


		for(int i=0;i<*nodes;i++){
			delete [] weight[i];
			delete [] dist[i];
			delete [] pred[i];
		}
		delete [] weight;
		delete [] dist;
		delete [] pred;


		cout << "\n\nPrint Info and exit....\n";

		ofstream info_SUMMARY("info_graphFFs.txt", ios::app);
		info_SUMMARY << fixed
				<< inFile << "\t"
				<< OK_undirected << "\t"
				<< OK_connected << "\t"
				<< *nodes << "\t"
				<< *edges << "\n";
		info_SUMMARY.close();


		exit(-1);
	}
	/////////////////////////////////////////////////////////////////////////////////

	delete[] buf_1;


	in.close();

	cout << "read!\n";

	return;

}

/********************************************************************************************************************************/
void funzione_crea_figura_grafo_dot(char *nome_file,int n_vertex,int n_layer,int *n_out,int *nxt_out,int *head)
/********************************************************************************************************************************/
{

	char *probname=new char[1000];
	sprintf(probname, "%s_dot_graphFF.txt",nome_file);

	int i, j;
	int n;
	ofstream myfile;
	myfile.open (probname);
	myfile << "digraphFF G {\n";
	//myfile << "\tsize = \"75, 100\";\n";
	myfile << "\tranksep = 30; size = \"7.5, 10\";\n";
	myfile << "\trankdir = LR;\n";
	//////////////////////////////////////////////////////////

	int dummy=0;
	for (i = 0; i < n_vertex - 1; i++){
		n = n_out[i];
		if (n == 0) continue;
		myfile << "\t"; myfile << i;
		myfile << " -> { ";
		for (j = 0; j < n; j++){
			if(j==n-1){myfile << nxt_out[dummy] << "";}
			else{myfile << nxt_out[dummy] << "; ";}
			dummy++;
		}
		myfile << "}\n";
	}
	//////////////////////////////////////////////////////////
	for (i = 0; i < n_layer; i++){
		myfile << "\t{ rank = same";
		for (j = 0; j < n_vertex ; j++) 
			if (head[j] == i){
				myfile << "; "; myfile << j;
			}
		myfile << ";}\n";
	}
	//////////////////////////////////////////////////////////
	myfile << "}\n";
	myfile.close();
	delete []probname;
}


//vuole in pasto il nome del file, le coordinate dei punti e la matrice delle distanze se i punti sono collegati (flag_archi = 1 disegna anche gli archi)
/********************************************************************************************************************************/
void funzione_crea_figura_grafo_gml(char *nome_file,double *vet_x,double *vet_y,int numero_punti,double coef_scala,double **matrice_Dist,int flag_archi)
/********************************************************************************************************************************/
{

	int i,j;

	char *probname=new char[1000];
	sprintf(probname, "%s.gml",nome_file);
	//cout << probname << endl;
	ofstream grafo(probname);

	char *buf_=new char[1000];
	char *buf_COOR=new char[1000];

	grafo  << "graphFF  [ hierarchic  1  directed  1 \n\n" << endl;

	for (i = 0; i < numero_punti; ++i){
		//solo puntini
		sprintf(buf_COOR, " %d ",i);
		sprintf(buf_, "node  [ id  %d  graphFFics  [ x %f   y %f  w 11 h 11 type \"roundrectangle\" ]  LabelGraphFFics  [ text  %s fontSize  7 ]  ]  ",i,vet_x[i]*coef_scala,vet_y[i]*coef_scala,buf_COOR);
		//scatole con le coordinare
		//sprintf(buf_COOR, "\" %d(%.1f,%.1f)  \"",i,vet_x[i],vet_x[i]);
		//sprintf(buf_, "node  [ id  %d  graphFFics  [ x %f   y %f  w 28 h 8 type \"roundrectangle\" ]  LabelGraphFFics  [ text  %s fontSize  4 ]  ]  ",i,vet_x[i]*coef_scala,vet_y[i]*coef_scala,buf_COOR);
		grafo  << buf_ << endl;
	}

	if(flag_archi==1){
		for (i = 0; i < numero_punti; ++i){
			for (j = 0; j < numero_punti; ++j){
				if(matrice_Dist[i][j]>0){
					//cout << vet_x[i] << " " << vet_y[i] << " "<< vet_x[j] << " " << vet_x[j] << endl; 
					sprintf(buf_, " edge   [ source  %d  target  %d  graphFFics  [ targetArrow \"delta\" Line  [ point  [ x %.2f  y %.2f  ]  point  [ x %.2f  y %.2f  ]  ]  ]  ]  ",i,j,vet_x[i]*coef_scala,vet_y[i]*coef_scala,vet_x[j]*coef_scala,vet_y[j]*coef_scala);
					grafo  << buf_ << endl;
				}
			}
		}
	}

	grafo  << "\n] \n\n" << endl;

	grafo.close();
	delete [] probname;
	delete [] buf_;
	delete [] buf_COOR;

	//cin.get();


}

//vuole in pasto il nome del file, le coordinate dei punti e la matrice delle distanze se i punti sono collegati (flag_archi = 1 disegna anche gli archi)
/********************************************************************************************************************************/
void funzione_crea_figura_grafo_gml_2(char *nome_file,double *vet_x,double *vet_y,int numero_punti,double coef_scala,double **matrice_Dist,int flag_archi)
/********************************************************************************************************************************/
{

	int i,j;

	char *probname=new char[1000];
	sprintf(probname, "%s.gml",nome_file);
	//cout << probname << endl;
	ofstream grafo(probname);

	char *buf_=new char[1000];
	char *buf_COOR=new char[1000];

	grafo  << "graph  [ hierarchic  1  directed  1 \n\n" << endl;

	for (i = 0; i < numero_punti; ++i){
		//solo puntini
		sprintf(buf_COOR, " %d ",i+1);
		sprintf(buf_, "node  [ id  %d  graphics  [ x %f   y %f  w 11 h 11 type \"roundrectangle\" ]  LabelGraphics  [ text  %s fontSize  7 ]  ]  ",i,vet_x[i]*coef_scala,vet_y[i]*coef_scala,buf_COOR);
		//scatole con le coordinare
		//sprintf(buf_COOR, "\" %d(%.1f,%.1f)  \"",i,vet_x[i],vet_x[i]);
		//sprintf(buf_, "node  [ id  %d  graphFFics  [ x %f   y %f  w 28 h 8 type \"roundrectangle\" ]  LabelGraphFFics  [ text  %s fontSize  4 ]  ]  ",i,vet_x[i]*coef_scala,vet_y[i]*coef_scala,buf_COOR);
		grafo  << buf_ << endl;
	}

	if(flag_archi==1){
		for (i = 0; i < numero_punti; ++i){
			for (j = i+1; j < numero_punti; ++j){
				if(matrice_Dist[i][j] >0 || matrice_Dist[j][i] >0){
					//cout << vet_x[i] << " " << vet_y[i] << " "<< vet_x[j] << " " << vet_x[j] << endl; 
					//sprintf(buf_, " edge   [ source  %d  target  %d  graphFFics  [ targetArrow \"delta\" Line  [ point  [ x %.2f  y %.2f  ]  point  [ x %.2f  y %.2f  ]  ]  ]  ]  ",i,j,vet_x[i]*coef_scala,vet_y[i]*coef_scala,vet_x[j]*coef_scala,vet_y[j]*coef_scala);
					sprintf(buf_, " edge   [ source  %d  target  %d  graphics  [ fill	\"#000000\" ]  ]",i,j);
					grafo  << buf_ << endl;
				}
			}
		}
	}

	grafo  << "\n] \n\n" << endl;

	grafo.close();
	delete [] probname;
	delete [] buf_;
	delete [] buf_COOR;

	//cin.get();


}

int indexx(int i, int j)
{return j + k * i;}

/********************************************************************************************************************************/
void funzione_crea_figura_grafo_gml_3(char *nome_file,double *vet_x,double *vet_y,int numero_punti,double coef_scala,double **matrice_Dist,int flag_archi, vector<int> xSol)
/********************************************************************************************************************************/
{

	int i,j;

	char *probname=new char[1000];
	sprintf(probname, "%s_sol.gml",nome_file);
	//cout << probname << endl;
	ofstream grafo(probname);

	char *buf_=new char[1000];
	char *buf_COOR=new char[1000];

	grafo  << "graph  [ hierarchic  1  directed  1 \n\n" << endl;

	vector<int> punti;
	punti.resize(static_cast<unsigned long>(vertex_number));
	for(i=0; i<vertex_number; i++){
		punti[i]=-1;
	}
	for(i=0; i<vertex_number; i++){
		for(j=0; j<k; j++){
			if(xSol[indexx(i,j)]>0.5){
				punti[i] = j;
			}
		}
	}


	for (i = 0; i < numero_punti; ++i){
		//solo puntini
		sprintf(buf_COOR, " %d ",i+1);
		if(punti[i] == -1){
			sprintf(buf_,
					"node  [ id  %d  graphics  [ x %f   y %f  w 11 h 11 type \"roundrectangle\" fill \"#FF0000\" ]  LabelGraphics  [ text  %s fontSize  7 ]  ]  ",
					i, vet_x[i]*coef_scala, vet_y[i]*coef_scala, buf_COOR);
		}
		else{
			sprintf(buf_,
					"node  [ id  %d  graphics  [ x %f   y %f  w 11 h 11 type \"roundrectangle\" fill \"#CCCCFF\" ]  LabelGraphics  [ text  %s fontSize  7 ]  ]  ",
					i, vet_x[i]*coef_scala, vet_y[i]*coef_scala, buf_COOR);
		}
		//scatole con le coordinare
		//sprintf(buf_COOR, "\" %d(%.1f,%.1f)  \"",i,vet_x[i],vet_x[i]);
		//sprintf(buf_, "node  [ id  %d  graphFFics  [ x %f   y %f  w 28 h 8 type \"roundrectangle\" ]  LabelGraphFFics  [ text  %s fontSize  4 ]  ]  ",i,vet_x[i]*coef_scala,vet_y[i]*coef_scala,buf_COOR);
		grafo  << buf_ << endl;
	}

	if(flag_archi==1){
		for (i = 0; i < numero_punti; ++i){
			for (j = i+1; j < numero_punti; ++j){
				if(matrice_Dist[i][j] >0 || matrice_Dist[j][i] >0){
					if(punti[i]==-1 || punti[j]==-1)
						sprintf(buf_, " edge   [ source  %d  target  %d  graphics  [ fill	\"#000000\" style \"dotted\"]  ]",i,j);
					//cout << vet_x[i] << " " << vet_y[i] << " "<< vet_x[j] << " " << vet_x[j] << endl;
					//sprintf(buf_, " edge   [ source  %d  target  %d  graphFFics  [ targetArrow \"delta\" Line  [ point  [ x %.2f  y %.2f  ]  point  [ x %.2f  y %.2f  ]  ]  ]  ]  ",i,j,vet_x[i]*coef_scala,vet_y[i]*coef_scala,vet_x[j]*coef_scala,vet_y[j]*coef_scala);
					else
						sprintf(buf_, " edge   [ source  %d  target  %d  graphics  [ fill	\"#000000\" ]  ]",i,j);
					grafo  << buf_ << endl;
				}
			}
		}
	}

	grafo  << "\n] \n\n" << endl;

	grafo.close();
	delete [] probname;
	delete [] buf_;
	delete [] buf_COOR;

	//cin.get();


}

/********************************************************************************************************************************/
void read_weights(graphFF G,graphFF G_bar,char * inFile)
/********************************************************************************************************************************/

{


	ifstream in(inFile);
	if(!in)
	{
		cout << "File could not be opened. " << endl;
		exit(1);
	}

	for(int i=0;i<G->n;i++){
		in >> G->W[i];
		G_bar->W[i]=G->W[i];;
	}

	cout << "\nWEIGHTS\n";
	for(int i=0;i<G->n;i++){
		cout << "vertex\t" << i << "\t" << G->W[i]  << endl;
	}
	cout << endl;


	in.close();
}

