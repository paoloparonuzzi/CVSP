

#include "global_functions.h"


/***************************************************************************/
void SORT_NON_INCR_INT(int *item,int *score,int n)
/***************************************************************************/
{
	int salto,i,j,tempItem;
	int tempScore;

	for (salto = n/2; salto > 0; salto /=2)
		for (i = salto; i < n; i++)
			for (j = i-salto; j >= 0; j-=salto)
			{
				if (score[j] >= score[j+salto]) break;
				tempScore = score[j];
				score[j]=score[j+salto];
				score[j+salto]=tempScore;
				tempItem = item[j];
				item[j]=item[j+salto];
				item[j+salto]=tempItem;
			}
}

/***************************************************************************/
void kp_load_cplex(int n_item, double C, double *weights)
/***************************************************************************/
{

	int i;



	env_kp=CPXopenCPLEX(&status);
	if(status!=0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	lp_kp=CPXcreateprob(env_kp,&status,"KP");
	if(status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}


	CPXchgobjsen(env_kp,lp_kp,CPX_MAX);


	//variables
	ccnt=n_item;
	obj=(double*) calloc(ccnt,sizeof(double));
	lb=(double*) calloc(ccnt,sizeof(double));
	ub=(double*) calloc(ccnt,sizeof(double));
	char *_ctype=(char*) calloc(ccnt,sizeof(char));

	for(i=0; i<ccnt; i++)
	{
		obj[i]=(double)0.0;
		lb[i]=0.0;
		ub[i]=1.0;
		_ctype[i]='B';
	}

	status=CPXnewcols(env_kp,lp_kp,ccnt,obj,lb,ub,_ctype,NULL);
	if(status!=0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	//constraints
	rcnt=1;
	nzcnt=n_item;
	rhs=(double*) calloc(rcnt,sizeof(double));
	sense=(char*) calloc(rcnt,sizeof(double));

	rhs[0]=(double)C;		
	sense[0]='L';


	rmatbeg=(int*) calloc(rcnt,sizeof(int));
	rmatind=(int*) calloc(nzcnt,sizeof(int));
	rmatval=(double*) calloc(nzcnt,sizeof(double));



	for(i=0; i<n_item; i++)
	{
		rmatval[i]=(double)weights[i];
		rmatind[i]=i;
	}

	rmatbeg[0]=0;

	status=CPXaddrows(env_kp,lp_kp,0,rcnt,nzcnt,rhs,sense,rmatbeg,rmatind,rmatval,NULL,NULL);
	if(status!=0)
	{
		printf("error in CPXaddrows\n");
		exit(-1);
	}

	free(rmatbeg);
	free(rmatval);
	free(rmatind);
	free(rhs);
	free(sense);


	// writing the problem in a .lp format for control
	//	status=CPXwriteprob(env_kp,lp_kp,"KP.lp",NULL);
	//	if(status!=0) {printf("error in CPXwriteprob\n");	exit(-1);}
	//	cin.get();

	free(obj);
	free(lb);
	free(ub);
	free(_ctype);

}

/***************************************************************************/
double kp_solve_cplex(int n_item,double *profits, double *solution)
/***************************************************************************/
{

	int i;

#ifdef CPLEX_OUTPUT
	CPXsetintparam (env_kp, CPX_PARAM_SCRIND, CPX_ON);
#endif

	// * Set relative tolerance *
	status = CPXsetdblparam (env_kp, CPX_PARAM_EPAGAP, 0.0);
	if (status)
	{
		fprintf (stderr, "error for CPX_PARAM_EPAGAP\n");
	}

	// * Set relative tolerance *
	status = CPXsetdblparam (env_kp, CPX_PARAM_EPGAP, 0.0);
	if (status)
	{
		fprintf (stderr, "error for CPX_PARAM_EPGAP\n");
	}


	// * Set mip tolerances integrality *
	status = CPXsetdblparam (env_kp, CPX_PARAM_EPINT, 0.0);
	if (status)
	{
		fprintf (stderr, "error for CPX_PARAM_EPINTP\n");
	}

	// * Set mip tolerances integrality *
	status = CPXsetdblparam (env_kp, CPX_PARAM_EPINT, 0.0);
	if (status)
	{
		fprintf (stderr, "error for CPX_PARAM_EPINTP\n");
	}

	// * Set Feasibility tolerance *
	status = CPXsetdblparam (env_kp, CPX_PARAM_EPRHS, 1e-9);
	if (status)
	{
		fprintf (stderr, "error for CPX_PARAM_EPRHS\n");
	}

	// * Emphasizes precision in numerically unstable or difficult problems *
	//status = CPXsetdblparam (env_kp, CPX_PARAM_NUMERICALEMPHASIS, 1);
	//if (status)
	//{
	//	fprintf (stderr, "error for CPX_PARAM_NUMERICALEMPHASIS\n");
	//}


	int *ind = (int*) malloc(sizeof(int) * n_item);
	double *d = (double*) malloc(sizeof(double)* n_item);
	for (i = 0; i < n_item; i++)
	{
		d[i] = profits[i];
		ind[i] = i;
	}

	status = CPXchgobj(env_kp, lp_kp,n_item, ind, d);
	if (status != 0)
	{
		printf("error in CPXchgobj\n");
		exit(-1);
	}

	free(ind);
	free(d);

#ifdef writer_KP
	status=CPXwriteprob(env_kp,lp_kp,"kp.lp",NULL);
	if(status!=0) {
		printf("error in CPXwriteprob\n");
		exit(-1);

	}
#endif

	status=CPXmipopt(env_kp,lp_kp);
	if(status!=0)
	{
		printf("error in CPXmipopt\n");
		exit(-1);
	}


	//getting the solution
	double *_x=(double*) calloc(n_item,sizeof(double));
	double _objval_p=0;

	status=CPXgetmipx(env_kp,lp_kp,_x,0,n_item-1);
	if(status!=0)
	{
		printf("error in CPXgetmipx\n");
		exit(-1);
	}
	status=CPXgetmipobjval(env_kp,lp_kp,&_objval_p);
	if(status!=0)
	{
		printf("error in CPXgetmipobjval\n");
		exit(-1);
	}

#ifdef print_ist_sol	
	cout << "\n\nOBJ_KP01 ->\t " << _objval_p << endl;
#endif

	for (i = 0; i < n_item; i++)
	{
		solution[i]=(int)(_x[i]+0.5);
	}

	free(_x);

#ifdef print_ist_sol
	cout << "\n\nSolution\n";
	for (i = 0; i < n_item; i++)
	{
		cout << solution[i] << "";
	}
	cout << endl;
	cout << endl;
#endif


	return _objval_p;

}

// Clique algorithm for separation of clique path constraint (first disjoint clique)
/***************************************************************************/
double findGoodClique(int uIndex, int vIndex, graphFF G, IloNumArray zStar, IloIntArray goodClique)
/***************************************************************************/
{
	
	int in=uIndex;
	int currMax= 0;
	int predMax=-1;
	double best=-1.0;
	double cliqueWeight=0.0;

	for(int i=0; i<vertex_number; i++){
		goodVertex[i]=0;
		goodClique[i]=0;
	}
	goodVertex[vIndex]=-vertex_number;
	
	double randomSeed = 0;//((double) rand() / RAND_MAX);
	double check;
	while(currMax>predMax){
		goodClique[in]=1;
		cliqueWeight=cliqueWeight+zStar[in];
		predMax=currMax;

		for(int i=0; i<vertex_number; i++)
			goodVertex[i]=goodVertex[i]+G->AMatrix[in][i]+G->AMatrix[i][in];

		for(int i=0; i<vertex_number; i++){
			
			check = zStar[i]*randomSeed+G->DT[i]*(1-randomSeed);
			
			if(goodVertex[i]>predMax && check > best){				
				currMax = goodVertex[i];
				in = i;
				best = check;
			}
			//cout << best << endl;
		}
		best=-1.0;
	}
	return cliqueWeight;
}

//Procedure used to perform DFS and in the same time to build the constraint expression to be added to the model
/***************************************************************************/
void lazyArrayDFS(IloExpr *expr, IloIntVarArray x, int v, int nComp, int *dim)
/***************************************************************************/
{
	*dim=*dim+1;
	compVertices[nComp].push_back(v);
	*expr += cap*x[v] + (1.0 - cap)*x[v];

	//as soon as the capacity is violated, store the expr to be lifted and added later to the model 
	if(*dim > cap && statusDFS[v]==-1){
		int index =static_cast<int>(*dim - cap - 1);
		capConstrArray[index] = x[0] - x[0];		
		capConstrArray[index] += *expr;
	}

	statusDFS[v] = nComp;

	shuffle(&G->neighb[v][0], &G->neighb[v][G->DT[v]], default_random_engine());
	for(int i=0; i<G->DT[v]; i++){
		int vNext=G->neighb[v][i];
		if(statusDFS[vNext] == -1){
			*expr += cap*x[v];
			lazyArrayDFS(&*expr, x, vNext, nComp, &*dim);
		}
	}
}

//Procedure used to perform BFS and in the same time to build the constraint expression to be added to the model
/***************************************************************************/
void lazyArrayBFS(IloExpr *expr, IloIntVarArray x, int v, int nComp, int *dim)
/***************************************************************************/
{
	// Create a queue for BFS
	list<int> queue;

	// Mark the current node as visited and enqueue it
	visited[v] = true;
	*dim=*dim+1;
	*expr += cap*x[v] + (1.0 - cap)*x[v];
	compVertices[nComp].push_back(v);
	statusDFS[v] = nComp;
	queue.push_back(v);

	while(!queue.empty())
	{
		// Dequeue a vertex from queue and print it
		int s = queue.front();
		queue.pop_front();

		// Get all adjacent vertices of the dequeued
		// vertex s. If a adjacent has not been visited,
		// then mark it visited and enqueue it
		shuffle(&G->neighb[s][0], &G->neighb[s][G->DT[s]], default_random_engine());
		for(int ii=0; ii<G->DT[s]; ii++){
			int i = G->neighb[s][ii];
			if (!visited[i] && statusDFS[i]==-1)
			{
				visited[i] = true;
				parent[i] = s;
				*dim=*dim+1;

				*expr += cap*x[s];
				*expr += cap*x[i] + (1.0 - cap)*x[i];
				//as soon as the capacity is violated, store the expr to be lifted and added later to the model
				if(*dim > cap){
					int index =static_cast<int>(*dim - cap - 1);
					capConstrArray[index] = x[0] - x[0];
					capConstrArray[index] += *expr;
				}

				compVertices[nComp].push_back(i);
				statusDFS[i] = nComp;
				queue.push_back(i);
			}
		}
	}
}


// Clique algorithm for separation of clique in the capacity version
/***************************************************************************/
double findCapClique(int uIndex, graphFF G, IloNumArray xStar)
/***************************************************************************/
{
	
	int in=uIndex;
	int currMax= 0;
	int predMax=-1;
	double best=-1.0;
	double cliqueWeight=0.0;

	for(int i=0; i<vertex_number; i++){
		goodVertex[i]=0;
		goodCliqueU[i]=0;
	}
	
	double check;
	while(currMax>predMax){
		goodCliqueU[in]=1;
		cliqueWeight=cliqueWeight+(1-xStar[in]);
		predMax=currMax;

		for(int i=0; i<vertex_number; i++)
			goodVertex[i]=goodVertex[i]+G->AMatrix[in][i]+G->AMatrix[i][in];

		for(int i=0; i<vertex_number; i++){
			
			check = G->DT[i];
			
			if(goodVertex[i]>predMax && check>best && xStar[i]<0.9){				
				currMax = goodVertex[i];
				in = i;
				best = check;
			}
			//cout << best << endl;
		}
		best=-1.0;
	}
	return cliqueWeight;
}


// Returns true if there is a path from source 's' to sink 't' in residual graph. Also fills parent[] to store the path
/***************************************************************************/
bool bfs(int s, int t)
/***************************************************************************/
{

	//mark all vertices as not visited
	for(int i=0; i<vertex_number*2; i++)
		visited[i]=false;

	// Create a queue, enqueue source vertex and mark source vertex
	// as visited
	queue <int> q;
	q.push(s);
	visited[s] = true;
	parent[s] = -1;

	// Standard BFS Loop
	while (!q.empty())
	{
		int u = q.front();
		q.pop();

		for (int v=0; v<vertex_number*2; v++)
		{
			if (!visited[v] && rGraph[u][v]>0)
			{
				q.push(v);
				parent[v] = u;
				visited[v] = true;
			}
		}
	}

	// If we reached sink in BFS starting from source, then return
	// true, else false
	return visited[t];
}

// Returns the maximum flow from s to t in the given graph
/**************************************************************************/
int fordFulkerson(int s, int t, int dim, int nComp, int boundK)
/**************************************************************************/
{
	// Residual graph where rGraph[i][j] indicates residual capacity of edge from i to j (if there
	// is an edge. If rGraph[i][j] is 0, then there is not)
	for(int i=0; i<vertex_number*2; i++){
		parent[i]=-1;
		for(int j=0; j<vertex_number*2; j++){
			rGraph[i][j]=0;
		}
	}

	for(int ii=0; ii<dim; ii++){
		int i = compVertices[nComp][ii];

		if(i != s && i != t)
			rGraph[i+vertex_number][i]=1;

		for(int jj=ii+1; jj<dim; jj++){
			int j=compVertices[nComp][jj];

			if(G->AMatrix[i][j] == 1 || G->AMatrix[j][i] == 1){
				if(i == s)
					rGraph[i][j+vertex_number]=1;
				else if(i == t)
					rGraph[j][i]=1;
				else if(j == s)
					rGraph[j][i+vertex_number]=1;
				else if(j == t)
					rGraph[i][j]=1;
				else{
					rGraph[i][j+vertex_number]=1;
					rGraph[j][i+vertex_number]=1;
				}
			}
		}
	}
	int max_flow = 0;  // There is no flow initially
	// Augment the flow while there is path from source to sink
	while (bfs(s, t) && max_flow<boundK){
		// Find minimum residual capacity of the edges along the
		// path filled by BFS. Or we can say find the maximum flow
		// through the path found.
		int path_flow = INT_MAX;
		for (int v=t; v!=s; v=parent[v]){
			int u = parent[v];
			path_flow = min(path_flow, rGraph[u][v]);
		}

		// update residual capacities of the edges and reverse edges
		// along the path
		for (int v=t; v != s; v=parent[v]){
			int u = parent[v];
			rGraph[u][v] -= path_flow;
			rGraph[v][u] += path_flow;
		}

		// Add path flow to overall flow
		max_flow += path_flow;
	}

	// Return the overall flow
	return max_flow;
}

/***************************************************************************/
int kVertexConnect(int dim, int nComp, int boundK)
/***************************************************************************/
{
	int kVertex;
	kVertex=max(dim-1,1);
	int ii=0;
	while(ii<dim && kVertex>1){
		int i = compVertices[nComp][ii];
		int jj=ii+1;
		while(jj<dim && kVertex>1){
			int j = compVertices[nComp][jj];
			if(G->AMatrix[i][j] != 1 && G->AMatrix[j][i] != 1){
				int check=fordFulkerson(i, j, dim, nComp, boundK);
				kVertex=min(kVertex, check);
			}
			jj++;
		}
		ii++;
	}
	return kVertex;
}

//Simple DFS that enumerate the component and find its dimension
/***************************************************************************/
void simpleDFS(int v, int nComp, int *dim)
/***************************************************************************/
{
	*dim=*dim+1;
	compVertices[nComp].push_back(v);
	statusDFS[v] = nComp;

	for(int i=0; i<G->DT[v]; i++){
		int vNext=G->neighb[v][i];
		if(statusDFS[vNext] == -1){
			simpleDFS(vNext, nComp, &*dim);
		}
	}
}

//DFS to find and define a tree inside a component
/***************************************************************************/
void componentDFS(int v, int nComp, int *dim)
/***************************************************************************/
{
	*dim=*dim+1;
	compVertices[nComp].push_back(v);
	statusDFS[v] = nComp;

	shuffle(&G->neighb[v][0], &G->neighb[v][G->DT[v]], default_random_engine());
	for(int i=0; i<G->DT[v]; i++){
		int vNext=G->neighb[v][i];
		if(statusDFS[vNext] == -1){
			rGraph[v][vNext]=1;
			rGraph[vNext][v]=1;
			componentDFS(vNext, nComp, &*dim);
		}
	}
}

//DFS that moves only on the tree defined by the previous function (componentDFS)
/***************************************************************************/
void compTreeDFS(int v, int nComp, int *dim)
/***************************************************************************/
{
	*dim=*dim+1;
	statusDFS[v] = -1;

	for(int i=0; i<G->DT[v]; i++){
		int vNext=G->neighb[v][i];
		if(statusDFS[vNext] == nComp && rGraph[v][vNext])
			compTreeDFS(vNext, nComp, &*dim);
	}
}

/***************************************************************************/
IloInt findCompCoeff(int dim, int vIndexRemoved, int nComp)
/***************************************************************************/
{
	if(dim-cap==1)
		return 1;

	//See what happen if I remove vIndexComp vertex
	statusDFS[vIndexRemoved]=vertex_number+1;

	int newDim=0;
	int maxDim=0;
	//Run DFS on the component with vIndexComp removed.
	for(int i=0; i<dim; i++){
		int otherVertexIndex=compVertices[nComp][i];

		if(statusDFS[otherVertexIndex]==nComp){
			compTreeDFS(otherVertexIndex, nComp, &newDim);

			//When a new component is found, compare its dimension to see if it is the largest
			if(newDim>maxDim)
				maxDim=newDim;

			newDim=0;
		}
	}

	//Before end, restore the initial situation
	for(int i=0; i<dim; i++)
		statusDFS[compVertices[nComp][i]]=nComp;

	//Return the maximum dimension of components
	return maxDim;
}


//DFS that (try to) find a biconnected component and, as usual, define connected component
/***************************************************************************/
void connectivityDFS(int v, int nComp, int *dim, int *connValue, int pred, int d)
/***************************************************************************/
{
	depth[v]=d;
	low[v]=d;
	int children = 0;
	*dim=*dim+1;
	localDim[v]=*dim;
	statusDFS[v] = nComp;
	compVertices[nComp].push_back(v);

	for(int i=0; i<G->DT[v]; i++){
		int vNext=G->neighb[v][i];
		if(statusDFS[vNext] == -1){
			connectivityDFS(vNext, nComp, &*dim, &*connValue, v, d+1);
			children++;
			if((localDim[v]<=cap+2 && localDim[vNext]<=cap+2 && low[vNext] >= depth[v] && pred!=-1) || (pred==-1 && children>1))
				*connValue=1;
			if(localDim[v]<=cap+2 && localDim[vNext]<=cap+2)
				low[v]=min(low[v], low[vNext]);
		}
		else if(localDim[v]<=cap+2 && localDim[vNext]<=cap+2 && vNext != pred)
			low[v]=min(low[v], depth[vNext]);
	}
}

//DFS that (try to) find a biconnected component and, as usual, define connected component
/***************************************************************************/
void newConnectivityDFS(int v, int nComp, int *dim, int *connValue, int pred, int d, IloInt lazyLambda)
/***************************************************************************/
{
	depth[v]=d;
	low[v]=d;
	int children = 0;
	*dim=*dim+1;
	localDim[v]=*dim;
	statusDFS[v] = nComp;
	compVertices[nComp].push_back(v);

	for(int i=0; i<G->DT[v]; i++){
		int vNext=G->neighb[v][i];
		if(statusDFS[vNext] == -1){
			newConnectivityDFS(vNext, nComp, &*dim, &*connValue, v, d+1, lazyLambda);
			children++;
			if((localDim[v]<=lazyLambda+2 && localDim[vNext]<=lazyLambda+2 && low[vNext] >= depth[v] && pred!=-1) || (pred==-1 && children>1))
				*connValue=1;
			if(localDim[v]<=lazyLambda+2 && localDim[vNext]<=lazyLambda+2)
				low[v]=min(low[v], low[vNext]);
		}
		else if(localDim[v]<=lazyLambda+2 && localDim[vNext]<=lazyLambda+2 && vNext != pred)
			low[v]=min(low[v], depth[vNext]);
	}
}

//DFS that (try to) find a biconnected component and, as usual, define connected component
/***************************************************************************/
void mixedDFS(int v, int nComp, int *dim, int *connValue, int pred, int d)
/***************************************************************************/
{
	depth[v]=d;
	low[v]=d;
	int children = 0;
	*dim=*dim+1;
	localDim[v]=*dim;
	statusDFS[v] = nComp;
	compVertices[nComp].push_back(v);

	for(int i=0; i<G->DT[v]; i++){
		int vNext=G->neighb[v][i];
		if(statusDFS[vNext] == -1){
			//rGraph[v][vNext]=1;
			//rGraph[vNext][v]=1;
			mixedDFS(vNext, nComp, &*dim, &*connValue, v, d+1);
			children++;
			if((localDim[v]<=cap+2 && localDim[vNext]<=cap+2 && low[vNext] >= depth[v] && pred!=-1) || (pred==-1 && children>1))
				*connValue=1;
			if(localDim[v]<=cap+2 && localDim[vNext]<=cap+2)
				low[v]=min(low[v], low[vNext]);
		}
		else if(localDim[v]<=cap+2 && localDim[vNext]<=cap+2 && vNext != pred)
			low[v]=min(low[v], depth[vNext]);
	}
}

//Benders BFS that enumerate the component
/***************************************************************************/
void bendersBFS(int v, int nComp, int *dim)
/***************************************************************************/
{
	// Create a queue for BFS
	list<int> queue;

	// Mark the current node as visited and enqueue it
	visited[v] = true;
	*dim=*dim+1;
	compVertices[nComp].push_back(v);
	statusDFS[v] = nComp;
	queue.push_back(v);

	while(!queue.empty())
	{
		// Dequeue a vertex from queue and print it
		int s = queue.front();
		queue.pop_front();

		// Get all adjacent vertices of the dequeued
		// vertex s. If a adjacent has not been visited,
		// then mark it visited and enqueue it
		shuffle(&G->neighb[s][0], &G->neighb[s][G->DT[s]], default_random_engine());
		for(int ii=0; ii<G->DT[s]; ii++){
			int i = G->neighb[s][ii];
			if (!visited[i] && statusDFS[i]==-1)
			{
				visited[i] = true;
				parent[i] = s;
				*dim=*dim+1;
				compVertices[nComp].push_back(i);
				statusDFS[i] = nComp;
				queue.push_back(i);
			}
		}
	}
}

//Component BFS that enumerate the component and find and define a tree inside a component
/***************************************************************************/
void componentBFS(int v, int nComp, int *dim)
/***************************************************************************/
{
	// Create a queue for BFS
	list<int> queue;

	// Mark the current node as visited and enqueue it
	visited[v] = true;
	*dim=*dim+1;
	compVertices[nComp].push_back(v);
	statusDFS[v] = nComp;
	queue.push_back(v);

	while(!queue.empty())
	{
		// Dequeue a vertex from queue and print it
		int s = queue.front();
		queue.pop_front();

		// Get all adjacent vertices of the dequeued
		// vertex s. If a adjacent has not been visited,
		// then mark it visited and enqueue it
		shuffle(&G->neighb[s][0], &G->neighb[s][G->DT[s]], default_random_engine());
		for(int ii=0; ii<G->DT[s]; ii++){
			int i = G->neighb[s][ii];
			if (!visited[i] && statusDFS[i]==-1)
			{
				visited[i] = true;
				parent[i] = s;
				*dim=*dim+1;
				compVertices[nComp].push_back(i);
				statusDFS[i] = nComp;
				queue.push_back(i);
				rGraph[s][i]=1;
				rGraph[i][s]=1;
			}
		}
	}
}

//Mixed BFS find and define a tree inside a (pre-defined) component
/***************************************************************************/
void mixedBFS(int v, int nComp)
/***************************************************************************/
{
	// Create a queue for BFS
	list<int> queue;

	// Mark the current node as visited and enqueue it
	visited[v] = true;
	queue.push_back(v);

	while(!queue.empty())
	{
		// Dequeue a vertex from queue and print it
		int s = queue.front();
		queue.pop_front();

		// Get all adjacent vertices of the dequeued
		// vertex s. If a adjacent has not been visited,
		// then mark it visited and enqueue it
		shuffle(&G->neighb[s][0], &G->neighb[s][G->DT[s]], default_random_engine());
		for(int ii=0; ii<G->DT[s]; ii++){
			int i = G->neighb[s][ii];
			if (!visited[i] && statusDFS[i]==nComp)
			{
				visited[i] = true;
				parent[i] = s;
				queue.push_back(i);
				rGraph[s][i]=1;
				rGraph[i][s]=1;
			}
		}
	}
}
