//
// Created by paolo on 15/11/18.
//

#include "ConnectivityCut.h"
#include <list>
#include "binPackingProblem.h"
#include "Graph_v4.h"

ILOLAZYCONSTRAINTCALLBACK1(LazyConnectivityCallback, IloIntVarArray, x)
{
	lazyCall++;
	IloEnv env = getEnv();
	getValues(lazySolx,x);
	int nComp, dim, randomSeed;
	bool lazyFound=false;
	nComp=0;
	dim=0;
	for(int i=0; i<vertex_number; i++){
		vindex[i]=i;
		if(lazySolx[i] <= 0.1)
			statusDFS[i] = -1;
		else
			statusDFS[i] = vertex_number+1;
		low[i] = vertex_number+1;
		depth[i] = vertex_number+1;
		localDim[i] = vertex_number+1;
	}

	randomSeed = rand() % 500 + 1;
	shuffle(vindex, static_cast<size_t>(vertex_number), randomSeed);

	compVertices.resize(static_cast<unsigned long>(vertex_number));
	weights.clear();
	////////////////////Connectivity Cuts//////////////
	IloExpr connConstr(env);
	for(int i=0; i<vertex_number; i++){
		int v=vindex[i];
		if(statusDFS[v] == -1){
			compVertices[nComp].reserve(static_cast<unsigned long>(vertex_number));
			int connValue=2;
			connectivityDFS(v, nComp, &dim, &connValue, -1, 0);
			weights.push_back(dim);
			int diff=static_cast<int>(dim - cap);
			connValue=min(connValue, diff);
			if(dim>cap){
				lazyFound=true;
				//At least connValue vertices in the sub-component with cardinality cap+connValue
				for(int j=0; j<cap+connValue; j++){
					int vIndexComp=compVertices[nComp][j];
					connConstr+=x[vIndexComp];
				}
				add(connConstr>=connValue);
				kConnCuts++;
				connConstr.clear();
			}
			nComp++;
			dim=0;
		}
	}
	connConstr.end();
	////////////////////////////////////////////////////////////////////////////////

	if(lazyFound){
		weights.clear();
		compVertices.clear();
		return;
	}

	bool checkBPP=true;
	if(nComp>0)
		checkBPP=lazyBinPackingModel(env, nComp);

	if(!checkBPP){
		IloExpr bppConstr(env);
		bppConstr=x[0]-x[0];
		for(int i=0; i<vertex_number; i++){
			if(statusDFS[i]<vertex_number+1)
				bppConstr += x[i];
		}
		add(bppConstr>=1);
		bppCuts++;
		bppConstr.end();
	}
	compVertices.clear();
	weights.clear();
}

//DFS that (try to) find a biconnected component and, as usual, define connected component
/***************************************************************************/
void userConnectivityDFS(int v, int nComp, int *dim, double *size, int *connValue, int pred, int d, IloExpr *expr, IloIntVarArray x)
/***************************************************************************/
{
	depth[v]=d;
	low[v]=d;
	int children = 0;
	*dim=*dim+1;

	if(*dim<=cap+*connValue){
		*size=*size + userSolx[v];
		*expr+=x[v];
	}

	localDim[v]=*dim;
	statusDFS[v] = nComp;


	for(int i=0; i<vertex_number; i++){
		int vNext=vindex[i];
		if(statusDFS[vNext] == -1 && vNext!=v && symmAdjMatrix[v][vNext]){
			userConnectivityDFS(vNext, nComp, &*dim, &*size, &*connValue, v, d+1, expr, x);
			children++;
			if((localDim[v]<=cap+2 && localDim[vNext]<=cap+2 && low[vNext] >= depth[v] && pred!=-1) || (pred==-1 && children>1))
				*connValue=1;
			if(localDim[v]<=cap+2 && localDim[vNext]<=cap+2)
				low[v]=min(low[v], low[vNext]);
		}
		else if(localDim[v]<=cap+2 && localDim[vNext]<=cap+2 && vNext != pred && vNext!=v && symmAdjMatrix[v][vNext])
			low[v]=min(low[v], depth[vNext]);
	}
	if(pred==-1 && children>1)
		*connValue=1;
}

ILOUSERCUTCALLBACK1(UserConnectivityCallback, IloIntVarArray, x){
	IloInt nNodes = getNnodes();
	//if(nNodes % userFreq != 0) return;
	//if(nNodes > 0) return;
	getValues(userSolx,x);
	IloEnv env = getEnv();

	int nComp, dim;
	double size;
	nComp=0;
	dim=0;
	size=0;
	for(int i=0; i<vertex_number; i++){
		vindex[i]=i;
		vweight[i] = userSolx[i];
	}
	SORT_NON_DECR(vindex,vweight,vertex_number);

	////////////////////Connectivity Cuts//////////////
	IloExpr connConstr(env);
	int i=0;

	while(i<vertex_number && userSolx[vindex[i]]<0.5){
		int v=vindex[i];

		for(int j=0; j<vertex_number; j++){
			if(userSolx[j] < 0.9)
				statusDFS[j] = -1;
			else
				statusDFS[j] = vertex_number+1;
			low[j] = -1;
			depth[j] = -1;
			localDim[j] = vertex_number;
		}

		int connValue=2;
		connConstr=x[0]-x[0];
		userConnectivityDFS(v, nComp, &dim, &size, &connValue, -1, 0, &connConstr, x);
		int diff=static_cast<int>(dim - cap);
		connValue=min(connValue, diff);
		if(size<connValue-tolerance){
			//At least connValue vertices in the sub-component with cardinality cap+connValue
			add(connConstr>=connValue);
			userCuts++;
			connConstr.end();
			return;
		}
		connConstr.clear();
		nComp++;
		dim=0;
		size=0;
		i++;
	}
	connConstr.end();
}

void ConnectivityCutModel()
{

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	parent = new int[vertex_number*2];
	statusDFS = new int[vertex_number];
	wweight = new double[edge_number];
	iindex = new int[edge_number];
	vindex = new int[vertex_number];
	vweight = new double[vertex_number];
	visited = new bool[vertex_number*2]();
	depth = new int[vertex_number];
	low = new int[vertex_number];
	localDim = new int[vertex_number];

	rGraph = new int*[vertex_number*2];
	for(int i=0; i<vertex_number*2; i++)
		rGraph[i] = new int[vertex_number*2]();

	///////////////////////////////////////////////////////////////////
	//Create a symmetric adjacent matrix
	symmAdjMatrix = new int*[vertex_number];
	for(int i=0; i<vertex_number; i++)
		symmAdjMatrix[i]=new int[vertex_number];

	for(int i=0; i<vertex_number; i++)
		for(int j=0; j<vertex_number; j++)
			symmAdjMatrix[i][j] = G->AMatrix[i][j] + G->AMatrix[j][i];
	///////////////////////////////////////////////////////////////////

	try{

		IloIntVarArray x(env, vertex_number, 0, 1);
		userSolx = IloNumArray(env, vertex_number);
		lazySolx = IloNumArray(env, vertex_number);
		degreeTree = IloIntArray(env, vertex_number);
		capConstrArray = IloExprArray(env, static_cast<IloInt>(vertex_number - cap));

		//Objective Function: minimize the number of vertices in the separator
		IloObjective obj(env, IloSum(x), IloObjective::Minimize);
		model.add(obj);

		//Impose priority between vertices
		for(int i=0; i<vertex_number; i++){
			for(int j=i+1; j<vertex_number; j++){
				int a=0;
				bool go=true;
				while(a<vertex_number && go){
					if(a != i && a != j){
						go = (symmAdjMatrix[i][a] >= symmAdjMatrix[j][a]);
					}
					a++;
				}
				if(go){
					model.add(x[i]>=x[j]);
				}
			}
		}
		for(int i=vertex_number-1; i>=0; i--){
			for(int j=i-1; j>=0; j--){
				int a=0;
				bool go=true;
				while(a<vertex_number && go){
					if(a != i && a!= j){
						go = (symmAdjMatrix[i][a] >= symmAdjMatrix[j][a]);
					}
					a++;
				}
				if(go && (G->DT[i] > G->DT[j])){
					model.add(x[i]>=x[j]);
				}
			}
		}

		//New Degree Constraint: sum_(v in V : v in Adj(i)) x_v >= (deg(i) +1 -u)*(1-x)
		IloExpr degreeConstr(env);
		for(int i=0; i<vertex_number; i++){
			if(G->DT[i] >= cap){
				degreeConstr = (G->DT[i] +1 -cap)*(x[i]-1);
				for(int j=0; j<vertex_number; j++){
					if(G->AMatrix[i][j] == 1 || G->AMatrix[j][i] == 1)
						degreeConstr += x[j];
				}
				model.add(degreeConstr >= 0.0);
				degreeConstr.clear();
			}
		}
		degreeConstr.end();

		IloCplex::Callback sec;
		sec = cplex.use(LazyConnectivityCallback(env, x));

		IloCplex::Callback cut;
		//cut = cplex.use(UserConnectivityCallback(env, x));

		IloCplex::Callback inc;
		//inc = cplex.use(IncumbentCallback(env, x));

		//cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, treeLimit);
		//cplex.setParam(IloCplex::Param::Emphasis::MIP, 1);
		//cplex.setParam(IloCplex::Param::MIP::Strategy::FPHeur, 1);
		//cplex.setParam(IloCplex::Param::MIP::Strategy::Dive, 2);
		//cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
		//cplex.setParam(IloCplex::NodeLim, 0);
		//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
		//cplex.setParam(IloCplex::PreInd, IloFalse);

		clock_t time1 = clock(); //cplex.getTime();
		cplex.solve();
		clock_t time2 = clock(); //cplex.getTime();
		double time = (double) (time2 - time1)/CLOCKS_PER_SEC;

		//cplex.exportModel("DegreeCutModel.lp");

		/*for(int i=0; i<vertex_number; i++)
			if(cplex.getValue(x[i])>=0.8)
				cout << i << ", ";
		cout << endl;*/

		// Check if the corresponding instance of Bin Packing Problem is feasible. If yes, then the solution is optimal
		//bool bbpSolved=binPackingModel(G, x, env, cplex);

		if(cplex.getStatus() == IloAlgorithm::Unknown || cplex.getStatus() == IloAlgorithm::Infeasible){
			IloNum bound = cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();
			IloInt status = cplex.getStatus();
			IloInt nVar = cplex.getNcols();
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << "null" << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t"
					   << nConst << "\t" << nNode << "\t" << userCuts << "\t" << compCuts << "\t" << kConnCuts << "\t" <<  bppCuts << "\t" << lazyCall;
			}
		}
		else{
			IloNum objective = cplex.getObjValue();
			IloNum bound = cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();
			IloInt status = cplex.getStatus();
			IloInt nVar = cplex.getNcols()-1;
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << objective << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t"
					   << nConst << "\t" << nNode << "\t" << userCuts << "\t" << compCuts << "\t" << kConnCuts << "\t" << bppCuts << "\t" << lazyCall;
			}
		}

		//cout << cplex.getStatus() << endl;
		sec.end();
		cut.end();
		obj.end();
		x.end();
		userSolx.end();
		lazySolx.end();
		degreeTree.end();
		capConstrArray.end();
	}

	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;

		IloInt nNode = cplex.getNnodes();
		IloNum bound = cplex.getBestObjValue();
		IloInt status = cplex.getStatus();
		IloInt nVar = cplex.getNcols();
		IloInt nConst = cplex.getNrows();

		if(output.is_open())
		{
			output.precision(10);
			output << "null" << "\t" << bound << "\t"  << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t"
				   << nConst << "\t" << nNode << "\t" << userCuts << "\t" << compCuts << "\t" << kConnCuts << "\t" << bppCuts << "\t" << lazyCall;
		}
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	cplex.end();
	model.end();
	env.end();

	for(int i=0; i<vertex_number*2; i++)
		delete [] rGraph[i];
	delete [] rGraph;

	for(int i=0; i<vertex_number; i++)
		delete [] symmAdjMatrix[i];
	delete [] symmAdjMatrix;

	delete [] parent;
	delete [] wweight;
	delete [] vweight;
	delete [] iindex;
	delete [] vindex;
	delete [] statusDFS;
	delete [] visited;
	delete [] depth;
	delete [] low;
	delete [] localDim;
	compVertices.end();
	weights.end();
}


