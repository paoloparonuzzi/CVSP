//
// Created by paolo on 22/11/18.
//

#include "MixedModel.h"
#include <list>
#include "binPackingProblem.h"
#include "Graph_v4.h"

double bppTime=0.0;

//Callback to run the heuristic in (each) node to find a feasible solution and maybe help cplex
ILOHEURISTICCALLBACK1(HeuristicCallback, IloIntVarArray, x){
	
	IloInt nNodes = getNnodes();
	if(nNodes % userFreq != 0) return;
	
	IloEnv env = getEnv();
	getValues(userSolx,x);
	
	double objHeur=0;
	bool feasible=false;
	IloNumArray scoreHeur(env, vertex_number);
	
	for(int i=0; i<vertex_number; i++){
		removed[i]=false;
		scoreHeur[i]=(1-userSolx[i])*G->DT[i];
	}
	
	while(!feasible){
		compVertices.resize(static_cast<unsigned long>(vertex_number));
		int nComp=0;
		int dim=0;
		feasible=true;

		for(int i=0; i<vertex_number; i++){
			if(removed[i])	statusDFS[i]=vertex_number+1;
			else			statusDFS[i]=-1;
		}

		for(int v=0; v<vertex_number; v++){
			if(statusDFS[v] == -1){
				compVertices[nComp].reserve(static_cast<unsigned long>(vertex_number));
				simpleDFS(v, nComp, &dim);
				if(dim>cap){
					feasible=false;
					int best=compVertices[nComp][0];
					for(int u=1; u<dim; u++){
						int check=compVertices[nComp][u];
						if(scoreHeur[check] > scoreHeur[best])
							best=check;
					}
					for(int n=0; n<G->DT[best]; n++)
						scoreHeur[G->neighb[best][n]]=scoreHeur[G->neighb[best][n]]-(1-userSolx[G->neighb[best][n]]);					
					removed[best]=true;
					objHeur++;
				}
				nComp++;
				dim=0;
			}
		}
		compVertices.clear();
	}

	setSolution(x, scoreHeur, objHeur);
	compVertices.clear();
	scoreHeur.end();	
}

ILOLAZYCONSTRAINTCALLBACK1(LazyMixedCallback, IloIntVarArray, x)
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
		low[i] = -1;
		depth[i] = -1;
		localDim[i] = vertex_number;
	}

	//Used to store the trees for the Component Cuts
	for(int i=0; i<vertex_number; i++)
		for(int j=0; j<vertex_number; j++)
			rGraph[i][j]=0;

	randomSeed = rand() % 500 + 1;
	shuffle(vindex, static_cast<size_t>(vertex_number), randomSeed);

	compVertices.resize(static_cast<unsigned long>(vertex_number));
	weights.clear();
	vector<int> connectivityValues;
	connectivityValues.reserve(static_cast<unsigned long>(vertex_number));
	////////////////////Connectivity Cuts//////////////
	IloExpr connConstr(env);
	for(int i=0; i<vertex_number; i++){
		int v=vindex[i];
		if(statusDFS[v] == -1){
			compVertices[nComp].reserve(static_cast<unsigned long>(vertex_number));
			int connValue=2;
			mixedDFS(v, nComp, &dim, &connValue, -1, 0);
			weights.push_back(dim);
			connectivityValues.push_back(connValue);
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

	//////////////////////////////Component Cuts////////////////////////////////
	IloExpr compConstr(env);
	for(int i=0; i<nComp; i++){
		dim=weights[i];
		if((dim==cap+2 && connectivityValues[i]==1) || (dim>cap+2)){
			//Component Cut
			compConstr=x[0] - x[0];
			for(int j=0; j<dim; j++){
				int vIndexComp=compVertices[i][j];
				IloInt maxDim = findCompCoeff(dim, vIndexComp, nComp);
				IloInt compCoeff=(IloInt) min((double)dim - cap, (double)dim - maxDim);
				compConstr+=compCoeff*x[vIndexComp];
			}
			add(compConstr>=dim - cap);
			compCuts++;
			compConstr.clear();
		}
	}
	compConstr.end();

	/////////////////////////////////////////////////////////////////////////////

	
	if(lazyFound){
		weights.clear();
		compVertices.clear();
		return;
	}

	//BBP constraints
	bool checkBPP=true;
	if(nComp>0){
		double before=clock();
		checkBPP=lazyBinPackingModel(env, nComp);
		double after=clock();
		bppTime+= (after-before);
	}
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

ILOUSERCUTCALLBACK1(UserMixedCallback, IloIntVarArray, x){

	//IloInt nNodes = getNnodes();
	//if(nNodes % userFreq != 0) return;
	//if(nNodes > 0) return;
	getValues(userSolx,x);
	IloEnv env = getEnv();

	for(int i=0; i<vertex_number; i++){
		vindex[i]=i;
		vweight[i]=userSolx[i];
	}

	SORT_NON_DECR(vindex,vweight,vertex_number);

	IloExpr connConstr(env);
	vector<int> tree;
	tree.resize(vertex_number);

	int i=0;
	while(i<vertex_number && userSolx[vindex[i]]<1){
		for(int j=0; j<vertex_number; j++)
			visited[j]=false;

		int v=vindex[i];
		connConstr=x[v];
		int conn=2;
		visited[v]=true;
		double size=userSolx[v];
		int dim=1;
		tree.push_back(v);

		int j=i+1;
		while(j<vertex_number && userSolx[vindex[j]]<1){
			int u=vindex[j];

			if(visited[u] || u==v){
				j++;
				continue;
			}

			int neighbors=0;
			int n=0;
			while(neighbors<=1 && n<G->DT[u]){
				if(visited[G->neighb[u][n]])
					neighbors++;
				n++;
			}

			if(neighbors==0){
				j++;
				continue;
			}
			if(dim==cap+1 && neighbors==1){
				j++;
				continue;
			}
			if(dim >= 2 && neighbors==1) conn=1;

			visited[u]=true;
			tree.push_back(u);
			size+=userSolx[u];
			connConstr+=x[u];
			dim++;
			if(dim>=cap+conn){
				if(size<conn-tolerance){
					add(connConstr>=conn);
					//cout << connConstr << " >= " << conn << endl;
					userCuts++;
				}
				connConstr.end();
				return;
			}
			j=i+1;
		}
		tree.clear();
		i++;
	}

}

void MixedModel()
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
		userSolx = IloNumArray(env, vertex_number);
		lazySolx = IloNumArray(env, vertex_number);
		degreeTree = IloIntArray(env, vertex_number);
		capConstrArray = IloExprArray(env, static_cast<IloInt>(vertex_number - cap));

		IloIntVarArray x(env, vertex_number, 0, 1);

		//Objective Function: minimize the number of vertices in the separator
		IloObjective obj(env, IloSum(x), IloObjective::Minimize);
		model.add(obj);

		if(option>0){
			IloNumArray startVal(env);
			IloNumVarArray startVar(env);
			for(int i=0; i<vertex_number; i++){
				startVar.add(x[i]);
				startVal.add(removed[i]);
			}
			cplex.addMIPStart(startVar, startVal);
			startVal.end();
			startVar.end();

			//Save Leaves and isolated vertices
			/*for(int i=0; i<vertex_number; i++)
				if(G->DT[i]<=1) x[i].setUB(0);
			*/
		
		}

		//Priority constraints between equivalent and dominated vertices
		
		for(int i=0; i<vertex_number; i++){
			for(int j=i + 1; j<vertex_number; j++){
				int a=0;
				bool go=true;
				while(a<vertex_number && go){
					if(a != i && a != j){
						go=(symmAdjMatrix[i][a]>=symmAdjMatrix[j][a]);
					}
					a++;
				}
				if(go){
					model.add(x[i]>=x[j]);
				}
			}
		}
		for(int i=vertex_number - 1; i>=0; i--){
			for(int j=i - 1; j>=0; j--){
				int a=0;
				bool go=true;
				while(a<vertex_number && go){
					if(a != i && a != j){
						go=(symmAdjMatrix[i][a]>=symmAdjMatrix[j][a]);
					}
					a++;
				}
				if(go && (G->DT[i]>G->DT[j])){
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
		sec = cplex.use(LazyMixedCallback(env, x));

		IloCplex::Callback cut;
		if(userFreq>0)
			cut = cplex.use(UserMixedCallback(env, x));
		
		IloCplex::Callback heur;
		//heur = cplex.use(HeuristicCallback(env,x));

		cplex.setOut(env.getNullStream());
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

		/*for(int i=0; i<vertex_number; i++)
			if(cplex.getValue(x[i])>=0.8)
				cout << i << ", ";
		cout << endl;*/

		if(cplex.getStatus() == IloAlgorithm::Unknown || cplex.getStatus() == IloAlgorithm::Infeasible){
			IloNum bound = cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();
			IloInt status = cplex.getStatus();
			IloInt nVar = cplex.getNcols();
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << "null" << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t"  << nConst << "\t" << nNode
					   << "\t" << userCuts << "\t" << compCuts << "\t" << kConnCuts << "\t" <<  bppCuts << "\t" << bppTime/CLOCKS_PER_SEC  << "\t" << lazyCall;
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
				output << objective << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" << nConst << "\t" << nNode
					   << "\t" << userCuts << "\t" << compCuts << "\t" << kConnCuts << "\t" << bppCuts << "\t" << bppTime/CLOCKS_PER_SEC << "\t" << lazyCall;
			}
		}

		//cout << cplex.getStatus() << endl;
		sec.end();
		cut.end();
		heur.end();
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
			output << "null" << "\t" << bound << "\t"  << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t"  << nConst << "\t" << nNode
				   << "\t" << userCuts << "\t" << compCuts << "\t" << kConnCuts << "\t" << bppCuts << "\t" << bppTime/CLOCKS_PER_SEC  << "\t" << lazyCall;
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

	if(option>0)
		delete [] removed;
	
	compVertices.end();
	weights.end();
}


