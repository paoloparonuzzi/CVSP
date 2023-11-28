
#include "DegreeCut.h"


ILOLAZYCONSTRAINTCALLBACK1(LazyDegreeCallback, IloIntVarArray, x)
//ILOINCUMBENTCALLBACK1(IncumbentCallback, IloIntVarArray, x)
{
	lazyCall++;
	IloEnv env = getEnv();
	getValues(lazySolx,x);
	int nComp, dim, randomSeed;
	bool lazyFound=false;

//////////////////////////////Degree Cuts////////////////////////////////
	compVertices.resize(static_cast<unsigned long>(vertex_number));
	weights.clear();
	nComp=0;
	dim=0;
	//clock_t time11 = clock();
	for(int i=0; i<vertex_number; i++){
		vindex[i]=i;
		parent[i]=-1;
		visited[i] = false;
		if(lazySolx[i] <= 0.5)
			statusDFS[i] = -1;
		else
			statusDFS[i] = vertex_number+1;
	}


	randomSeed = rand() % 1000 + 1;
	shuffle(vindex, static_cast<size_t>(vertex_number), randomSeed);

	IloNum coef, rhs;
	IloNumVar xxx;
	IloExpr capConstr(env);
	IloExpr::LinearIterator it;

	for(int i=0; i<vertex_number; i++){
		int v=vindex[i];
		if(statusDFS[v] == -1){

			capConstr = -cap*x[v]; //I must inizialize with "minus first vertex" to don't count it twice.

			//lazyArrayDFS(&capConstr, x, v, nComp, &dim);
			lazyArrayBFS(&capConstr, x, v, nComp, &dim);
			weights.push_back(dim);
			if(dim > cap){
				lazyFound=true;
				for(int j=0; j<(dim-cap); j++){
					capConstr = capConstrArray[j];
					rhs = (IloNum) (j+1);
					it = capConstr.getLinearIterator();
					while(it.ok()){
						coef = it.getCoef();
						xxx = it.getVar();
						capConstr.setLinearCoef(xxx, (IloNum) min(coef, rhs));
						++it;
					}
					add(capConstr >= rhs).end();
					//cout << "Cap Contr: " << capConstr << " >= " << rhs << endl;

					capConstr.clear();
					lazyCuts++;
					capConstrArray[j]=nullptr;
				}
			}
			nComp++;
			dim=0;
		}
	}
	capConstr.end();

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

/*
ILOUSERCUTCALLBACK1(UserCapSubgraphCallback, IloIntVarArray, x)
{
	IloInt nNodes = getNnodes();
	if(nNodes % userFreq != 0) return;

	IloEnv env = getEnv();
	getValues(userSolx,x);

	int iMaxWeight=0;
	for(int i=0;i<edge_number;i++){
		wweight[i]=1 -userSolx[G->H[i]] -userSolx[G->T[i]];
		if(wweight[i]>wweight[iMaxWeight])
			iMaxWeight=i;
		iindex[i]=i;
	}

	//SORT_NON_INCR(iindex,wweight,edge_number);

	int a, b, next=0, sizeTree;
	IloNum coef;
	IloNumVar xxx;
	IloExpr capConstr(env);
	IloExpr::LinearIterator it;
	IloNum rhs;

	for(int i=0; i<vertex_number; i++){
		visited[i]=false;
		degreeTree[i]=0;
	}
	visited[G->H[iMaxWeight]]=true;
	visited[G->T[iMaxWeight]]=true;
	degreeTree[G->T[iMaxWeight]]=1;
	degreeTree[G->H[iMaxWeight]]=1;

	double obj = wweight[iMaxWeight] -(1.0-userSolx[G->T[iMaxWeight]])*(1.0-1.0/cap) -(1.0-userSolx[G->H[iMaxWeight]])*(1.0-1.0/cap);
	sizeTree=2;

	int randomSeed = rand() % 1000 + 1;
	shuffle(iindex, static_cast<size_t>(edge_number), randomSeed);

	next=0;
	while(sizeTree<vertex_number-1 && next<edge_number)
	{

		a=G->T[iindex[next]];
		b=G->H[iindex[next]];

		if(visited[a] && !visited[b] && (wweight[iindex[next]]-(1-userSolx[b])*(1.0-1.0/cap))>0){
			visited[b]=true;
			degreeTree[b]=1;
			degreeTree[a]=static_cast<IloInt>(degreeTree[a] + cap);
			sizeTree++;
			obj=obj+wweight[iindex[next]]-(1.0-userSolx[b])*(1.0-1.0/cap);
			next=0;
		}
		else if(!visited[a] && visited[b] && (wweight[iindex[next]]-(1-userSolx[a])*(1.0-1.0/cap))>0){
			visited[a]=true;
			degreeTree[a]=1;
			degreeTree[b]=static_cast<IloInt>(degreeTree[b] + cap);
			sizeTree++;
			obj=obj+wweight[iindex[next]]-(1.0-userSolx[a])*(1.0-1.0/cap);
			next=0;
		}
		else	next++;
	}
	if(obj>tolerance){
		capConstr = IloScalProd(degreeTree, x);
		rhs = ((IloNum) sizeTree)- cap;
		it = capConstr.getLinearIterator();
		while(it.ok()){
			coef = it.getCoef();
			xxx = it.getVar();
			capConstr.setLinearCoef(xxx, min(coef, rhs));
			++it;
		}
		add(capConstr >= rhs).end();
		//cout << capConstr << " >= " << rhs << endl;
		userCuts++;
		capConstr.clear();
	}
	capConstr.end();
}


ILOUSERCUTCALLBACK1(UserDegreeCallback, IloIntVarArray, x){
	IloInt nNodes = getNnodes();
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

		int j=0;
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
					userCuts++;
				}
				connConstr.end();
				return;
			}
			j++;
		}
		tree.clear();
		i++;
	}
}
*/

void DegreeCutModel()
{

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	parent = new int[vertex_number];
	statusDFS = new int[vertex_number];
	vindex = new int[vertex_number];
	visited = new bool[vertex_number]();

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
		sec = cplex.use(LazyDegreeCallback(env, x));
		
		IloCplex::Callback cut;
		/*if(userFreq>0)
			cut = cplex.use(UserDegreeCallback(env, x));*/

		IloCplex::Callback inc;
		//inc = cplex.use(IncumbentCallback(env, x));

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
		
		//cplex.exportModel("DegreeCutModel.lp");

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
				output << "null" << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t"
				<< nConst << "\t" << nNode << "\t" << userCuts << "\t" << lazyCuts << "\t" << compCuts << "\t" << kConnCuts << "\t" <<  bppCuts << "\t" << lazyCall;
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
				<< nConst << "\t" << nNode << "\t" << userCuts << "\t" << lazyCuts << "\t" << compCuts << "\t" << kConnCuts << "\t" << bppCuts << "\t" << lazyCall;
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
				<< nConst << "\t" << nNode << "\t" << userCuts << "\t" << lazyCuts<< "\t" << compCuts << "\t" << kConnCuts << "\t" << bppCuts << "\t" << lazyCall;
		}
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	cplex.end();
	model.end();
	env.end();

	for(int i=0; i<vertex_number; i++)
		delete [] symmAdjMatrix[i];
	delete [] symmAdjMatrix;

	delete [] parent;
	delete [] vindex;
	delete [] statusDFS;
	delete [] visited;
	compVertices.end();
	weights.end();
}
