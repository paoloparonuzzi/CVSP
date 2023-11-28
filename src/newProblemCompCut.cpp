//
// Created by paolo on 27/02/20.
//

#include "newProblemCompCut.h"

ILOLAZYCONSTRAINTCALLBACK2(LazyNewProblemCompCallback, IloIntVarArray, x, IloNumVar, lambda)
{
	lazyCall++;
	IloEnv env = getEnv();
	getValues(lazySolx,x);
	IloNum localLambda=getValue(lambda);

	int nComp, dim, randomSeed;
	//////////////////////////////Component Cuts////////////////////////////////

	compVertices.resize(static_cast<unsigned long>(vertex_number));
	nComp=0;
	dim=0;
	for(int i=0; i<vertex_number; i++){
		vindex[i]=i;
		visited[i] = false;
		if(lazySolx[i] <= 0.5)
			statusDFS[i] = -1;
		else
			statusDFS[i] = vertex_number+1;
	}

	randomSeed = rand() % 500 + 1;
	shuffle(vindex, static_cast<size_t>(vertex_number), randomSeed);

	//Used to store the trees
	for(int i=0; i<vertex_number; i++)
		for(int j=0; j<vertex_number; j++)
			rGraph[i][j]=0;

	IloExpr compConstr(env);
	for(int i=0; i<vertex_number; i++){
		int v=vindex[i];
		if(statusDFS[v] == -1){
			compVertices[nComp].reserve(static_cast<unsigned long>(vertex_number));
			//componentDFS(v, nComp, &dim);
			componentBFS(v, nComp, &dim);
			if(dim>localLambda){
				compConstr=x[0] - x[0];
				for(int j=0; j<dim; j++){
					int vIndexComp=compVertices[nComp][j];
					IloInt maxDim = findCompCoeff(dim, vIndexComp, nComp);
					IloInt compCoeff = dim - maxDim;
					compConstr+=compCoeff*x[vIndexComp];
				}
				add(compConstr>=dim - lambda);
				compCuts++;
				compConstr.clear();
			}
			nComp++;
			dim=0;
		}
	}
	compConstr.end();
	compVertices.clear();
}

void minMaxC_CompCutModel()
{

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	parent = new int[vertex_number];
	statusDFS = new int[vertex_number];
	vindex = new int[vertex_number];
	vweight = new double[vertex_number];
	visited = new bool[vertex_number]();

	rGraph = new int*[vertex_number];
	for(int i=0; i<vertex_number; i++)
		rGraph[i] = new int[vertex_number]();

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
		IloNumVar lambda(env, 0, vertex_number-k);
		userSolx = IloNumArray(env, vertex_number);
		lazySolx = IloNumArray(env, vertex_number);
		degreeTree = IloIntArray(env, vertex_number);
		capConstrArray = IloExprArray(env, static_cast<IloInt>(vertex_number - cap));

		//Objective Function: minimize the size of the largest component (lambda)
		IloObjective obj(env, lambda, IloObjective::Minimize);
		model.add(obj);

		//No more than k vertices ca be removed
		model.add(IloSum(x) <= k);

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

		IloCplex::Callback sec;
		sec = cplex.use(LazyNewProblemCompCallback(env, x, lambda));

		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, treeLimit);
		//////////////////////////////////////////////////////////////
        cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 0.9);

		clock_t time1 = clock();
		IloBool solved = cplex.solve();
		clock_t time2 = clock();
		double time = (double) (time2 - time1)/CLOCKS_PER_SEC;

		if(solved){
			cout << "Remove: ";
		for(int i=0; i<vertex_number; i++)
			if(cplex.getValue(x[i])>=0.8)
				cout << i << " ";
		cout << endl;
		cout << "Max size: " << cplex.getValue(lambda) << endl;
		}
		
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
					   << nConst << "\t" << nNode << "\t" << userCuts  << "\t" << compCuts << "\t" << kConnCuts << "\t" <<  bppCuts << "\t" << lazyCall;
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
					   << nConst << "\t" << nNode << "\t" << userCuts  << "\t" << compCuts << "\t" << kConnCuts << "\t" << bppCuts << "\t" << lazyCall;
			}
		}

		//cout << cplex.getStatus() << endl;
		sec.end();
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

	for(int i=0; i<vertex_number; i++)
		delete [] rGraph[i];
	delete [] rGraph;

	for(int i=0; i<vertex_number; i++)
		delete [] symmAdjMatrix[i];
	delete [] symmAdjMatrix;

	delete [] parent;
	delete [] vweight;
	delete [] vindex;
	delete [] statusDFS;
	delete [] visited;

	compVertices.end();
	weights.end();
}

