//
// Created by paolo on 09/07/19.
//

#include "BendersCut.h"
#include "Graph_v4.h"

ILOLAZYCONSTRAINTCALLBACK1(LazyBendersCallback, IloIntVarArray, x)
//ILOINCUMBENTCALLBACK1(IncumbentCallback, IloIntVarArray, x)
{
	lazyCall++;
	IloEnv env = getEnv();
	getValues(lazySolx,x);
	int nComp, dim, randomSeed;
	bool lazyFound=false;

//////////////////////////////Benders Cuts////////////////////////////////
	compVertices.resize(static_cast<unsigned long>(vertex_number));
	weights.clear();
	nComp=0;
	dim=0;
	//clock_t time11 = clock();
	for(int i=0; i<vertex_number; i++){
		bondersCoeff[i]=0;
		vindex[i]=i;
		vweight[i]=G->DT[i];
		parent[i]=-1;
		visited[i] = false;
	}

	//cout << "deleted: ";
	for(int i=0; i<vertex_number; i++){
		if(lazySolx[i] <= 0.5)
			statusDFS[i] = -1;
		else{
			//cout << i << " ";
			statusDFS[i]=vertex_number + 1;
			for(int j=0; j<G->DT[i]; j++)
				vweight[G->neighb[i][j]]--;
		}
	}
	//cout << endl;

	SORT_NON_INCR(vindex,vweight,vertex_number);

	IloExpr bendersConstr(env);
	for(int i=0; i<vertex_number; i++){
		int v=vindex[i];
		if(statusDFS[v] == -1){

			bendersBFS(v, nComp, &dim);
			weights.push_back(dim);
			if(dim > cap){
				lazyFound=true;

				for(int j=0; j<dim; j++){
					int vIndexComp=compVertices[nComp][j];
					bondersCoeff[vIndexComp]=1;
					int p = parent[vIndexComp];
					while(p != -1){
						bondersCoeff[p]=bondersCoeff[p]+2;
						p=parent[p];
					}
				}
				bondersCoeff[v]=dim - 1;

				bendersConstr=x[0] - x[0];
				for(int j=0; j<dim; j++){
					int vIndexComp=compVertices[nComp][j];
					IloInt compCoeff=static_cast<IloInt>(IloMin(bondersCoeff[vIndexComp], dim - cap));
					//IloInt compCoeff=bondersCoeff[vIndexComp];
					bendersConstr+=compCoeff*x[vIndexComp];
					//cout << compCoeff << "x" << vIndexComp << "+ ";
				}
				//cout << endl;
				add(bendersConstr>=dim - cap);
				compCuts++;
				bendersConstr.clear();
			}

			nComp++;
			dim=0;
		}
	}
	bendersConstr.end();
	////////////////////////////////////////////////////////////////////////////////////////////////

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

void BendersCutModel()
{

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	statusDFS = new int[vertex_number];
	vindex = new int[vertex_number];
	vweight = new double[vertex_number];
	bondersCoeff = new int[vertex_number];
	parent = new int[vertex_number];
	visited = new bool[vertex_number];

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
		lazySolx = IloNumArray(env, vertex_number);

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

		/*for(int i=0; i<vertex_number; i++){
			for(int j=0; j<vertex_number; j++){
				cout << symmAdjMatrix[i][j] << ", ";
			}
			cout << endl;
		}*/

		IloCplex::Callback sec;
		sec = cplex.use(LazyBendersCallback(env, x));

		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, treeLimit);
		//cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);

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
		obj.end();
		x.end();
		lazySolx.end();
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

	delete [] statusDFS;
	delete [] vindex;
	delete [] vweight;
	delete [] bondersCoeff;
	delete [] parent;
	delete [] visited;

	compVertices.end();
	weights.end();
}