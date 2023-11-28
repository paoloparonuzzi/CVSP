#include "CompactModel.h"
#include "Graph_v4.h"

void LP_CompactModel()
{
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try{
		int i,j,e;
		int nVar = vertex_number*k;
		IloNumVarArray x(env, nVar, 0.0, 1.0);

		//Objective Function: maximize the number of vertices assigned to some component
		IloObjective obj(env, IloSum(x), IloObjective::Maximize);
		model.add(obj);

		//Contstraint: each vertex can be assigned at most one component
		IloExpr vertexConstr(env);
		for(i=0; i<vertex_number; i++){		
			for(j=i; j<nVar; j=j+vertex_number)
				vertexConstr = vertexConstr + x[j];
			model.add(vertexConstr <= 1);
			vertexConstr.clear();
		}
		vertexConstr.end();

		//Constraint: for each component, the capacity constraint must be respected		
		IloExpr capConstr(env);
		for(i=0; i<nVar; i=i+vertex_number){		
			for(j=i; j<i+vertex_number; j++)
				capConstr = capConstr + x[j];
			model.add(capConstr <= cap);
			capConstr.clear();
		}
		capConstr.end();


		//Constraint: each pair of vertices defining an edge can be assigned to at most one component
		IloExpr cliqueConstr(env);
		for(e=0; e<edge_number; e++){
			for(i=0; i<k; i++){
				cliqueConstr = x[G->T[e]+i*vertex_number];
				for(j=0; j<k; j++){
					if(j != i)
						cliqueConstr += x[G->H[e]+j*vertex_number];
				}
				model.add(cliqueConstr <= 1);
				cliqueConstr.clear();
			}
		}
		cliqueConstr.end();		

		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, treeLimit);
		cplex.setParam(IloCplex::NodeLim, 1);
		//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
		//cplex.setParam(IloCplex::PreInd, IloFalse);
		
		clock_t time1 = clock(); //cplex.getTime();
		cplex.solve();
		clock_t time2 = clock(); //cplex.getTime();
		double time = (double) (time2 - time1)/CLOCKS_PER_SEC;

		/*for(int i = 0; i<vertex_number; i++)
			cout << "x" << i << ": " << cplex.getValue(x[i]) << "\t";
		cout << endl;
		cout << "obj: " << cplex.getObjValue() << endl;
		cout << "user cuts: " << userCuts << endl;
		cout << "lazy cuts: " << lazyCuts << endl;
		cout << "dual obj: " << vertex_number - cplex.getObjValue() << endl;*/
		
		if(cplex.getStatus() == IloAlgorithm::Unknown || cplex.getStatus() == IloAlgorithm::Infeasible){		
			IloNum bound = vertex_number - cplex.getBestObjValue();			
			IloInt nNode = cplex.getNnodes();		
			IloInt status = cplex.getStatus();
			IloInt nVar = cplex.getNcols();
			IloInt nConst = cplex.getNrows();		
			if(output.is_open())
			{
				output.precision(10);
				output << "null" << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" 
				<< nConst << "\t" << nNode;
			}
		}			
		else{
			IloNum objective = vertex_number - cplex.getObjValue();
			IloNum bound = vertex_number - cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();		
			IloInt status = cplex.getStatus();
			IloInt nVar = cplex.getNcols();
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << objective << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" 
				<< nConst << "\t" << nNode;
			}		
		}

		obj.end();
		x.end();	
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
			output << "null" << "\t" << bound << "\t"  // << rootObjective << "\t" 
				<< time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" << nConst << "\t" << nNode;
		}
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	cplex.end();
	model.end();
	env.end();
}



int index(int i, int j)
    {return j + k * i;}

int newIndex(int i, int j)
	{return j + (vertex_number-k) * i;}

void compactEdgeModel()
{
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try{
		int nX = vertex_number*k;
		int nY = edge_number*k;
		IloIntVarArray x(env, nX, 0, 1);
		IloIntVarArray y(env, nY, 0, 1);

		//Objective Function: maximize the number of vertices out of the separator
		IloObjective obj(env, IloSum(x), IloObjective::Maximize);
		model.add(obj);

		//Contstraint: each vertex can be assigned at most one component
		IloExpr vertexConstr(env);
		for(int i=0; i<vertex_number; i++){
			for(int j=0; j<k; j++)
				vertexConstr = vertexConstr + x[index(i,j)];
			model.add(vertexConstr <= 1);
			vertexConstr.clear();
		}
		vertexConstr.end();

		//Contstraint: each edge can be assigned at most one component
		IloExpr edgeConstr(env);
		for(int i=0; i<edge_number; i++){
			for(int j=0; j<k; j++)
				edgeConstr += y[index(i,j)];
			model.add(edgeConstr <= 1);
			edgeConstr.clear();
		}
		edgeConstr.end();
		
	
		//Constraint: for each component, the capacity constraint must be respected		
		IloExpr capConstr(env);
		for(int j=0; j<k; j++){
			for(int i=0; i<vertex_number; i++)
				capConstr += x[index(i,j)];
			model.add(capConstr <= cap);
			capConstr.clear();
		}
		capConstr.end();

		//Constraint: each edge can be assigned to at most one component
		for(int i=0; i<edge_number; i++){
			for(int j=0; j<k; j++){
				model.add(x[index(G->T[i],j)] - y[index(i,j)] <= 0);
				model.add(x[index(G->H[i],j)] - y[index(i,j)] <= 0);
			}
		}
		
		
		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, treeLimit);
		//cplex.setParam(IloCplex::NodeLim, 1);
		//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
		//cplex.setParam(IloCplex::PreInd, IloFalse);
		
		clock_t time1 = clock(); //cplex.getTime();
		cplex.solve();
		clock_t time2 = clock(); //cplex.getTime();
		double time = (double) (time2 - time1)/CLOCKS_PER_SEC;
		
		if(cplex.getStatus() == IloAlgorithm::Unknown || cplex.getStatus() == IloAlgorithm::Infeasible){		
			IloNum bound = vertex_number - cplex.getBestObjValue();			
			IloInt nNode = cplex.getNnodes();		
			IloInt status = cplex.getStatus();
			IloInt nVar = cplex.getNcols();
			IloInt nConst = cplex.getNrows();		
			if(output.is_open())
			{
				output.precision(10);
				output << "null" << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" 
				<< nConst << "\t" << nNode;
			}
		}			
		else{
			IloNum objective = vertex_number - cplex.getObjValue();
			IloNum bound = vertex_number - cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();		
			IloInt status = cplex.getStatus();
			IloInt nVar = cplex.getNcols();
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << objective << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" 
				<< nConst << "\t" << nNode;
			}		
		}
		
		//cout << cplex.getStatus() << endl;
		//cplex.exportModel("model.lp");

		obj.end();
		x.end();	
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
			output << "null" << "\t" << bound << "\t"  // << rootObjective << "\t" 
				<< time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" << nConst << "\t" << nNode;
		}
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	cplex.end();
	model.end();
	env.end();
}

void CompactModelToCompare()
{
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try{
		int nVar = vertex_number*k;
		IloIntVarArray x(env, nVar, 0, 1);

		//Objective Function: maximize the number of vertices assigned to some component
		IloObjective obj(env, vertex_number - IloSum(x), IloObjective::Minimize);
		model.add(obj);

		//Contstraint: each vertex can be assigned to at most one component
		IloExpr vertexConstr(env);
		for(int i=0; i<vertex_number; i++){
			for(int j=0; j<k; j++)
				vertexConstr = vertexConstr + x[index(i,j)];
			model.add(vertexConstr <= 1);
			vertexConstr.clear();
		}
		vertexConstr.end();

		//Constraint: for each component, the capacity constraint must be respected
		IloExpr capConstr(env);
		for(int j=0; j<k; j++){
			for(int i=0; i<vertex_number; i++)
				capConstr = capConstr + x[index(i,j)];
			model.add(capConstr <= cap);
			//model.add(capConstr >= 1);
			capConstr.clear();
		}
		capConstr.end();


		//Constraint: each pair of vertices defining an edge must be assigned to at most one component
		IloExpr shoreConstr(env);
		for(int e=0; e<edge_number; e++){
			for(int i=0; i<k; i++){
				shoreConstr = x[index(G->T[e],i)];
				for(int j=0; j<k; j++){
					if(j != i)
						shoreConstr += x[index(G->H[e],j)];
				}
				model.add(shoreConstr <= 1);
				shoreConstr.clear();
			}
		}
		shoreConstr.end();

		string modelName = istname+string("_E_k")+to_string(k)+string(".mps");
		cplex.exportModel(modelName.c_str());
		return;

		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, treeLimit);
		//cplex.setParam(IloCplex::NodeLim, 1);
		//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
		//cplex.setParam(IloCplex::PreInd, IloFalse);

		clock_t time1 = clock(); //cplex.getTime();
		cplex.solve();
		clock_t time2 = clock(); //cplex.getTime();
		double time = (double) (time2 - time1)/CLOCKS_PER_SEC;

		cout << "OBJ=" << cplex.getObjValue() << ";" << endl;
		cout << "Separator: ";
		for(int i=0; i<vertex_number; i++){
			int count=0;
			for(int j=0; j<k; j++){
				if(cplex.getValue(x[index(i,j)])>=0.8)
					count++;
			}
			if(count<=0.8)
				cout << i << " ";
		}
		cout << endl;
		for(int j=0; j<k; j++){
			cout << "c" << j+1 << ": ";
			for(int i=0; i<vertex_number; i++){
				if(cplex.getValue(x[index(i,j)])>=0.8)
					cout << i << " ";
			}
			cout << endl;
		}

		vector<int> xSol;
		for(int i=0; i<vertex_number; i++){
			for(int j=0; j<k; j++){
				xSol.push_back(static_cast<int &&>(cplex.getValue(x[index(i, j)])+0.5));
			}
		}
		//printSolution(G, xSol);

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
					   << nConst << "\t" << nNode;
			}
		}
		else{
			IloNum objective = cplex.getObjValue();
			IloNum bound = cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();
			IloInt status = cplex.getStatus();
			IloInt nVar = cplex.getNcols();
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << objective << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t"
					   << nConst << "\t" << nNode;
			}
		}

		obj.end();
		x.end();
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
			output << "null" << "\t" << bound << "\t"  // << rootObjective << "\t"
				   << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" << nConst << "\t" << nNode;
		}
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	cplex.end();
	model.end();
	env.end();
}

void CompactModelWithClique()
{

	clock_t time_start_clique=clock();
	int nClique = buildCliqueEdgeCover(G);
	clock_t time_end_clique=clock();
	double time_clique = (double) (time_end_clique - time_start_clique)/CLOCKS_PER_SEC;
	//output << time_clique << "\t";

	cout << "number of clique: " << nClique << endl;
	for(int c=0; c<nClique; c++){
		cout << c << ": ";
		for(int s : G->edgeClique[c])
			cout << s << " ";
		cout << endl;
	}
	//output << nClique;

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try{
		int nVar = vertex_number*k;
		int nNewVar = nClique*k;
		IloIntVarArray x(env, nVar, 0, 1);
		IloIntVarArray y(env, nNewVar, 0, 1);

		//Objective Function: maximize the number of vertices assigned to some component
		IloObjective obj(env, vertex_number - IloSum(x), IloObjective::Minimize);
		model.add(obj);

		//Contstraint: each clique can be assigned to at most one component
		IloExpr cliqueConstr(env);
		for(int i=0; i<nClique; i++){
			for(int j=0; j<k; j++)
				cliqueConstr = cliqueConstr + y[index(i,j)];
			model.add(cliqueConstr <= 1);
			cliqueConstr.clear();
		}
		cliqueConstr.end();

		//Isolated vertices can be assigned to at most one component
		IloExpr isolatedConstr(env);
		for(int i=0; i<vertex_number; i++){
			if(G->DT[i]==0){
				for(int j=0; j<k; j++)
					isolatedConstr = isolatedConstr + x[index(i,j)];
				model.add(isolatedConstr <= 1);
				isolatedConstr.clear();
			}
		}
		isolatedConstr.end();

		//Constraint: for each component, the capacity constraint must be respected
		IloExpr capConstr(env);
		for(int j=0; j<k; j++){
			for(int i=0; i<vertex_number; i++)
				capConstr = capConstr + x[index(i,j)];
			model.add(capConstr <= cap);
			capConstr.clear();
		}
		capConstr.end();


		//Constraint: each vertex can be assignaed to some component only if its clique is assigned to the same component
		for(int j=0; j<k; j++){
			for(int q=0; q<nClique; q++){
				for(int i : G->edgeClique[q]){
					model.add(x[index(i,j)] <= y[index(q, j)]);
				}
			}
		}

        string modelName = istname+string("_C_k")+to_string(k)+string(".mps");
        cplex.exportModel(modelName.c_str());
        return;

		//cplex.setOut(env.getNullStream());
		//cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, treeLimit);
		//cplex.setParam(IloCplex::NodeLim, 1);
		//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
		//cplex.setParam(IloCplex::PreInd, IloFalse);

		clock_t time1 = clock(); //cplex.getTime();
		cplex.solve();
		clock_t time2 = clock(); //cplex.getTime();
		double time = (double) (time2 - time1)/CLOCKS_PER_SEC;

		cout << "OBJ=" << cplex.getObjValue() << ";" << endl;
		cout << "Separator: ";
		for(int i=0; i<vertex_number; i++){
			int count=0;
			for(int j=0; j<k; j++){
				if(cplex.getValue(x[index(i,j)])>=0.8)
					count++;
			}
			if(count<=0.8)
				cout << i << " ";
		}
		cout << endl;
		for(int j=0; j<k; j++){
			cout << "c" << j+1 << ": ";
			for(int i=0; i<vertex_number; i++){
				if(cplex.getValue(x[index(i,j)])>=0.8)
					cout << i << " ";
			}
			cout << endl;
		}

		if(cplex.getStatus() == IloAlgorithm::Unknown || cplex.getStatus() == IloAlgorithm::Infeasible){
			IloNum bound = vertex_number - cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();
			IloInt status = cplex.getStatus();
			IloInt nVars = cplex.getNcols();
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << "null" << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVars << "\t"
					   << nConst << "\t" << nNode;
			}
		}
		else{
			IloNum objective = cplex.getObjValue();
			IloNum bound = cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();
			IloInt status = cplex.getStatus();
			IloInt nVars = cplex.getNcols();
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << objective << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVars << "\t"
					   << nConst << "\t" << nNode;
			}
		}

		obj.end();
		x.end();
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
			output << "null" << "\t" << bound << "\t"  // << rootObjective << "\t"
				   << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" << nConst << "\t" << nNode;
		}
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	cplex.end();
	model.end();
	env.end();
}

void minMaxC_CompactClique()
{

	clock_t time_start_clique=clock();
	int nClique = buildCliqueEdgeCover(G);
	clock_t time_end_clique=clock();
	double time_clique = (double) (time_end_clique - time_start_clique)/CLOCKS_PER_SEC;
	//output << time_clique << "\t";

	/*cout << "number of clique: " << nClique << endl;
	for(int c=0; c<nClique; c++){
		cout << c << ": ";
		for(int s : G->edgeClique[c])
			cout << s << " ";
		cout << endl;
	}*/
	//output << nClique;

	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try{
		int nVar = vertex_number*(vertex_number-k);
		int nNewVar = nClique*(vertex_number-k);
		IloIntVarArray x(env, nVar, 0, 1);
		IloIntVarArray y(env, nNewVar, 0, 1);
		IloNumVar lambda(env, 0, vertex_number-k);

		//Objective Function: minimize the size of the largest component (lambda)
		IloObjective obj(env, lambda, IloObjective::Minimize);
		model.add(obj);

		//Contstraint: the budget must be respected
		model.add(IloSum(x) >= vertex_number - k);

		//Contstraint: each clique can be assigned to at most one component
		IloExpr cliqueConstr(env);
		for(int i=0; i<nClique; i++){
			for(int j=0; j<vertex_number-k; j++)
				cliqueConstr = cliqueConstr + y[newIndex(i,j)];
			model.add(cliqueConstr <= 1);
			cliqueConstr.clear();
		}
		cliqueConstr.end();

		//Isolated vertices can be assigned to at most one component
		IloExpr isolatedConstr(env);
		for(int i=0; i<vertex_number; i++){
			if(G->DT[i]==0){
				for(int j=0; j<vertex_number-k; j++)
					isolatedConstr = isolatedConstr + x[newIndex(i,j)];
				model.add(isolatedConstr <= 1);
				isolatedConstr.clear();
			}
		}
		isolatedConstr.end();

		//Constraint: for each component, lambda must be greater or equal to its size
		IloExpr capConstr(env);
		for(int j=0; j<vertex_number-k; j++){
			for(int i=0; i<vertex_number; i++)
				capConstr = capConstr + x[newIndex(i,j)];
			model.add(capConstr <= lambda);
			capConstr.clear();
		}
		capConstr.end();


		//Constraint: each vertex can be assignaed to some component only if its clique is assigned to the same component
		for(int j=0; j<vertex_number-k; j++){
			for(int q=0; q<nClique; q++){
				for(int i : G->edgeClique[q]){
					model.add(x[newIndex(i,j)] <= y[newIndex(q, j)]);
				}
			}
		}

		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, treeLimit);
		//cplex.setParam(IloCplex::NodeLim, 1);
		//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
		//cplex.setParam(IloCplex::PreInd, IloFalse);
        //////////////////////////////////////////////////////////////
        cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 1);

		clock_t time1 = clock(); //cplex.getTime();
		cplex.solve();
		clock_t time2 = clock(); //cplex.getTime();
		double time = (double) (time2 - time1)/CLOCKS_PER_SEC;

		cout << "OBJ=" << cplex.getObjValue() << ";" << endl;
		cout << "Separator: ";
		for(int i=0; i<vertex_number; i++){
			int count=0;
			for(int j=0; j<vertex_number-k; j++){
				if(cplex.getValue(x[newIndex(i,j)])>=0.8)
					count++;
			}
			if(count<=0.8)
				cout << i << " ";
		}
		cout << endl;
		for(int j=0; j<vertex_number-k; j++){
			cout << "c" << j+1 << ": ";
			for(int i=0; i<vertex_number; i++){
				if(cplex.getValue(x[newIndex(i,j)])>=0.8)
					cout << i << " ";
			}
			cout << endl;
		}

		if(cplex.getStatus() == IloAlgorithm::Unknown || cplex.getStatus() == IloAlgorithm::Infeasible){
			IloNum bound = cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();
			IloInt status = cplex.getStatus();
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << "null" << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar+nNewVar << "\t"
					   << nConst << "\t" << nNode;
			}
		}
		else{
			IloNum objective = cplex.getObjValue();
			IloNum bound = cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();
			IloInt status = cplex.getStatus();
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << objective << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar+nNewVar << "\t"
					   << nConst << "\t" << nNode;
			}
		}

		obj.end();
		x.end();
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
			output << "null" << "\t" << bound << "\t"  // << rootObjective << "\t"
				   << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" << nConst << "\t" << nNode;
		}
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	cplex.end();
	model.end();
	env.end();
}

void minMaxC_CompactEdge()
{
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try{
		int nVar = vertex_number*(vertex_number-k);
		IloIntVarArray x(env, nVar, 0, 1);
		IloIntVar lambda(env, 0, vertex_number-k);

		//Objective Function: minimize the size of the largest component (lambda)
		IloObjective obj(env, lambda, IloObjective::Minimize);
		model.add(obj);

		//Contstraint: the budget must be respected
		model.add(IloSum(x) >= vertex_number - k);

		//Contstraint: each vertex can be assigned to at most one component
		IloExpr vertexConstr(env);
		for(int i=0; i<vertex_number; i++){
			for(int j=0; j<(vertex_number-k); j++)
				vertexConstr = vertexConstr + x[newIndex(i,j)];
			model.add(vertexConstr <= 1);
			vertexConstr.clear();
		}
		vertexConstr.end();

		//Constraint: each component must have cardinality lower than lambda
		IloExpr capConstr(env);
		for(int j=0; j<(vertex_number-k); j++){
			for(int i=0; i<vertex_number; i++)
				capConstr = capConstr + x[newIndex(i,j)];
			model.add(capConstr <= lambda);
			capConstr.clear();
		}
		capConstr.end();


		//Constraint: each pair of vertices defining an edge must be assigned to at most one component
		IloExpr shoreConstr(env);
		for(int e=0; e<edge_number; e++){
			for(int i=0; i<(vertex_number-k); i++){
				shoreConstr = x[newIndex(G->T[e],i)];
				for(int j=0; j<(vertex_number-k); j++){
					if(j != i)
						shoreConstr += x[newIndex(G->H[e],j)];
				}
				model.add(shoreConstr <= 1);
				shoreConstr.clear();
			}
		}
		shoreConstr.end();

		cplex.setOut(env.getNullStream());
		cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, treeLimit);
		//cplex.setParam(IloCplex::NodeLim, 1);
		//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
		//cplex.setParam(IloCplex::PreInd, IloFalse);

		clock_t time1 = clock(); //cplex.getTime();
		cplex.solve();
		clock_t time2 = clock(); //cplex.getTime();
		double time = (double) (time2 - time1)/CLOCKS_PER_SEC;

		cout << "OBJ=" << cplex.getObjValue() << ";" << endl;
		cout << "Separator: ";
		for(int i=0; i<vertex_number; i++){
			int count=0;
			for(int j=0; j<(vertex_number-k); j++){
				if(cplex.getValue(x[newIndex(i,j)])>=0.8)
					count++;
			}
			if(count<=0.8)
				cout << i << " ";
		}
		cout << endl;
		for(int j=0; j<(vertex_number-k); j++){
			cout << "c" << j+1 << ": ";
			for(int i=0; i<vertex_number; i++){
				if(cplex.getValue(x[newIndex(i,j)])>=0.8)
					cout << i << " ";
			}
			cout << endl;
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
					   << nConst << "\t" << nNode;
			}
		}
		else{
			IloNum objective = cplex.getObjValue();
			IloNum bound = cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();
			IloInt status = cplex.getStatus();
			IloInt nVar = cplex.getNcols();
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << objective << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t"
					   << nConst << "\t" << nNode;
			}
		}

		obj.end();
		x.end();
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
			output << "null" << "\t" << bound << "\t"  // << rootObjective << "\t"
				   << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" << nConst << "\t" << nNode;
		}
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	cplex.end();
	model.end();
	env.end();
}

int sigmaIndex(int i, int j)
{return j + vertex_number * i;}

void minMaxC_CompactFollower()
{
	IloEnv env;
	IloModel model(env);
	IloCplex cplex(model);

	try{
		IloIntVarArray x(env, vertex_number, 0, 1);
		IloNumVarArray sigma(env, vertex_number*vertex_number, 0, 1);
		IloNumVar lambda(env, 0, vertex_number-k);

		//Objective Function: minimize the size of the largest component (lambda)
		IloObjective obj(env, lambda, IloObjective::Minimize);
		model.add(obj);

		//Contstraint: the budget must be respected
		model.add(IloSum(x) <= k);


		//Constraint: each component must have cardinality lower than lambda
		IloExpr capConstr(env);
		for(int l=0; l<vertex_number; l++){
			for(int v=0; v<vertex_number; v++)
				capConstr = capConstr + sigma[sigmaIndex(v,l)];
			model.add(lambda >= capConstr);
			capConstr.clear();
		}
		capConstr.end();

		//Constraint: if two neighboring vertices v and w are not interdicted,
		//and v is in the same component as l , then w must be in the same component as well
		for(int v=0; v<vertex_number; v++){
			for(int w=v+1; w<vertex_number; w++){
				if(G->AMatrix[v][w] || G->AMatrix[w][v]){
					for(int l=0; l<vertex_number; l++){
						model.add(x[w] + x[v]>=sigma[sigmaIndex(v, l)] - sigma[sigmaIndex(w, l)]);
						model.add(x[w] + x[v]>=sigma[sigmaIndex(w, l)] - sigma[sigmaIndex(v, l)]);
					}
				}
			}
		}

		for(int l=0; l<vertex_number; l++)
			model.add(sigma[sigmaIndex(l,l)]==1);


		//cplex.setOut(env.getNullStream());
		//cplex.setWarning(env.getNullStream());
		cplex.setParam(IloCplex::TiLim, timeLimit);
		cplex.setParam(IloCplex::Threads, 1);
		cplex.setParam(IloCplex::Param::MIP::Limits::TreeMemory, treeLimit);
		//cplex.setParam(IloCplex::NodeLim, 1);
		//cplex.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
		//cplex.setParam(IloCplex::PreInd, IloFalse);
        //////////////////////////////////////////////////////////////
        cplex.setParam(IloCplex::Param::MIP::Tolerances::AbsMIPGap, 1);

		clock_t time1 = clock(); //cplex.getTime();
		cplex.solve();
		clock_t time2 = clock(); //cplex.getTime();
		double time = (double) (time2 - time1)/CLOCKS_PER_SEC;

		cout << "OBJ=" << cplex.getObjValue() << ";" << endl;
		cout << "Separator: ";
		for(int i=0; i<vertex_number; i++){
			if(cplex.getValue(x[i])>=0.8)
				cout << i << " ";
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
					   << nConst << "\t" << nNode;
			}
		}
		else{
			IloNum objective = cplex.getObjValue();
			IloNum bound = cplex.getBestObjValue();
			IloInt nNode = cplex.getNnodes();
			IloInt status = cplex.getStatus();
			IloInt nVar = cplex.getNcols();
			IloInt nConst = cplex.getNrows();
			if(output.is_open())
			{
				output.precision(10);
				output << objective << "\t" << bound << "\t" << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t"
					   << nConst << "\t" << nNode;
			}
		}

		obj.end();
		x.end();
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
			output << "null" << "\t" << bound << "\t"  // << rootObjective << "\t"
				   << time << "\t" << cplex.getStatus() << "\t" << nVar << "\t" << nConst << "\t" << nNode;
		}
	}
	catch (...) {
		cerr << "Error" << endl;
	}
	cplex.end();
	model.end();
	env.end();
}