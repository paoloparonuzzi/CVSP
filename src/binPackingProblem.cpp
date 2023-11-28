//
// Created by paolo on 14/09/18.
//

#include <ilconcert/ilomodel.h>
#include <ilcplex/ilocplexi.h>
#include "binPackingProblem.h"
#include "global_variables.h"
#include "global_functions.h"

bool lazyBinPackingModel(IloEnv env, int nComp){

	bool bbpSolved=false;

	IloModel BPP(env);
	IloCplex cplexBPP(BPP);

	//Each partition is a Bin that I can activate or not
	//The maximum number of Bin that I could need is the number of component
	IloIntVarArray y(env, nComp, 0, 1);

	//Objective Function: minimize the number of Bin used
	IloObjective obj(env, IloSum(y), IloObjective::Minimize);
	BPP.add(obj);

	//A variable to assign the different component to different bins
	IloIntVarArray z(env, nComp*nComp, 0, 1);
	IloIntArray iloWeights = IloIntArray(env, nComp);
	for(int i=0; i<nComp; i++)
		iloWeights[i] = weights[i];

	//Capacity constraint for each one of the nComp Bins (partitions)
	IloExpr binConstr(env);
	binConstr = y[0]-y[0];
	for(int i=0; i<nComp; i++){
		for(int j=0; j<nComp; j++)
			binConstr += iloWeights[j]*z[j+i*nComp];
		BPP.add(binConstr <= y[i]*cap);
		binConstr.clear();
	}
	binConstr.end();

	//Each component must be assigned to exactly one partition
	IloExpr assComp(env);
	assComp = z[0]-z[0];
	for(int i=0; i<nComp; i++){
		for(int j=0; j<nComp; j++)
			assComp += z[i+j*nComp];
		BPP.add(assComp == 1);
		assComp.clear();
	}
	assComp.end();

	cplexBPP.setOut(env.getNullStream());
	cplexBPP.setWarning(env.getNullStream());
	cplexBPP.setParam(IloCplex::TiLim, timeLimit);
	cplexBPP.setParam(IloCplex::Threads, 1);
	cplexBPP.setParam(IloCplex::Param::MIP::Limits::TreeMemory, treeLimit);

	//cout << "-----> Solving Bin Packing Problem" << endl;
	cplexBPP.solve();
	if(cplexBPP.getBestObjValue() <= k)	bbpSolved=true;
	//cout << "Used partitions: " << cplexBPP.getBestObjValue() << " k: " << k << endl;

	y.end();
	z.end();
	cplexBPP.end();
	BPP.end();

	return bbpSolved;
}