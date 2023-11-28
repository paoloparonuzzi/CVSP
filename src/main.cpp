#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <iomanip> // for precision
#include <map>
#include "ilcplex/ilocplex.h"

using namespace std;

#include "Compact.h"
#include "DegreeCut.h"
#include "CompactModel.h"
#include "Graph_v4.h"
#include "global_functions.h"
#include "global_variables.h"
#include "CLIQUE.h"
#include "ComponentCut.h"
#include "ConnectivityCut.h"
#include "MixedModel.h"
#include "Heuristic.h"
#include "BendersCut.h"
#include "MixedModelBFS.h"
#include "newProblemCompCut.h"

int main(int argc, char** argv) {

	double density;
	const char* fileName;
	istname=new char[2000];
	if (argc == 12) {
		strcpy(istname, argv[1]);
		vertex_number=atoi(argv[2]);
		density=atof(argv[3]);
		algorithm=atoi(argv[4]);
		option=atoi(argv[5]);
		k=atoi(argv[6]);
		//q=atoi(argv[7]);
		//b=atoi(argv[8]);
		tolerance=atof(argv[7]);
		treeLimit=atoi(argv[8]);
		timeLimit=atof(argv[9]);
		//seed=atoi(argv[10]);
		userFreq=atoi(argv[10]);
		fileName = argv[11]; 
		//srand(seed);
	}
	else {cout << "ERROR NUMBER STANDARD PARAMETERS" << endl; exit(2);}

	cout << "\n\nINSTANCE: ->\t" <<  istname << endl ;
	cout << "algorithm: ->\t" <<  algorithm << endl ;
	cout << "k: ->\t" <<  k << endl ;
	cout << "q: ->\t" <<  q << endl ;
	cout << "b: ->\t" <<  b << endl ;
	cout << "timeLimit: ->\t" <<  timeLimit << endl ;

	int *heads=NULL;
	int *tails=NULL;
	double *weights_arcs=NULL;
	double *weights_nodes=NULL;
	string nameFile(istname);
	////////////////////////////////RANDOM////////////////////////////////

	if(!strcmp(argv[1],"NULL")){
		cout << "\n\nRANDOM INSTANCE\n\n";

		cout << "\n\nRANDOM INSTANCE\n\n";
		heads=new int[10000000];
		tails=new int[10000000];
		weights_arcs=new double[10000000];
		weights_nodes=new double[10000000];

		edge_number=random_undirected_graphFF_data_load_density(vertex_number,density,tails,heads,weights_nodes,weights_arcs,1,1);
	}

	////////////////////////////////DIMACS////////////////////////////////
	//the heads and tails cannot be dimensioned inside the functions and in this
	//point we still don't know the dimensions

	heads=new int[10000000];
	tails=new int[10000000];
	weights_arcs=new double[10000000];
	weights_nodes=new double[10000000];

	ReadDIMACSFile(istname,&vertex_number,&edge_number,tails,heads,false);

	for(int i=0;i<edge_number;i++){weights_arcs[i]=1.0;}
	for(int i=0;i<vertex_number;i++){weights_nodes[i]=1.0;}

	cout << "GRAPH BUILDING\n";
	G = buildGraphFF(vertex_number,edge_number,heads,tails,weights_nodes,weights_arcs,1);
	cout << "DONE\n";
	//printAM(G);
	/////////////////////////////////////////////////
	cout << "COMPLEMENTARY GRAPH BUILDING\n";
	graphFF G_bar = buildComplementaryGraphFF_undirected(G,1);
	cout << "DONE\n";
	/////////////////////////////////////////////////

	cout << "NODES\t" << G->n << endl;
	cout << "ARCS\t" << G->m << endl;

	int FSerr = CheckFS(G);
	printf("Errors FS: %1d\n", FSerr);

	int BSerr = CheckBS(G);
	printf("Errors BS: %1d\n", BSerr);

	delete[] heads;
	delete[] tails;
	delete[] weights_nodes;
	delete[] weights_arcs;

	//printGraph(G);

	//fix the capacity value
	cap=ceil( (double) vertex_number/k);
	if(cap<=1)
	    return -1;

	/*for(int i=0; i<vertex_number; i++){
		for(int j=0; j<vertex_number; j++)
			cout << G->AMatrix[j][i]+G->AMatrix[i][j] << ",";
		cout << endl;
	}
	cout << endl;*/

	cout << "capacity: " << cap << endl;
	output.open(fileName, std::ios_base::app);

	if(algorithm==9 || algorithm==10 || algorithm==11 || algorithm==12)
		k =static_cast<int>(ceil(vertex_number*((double)k/100.0)));

	cout << "k (or budget): " << k << endl;

	if(output.is_open()){
		output << istname << "\t" << G->n << "\t" 
			<< G->m << "\t" << algorithm << "\t" << option << "\t" << k << "\t" << cap << "\t" << tolerance << "\t" << userFreq << "\t";
	}

	////////////////////////////////////////////////////////////////////////
	if(algorithm==1)
	{
		cout << "\n---->Compact Formulation\n\n";
		clock_t time_start_CPLEX=clock();
		if(option == 0)
            CompactModelToCompare();
		if(option == 1)		
			LP_CompactModel();
		clock_t time_end_CPLEX=clock();
		double computation_time_CPLEX=(double)(time_end_CPLEX-time_start_CPLEX)/(double)CLOCKS_PER_SEC;
		cout << "\nTOT TIME:\t" << computation_time_CPLEX << endl << endl;
	}
	if(algorithm==2)
	{
		cout << "\n---->Capacity - Degree Formulation\n\n";
		clock_t time_start_CPLEX=clock();
		DegreeCutModel();
		clock_t time_end_CPLEX=clock();
		double computation_time_CPLEX=(double)(time_end_CPLEX-time_start_CPLEX)/(double)CLOCKS_PER_SEC;
		cout << "\nTOT TIME:\t" << computation_time_CPLEX << endl << endl;
	}
	if(algorithm==3)
	{
		cout << "\n---->Capacity - Component Formulation\n\n";
		clock_t time_start_CPLEX=clock();
		ComponentCutModel();
		clock_t time_end_CPLEX=clock();
		double computation_time_CPLEX=(double)(time_end_CPLEX-time_start_CPLEX)/(double)CLOCKS_PER_SEC;
		cout << "\nTOT TIME:\t" << computation_time_CPLEX << endl << endl;
	}
	if(algorithm==4)
	{
		cout << "\n---->Capacity - Connectivity Formulation\n\n";
		clock_t time_start_CPLEX=clock();
		ConnectivityCutModel();
		clock_t time_end_CPLEX=clock();
		double computation_time_CPLEX=(double)(time_end_CPLEX-time_start_CPLEX)/(double)CLOCKS_PER_SEC;
		cout << "\nTOT TIME:\t" << computation_time_CPLEX << endl << endl;
	}
	if(algorithm==5)
	{
		if(option>0)	heuristic(istname);
		cout << "\n---->Capacity - Mixed Formulation\n\n";
		clock_t time_start_CPLEX=clock();
		MixedModelBFS();
		clock_t time_end_CPLEX=clock();
		double computation_time_CPLEX=(double)(time_end_CPLEX-time_start_CPLEX)/(double)CLOCKS_PER_SEC;
		cout << "\nTOT TIME:\t" << computation_time_CPLEX << endl << endl;
	}
	int preProces=0;
	if(algorithm==6)
	{
		cout << "\n---->Capacity - Heuristic Algorithm\n\n";
		clock_t time_heur_start=clock();
		int bound = heuristic(istname);
		int preProces=0;
		for(int i=0; i<vertex_number; i++){
			if(G->DT[i] - bound > cap)
				preProces++;
		}

		clock_t time_heur_end=clock();
		double time_heur=(double)(time_heur_end-time_heur_start)/(double)CLOCKS_PER_SEC;
		cout << "\nTOT TIME:\t" << time_heur << endl << endl;
		output << time_heur;
	}
	if(algorithm==7)
	{
		cout << "\n---->Capacity - Benders Formulation\n\n";
		clock_t time_start_CPLEX=clock();
		BendersCutModel();
		clock_t time_end_CPLEX=clock();
		double computation_time_CPLEX=(double)(time_end_CPLEX-time_start_CPLEX)/(double)CLOCKS_PER_SEC;
		cout << "\nTOT TIME:\t" << computation_time_CPLEX << endl << endl;
	}
	if(algorithm==8)
	{
		cout << "\n---->Compare Models\n\n";
		clock_t time_start_CPLEX=clock();
        CompactModelWithClique();
		clock_t time_end_CPLEX=clock();
		double computation_time_CPLEX=(double)(time_end_CPLEX-time_start_CPLEX)/(double)CLOCKS_PER_SEC;
		cout << "\nTOT TIME:\t" << computation_time_CPLEX << endl << endl;
	}
	if(algorithm==9)
	{

		cout << "\n---->Min Max Size - Component Cut Formulation\n\n";
		clock_t time_start_CPLEX=clock();
		minMaxC_CompCutModel();
		clock_t time_end_CPLEX=clock();
		double computation_time_CPLEX=(double)(time_end_CPLEX-time_start_CPLEX)/(double)CLOCKS_PER_SEC;
		cout << "\nTOT TIME:\t" << computation_time_CPLEX << endl << endl;
	}
	if(algorithm==10)
	{
		cout << "\n---->Min Max Size - Compact Clique Formulation\n\n";
		clock_t time_start_CPLEX=clock();
		minMaxC_CompactClique();
		clock_t time_end_CPLEX=clock();
		double computation_time_CPLEX=(double)(time_end_CPLEX-time_start_CPLEX)/(double)CLOCKS_PER_SEC;
		cout << "\nTOT TIME:\t" << computation_time_CPLEX << endl << endl;
	}
	if(algorithm==11)
	{
		cout << "\n---->Min Max Size - Compact Edge Formulation\n\n";
		clock_t time_start_CPLEX=clock();
		minMaxC_CompactEdge();
		clock_t time_end_CPLEX=clock();
		double computation_time_CPLEX=(double)(time_end_CPLEX-time_start_CPLEX)/(double)CLOCKS_PER_SEC;
		cout << "\nTOT TIME:\t" << computation_time_CPLEX << endl << endl;
	}
	if(algorithm==12)
	{
		cout << "\n---->Min Max Size - Compact Follower Formulation\n\n";
		clock_t time_start_CPLEX=clock();
		minMaxC_CompactFollower();
		clock_t time_end_CPLEX=clock();
		double computation_time_CPLEX=(double)(time_end_CPLEX-time_start_CPLEX)/(double)CLOCKS_PER_SEC;
		cout << "\nTOT TIME:\t" << computation_time_CPLEX << endl << endl;
	}
	/////////////////////////////////////////////////////////////////////////



	deleteGraphFF(G);
	deleteGraphFF(G_bar);

	delete [] istname;
	printf("\nDONE!");

	if(output.is_open()){
		output << endl;
		output.close();
	}
	////////////////////////////////////////////////////////////////////////
	return 1;
}
