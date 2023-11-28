//
// Created by paolo on 19/12/18.
//

#include "Heuristic.h"
#include "global_variables.h"
#include "global_functions.h"
#include "Graph_v4.h"

double* score;

void readScore(char * inFile){
	char * pch;
	pch = strtok (inFile,"/");
	for(int i=0; i<2; i++)
		pch = strtok (NULL, "/");
	pch = strtok (pch, ".");

	string scoreString;
	switch(option){
		case 1 :
			scoreString = "R/"+string(pch)+".betweenness.txt";
			break;
		case 2 :
			scoreString = "R/"+string(pch)+".closeness.txt";
			break;
		case 3 :
			scoreString = "R/"+string(pch)+".degree.txt";
			break;
		case 4 :
			scoreString = "R/"+string(pch)+".eigen.txt";
			break;
		default:
			scoreString = "R/"+string(pch)+".betweenness.txt";
			break;
	}

	ifstream in(scoreString);
	if(!in)
	{
		cout << "File of Scores could not be opened. I use Degree. " << endl;
		//exit(1);
		for(int i=0; i< vertex_number; i++)
			score[i] = G->DT[i];
		option = 3;
		return;
	}
	string line;
	getline(in, line);

	for(int i=0; i<vertex_number; i++){
		getline(in, line, ';');
		in >> score[i];
	}
	in.close();
}


int heuristic(char * inFile){

	score = new double[vertex_number];
	readScore(inFile);

	int totRemoved=0;
	bool feasible=false;
	statusDFS = new int[vertex_number];
	removed = new bool[vertex_number];
	for(int i=0; i<vertex_number; i++)
		removed[i]=false;

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
						if(score[check] > score[best])
							best=check;
					}
					if(option!=1 && option!=2 && option!=4){
						for(int n=0; n<G->DT[best]; n++)
							score[G->neighb[best][n]]--;
					}
					removed[best]=true;
					totRemoved++;
				}
				nComp++;
				dim=0;
			}
		}
		compVertices.clear();
	}

	delete [] score;
	compVertices.clear();

	output << totRemoved  << "\t";

	return totRemoved;
}