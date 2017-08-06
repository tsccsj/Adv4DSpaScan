/**
 * mc.c
 * Authors: Yizhao Gao <yizhaotsccsj@gmail.com>
 * Date: {08/06/2017}
 */

#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <omp.h>
#include "scan.h"

using namespace std;

void randomLabel(int * indAll, int casCount, int allCount) {
	static std::random_device rd;
	static std::mt19937 rng(rd());
	static std::uniform_int_distribution<int> uni(0, allCount - 1);

	int casID;
	for(int i = 0; i < allCount; i++)
		indAll[i] = 0;

	for(int i = 0; i < casCount; i++) {
		casID = uni(rng);
		while(indAll[casID] == 1)
			casID = uni(rng);
		indAll[casID] = 1;
	}

	return;
}

int * monteCarlo(double * x1, double * y1, double * x2, double * y2, int * locEnding, int locCount, int casCount, int allCount, int windowShape, double wSize, int wCount, int elimIntersectOD, int highLow, double * clusterLL, int nClusters, int nSim, int pAsCenter, double xMin, double yMin, double cellSize, int nRow, int nCol) {

	int * nExtreme;

	if(NULL == (nExtreme = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < nClusters; i++)
		nExtreme[i] = 0;

	int totalWindows;
	if(windowShape == 2) {
		totalWindows = wCount * wCount;
	}
	else {
		totalWindows = wCount;
	}
	if(pAsCenter != 2) {
		totalWindows = totalWindows * locCount;
	}
	else {
		totalWindows = totalWindows * nRow * nCol * nRow * nCol;
	}


	int * indAll;
	int * simCass;
	int * simCons;

	if(NULL == (indAll = (int *) malloc (allCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simCass = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simCons = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	int indID, simCas, simCon;

	int * simCasInW;
	int * simConInW;
	double * simll;
	
	if(NULL == (simCasInW = (int *) malloc (totalWindows * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simConInW = (int *) malloc (totalWindows * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simll = (double *) malloc (totalWindows * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}


	double simMaxLL;

	for(int i = 0; i < nSim; i++) {
		randomLabel(indAll, casCount, allCount);
	
		indID = 0;
		for(int j = 0; j < locCount; j++) {
			simCas = 0;
			simCon = 0;
			for(; indID < locEnding[j]; indID ++) {
				if(indAll[indID] == 1) {
					simCas ++;
				}
				else {
					simCon ++;
				}
			}
			simCass[j] = simCas;
			simCons[j] = simCon;
		}

		if(windowShape == 0) {
			if(pAsCenter == 2) {
				getCCCount4DSph(x1, y1, x2, y2, simCass, simCons, locCount, wSize, wCount, xMin, yMin, cellSize, nCol, nRow, simCasInW, simConInW, elimIntersectOD);
			}
			else {
				getCCCount4DSph(x1, y1, x2, y2, simCass, simCons, locCount, wSize, wCount, simCasInW, simConInW, elimIntersectOD);	
			}
		}
		else if(windowShape == 1) {
			if(pAsCenter == 2) {
				getCCCount2DSame(x1, y1, x2, y2, simCass, simCons, locCount, wSize, wCount, xMin, yMin, cellSize, nCol, nRow, simCasInW, simConInW, elimIntersectOD);
			}
			else {
				getCCCount2DSame(x1, y1, x2, y2, simCass, simCons, locCount, wSize, wCount, simCasInW, simConInW, elimIntersectOD);	
			}
		}
		else {
			if(pAsCenter == 2) {
				getCCCount2DDiff(x1, y1, x2, y2, simCass, simCons, locCount, wSize, wCount, xMin, yMin, cellSize, nCol, nRow, simCasInW, simConInW, elimIntersectOD);
			}
			else {
				getCCCount2DDiff(x1, y1, x2, y2, simCass, simCons, locCount, wSize, wCount, simCasInW, simConInW, elimIntersectOD);	
			}
		}

		loglikelihood(simll, simCasInW, simConInW, locCount * wCount, casCount, allCount - casCount, highLow);

		simMaxLL = 1;
		int k = 0;
		for(; k < totalWindows; k++) {
			if(simll[k] < 0) {
				simMaxLL = simll[k];
				k++;
				break;
			}
		}

		for(; k < totalWindows; k++) {
			if(simll[k] < 0 && simll[k] > simMaxLL) {
				simMaxLL = simll[k];
			}
		}

//		printf("%d:\t%lf\n", i, simMaxLL);

		if(simMaxLL < 0) {
			for(int j = 0; j < nClusters; j++) {
				if(simMaxLL > clusterLL[j]) {
					nExtreme[j] ++;
				}			
			}
		}
	}


	free(simCasInW);
	free(simConInW);
	free(simll);

	free(indAll);
	free(simCass);
	free(simCons);

	return(nExtreme);

}

int * monteCarlo4DSphe(double * x1, double * y1, double * x2, double * y2, int * locEnding, int locCount, int casCount, int allCount, double wSize, int wCount, int elimIntersectOD, int highLow, double * clusterLL, int nClusters, int nSim) {
	int * nExtreme;

	if(NULL == (nExtreme = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	for(int i = 0; i < nClusters; i++)
		nExtreme[i] = 0;

	int * indAll;
	int * simCass;
	int * simCons;

	if(NULL == (indAll = (int *) malloc (allCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simCass = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simCons = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	int indID, simCas, simCon;

	int * simCasInW;
	int * simConInW;
	double * simll;
	
	if(NULL == (simCasInW = (int *) malloc (locCount * wCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simConInW = (int *) malloc (locCount * wCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (simll = (double *) malloc (locCount * wCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}


	double simMaxLL;

	for(int i = 0; i < nSim; i++) {
		randomLabel(indAll, casCount, allCount);
	
		indID = 0;
		for(int j = 0; j < locCount; j++) {
			simCas = 0;
			simCon = 0;
			for(; indID < locEnding[j]; indID ++) {
				if(indAll[indID] == 1) {
					simCas ++;
				}
				else {
					simCon ++;
				}
			}
			simCass[j] = simCas;
			simCons[j] = simCon;
		}

		getCCCount4DSph(x1, y1, x2, y2, simCass, simCons, locCount, wSize, wCount, simCasInW, simConInW, elimIntersectOD);	
		loglikelihood(simll, simCasInW, simConInW, locCount * wCount, casCount, allCount - casCount, highLow);
		
		simMaxLL = 1;
		int k = 0;
		for(; k < locCount * wCount; k++) {
			if(simll[k] < 0) {
				simMaxLL = simll[k];
				k++;
				break;
			}
		}

		for(; k < locCount * wCount; k++) {
			if(simll[k] < 0 && simll[k] > simMaxLL) {
				simMaxLL = simll[k];
			}
		}

//		printf("%d:\t%lf\n", i, simMaxLL);

		if(simMaxLL < 0) {
			for(int j = 0; j < nClusters; j++) {
				if(simMaxLL > clusterLL[j]) {
					nExtreme[j] ++;
				}			
			}
		}
	}

	free(simCasInW);
	free(simConInW);
	free(simll);

	free(indAll);
	free(simCass);
	free(simCons);

	return(nExtreme);
}


