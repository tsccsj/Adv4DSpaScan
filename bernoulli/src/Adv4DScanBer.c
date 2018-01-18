/**
 * Adv4DScanBer.c
 * Authors: Yizhao Gao <yizhaotsccsj@gmail.com>
 * Date: {08/06/2017}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "io.h"
#include "scan.h"
#include "mc.h"

int main(int argc, char ** argv) {

	if(argc != 10 && argc != 15) {
		
		printf("Incorrect of parameters: \n\tAdv4DScanBer [inputFile] [scanningWindowShape] [windowInc] [windowCount] [#ofClusters] [#ofMonteCarlo] [HighOrLowIndicator] [elimIntersectOD] [CenterLocation] ([cellSize] [xMin] [yMin] [xMax] [yMax])\n");
		printf("[scanningWindowShape]: the shape of scanning windows \n\t0: 4D Circular \n\t1: 2 2D cricles with the same radius\n\t2: 2 2D circles with different radius\n");
		printf("[HighOrLowIndicator]: high or low clusters to detect\n\t1: High Only\n\t-1: Low Only\n\t0: Both\n");
		printf("[elimIntersectOD]: can O and D intersect\n\t-1: Allow intersect\n\t1: No intersecting\n");
		printf("[CenterLocation]: cluster centers\n\t1: input points\n\t2: regular grids");
		printf("[cellSize]: the cell size when regular grids are used as cluster centers\n");
		printf("[xMin] [yMin] [xMax] [yMax]: the spatial extent of grids\n");
		exit(1);
	}

	int windowShape = atoi(argv[2]);
	double wSize = atof(argv[3]);
        int wCount = atoi(argv[4]);
        int nClusters = atoi(argv[5]);
        int nSim = atoi(argv[6]);
	int highLow = atoi(argv[7]);
	int elimIntersectOD = atoi(argv[8]);
	int pAsCenter = atoi(argv[9]);
	
//variables for grids as centers
	double cellSize;
	double xMin;
	double yMin;
	double xMax;
	double yMax;
	int nRow;
	int nCol;

	if(pAsCenter == 2) {
	
		cellSize = atof(argv[10]);
		xMin = atof(argv[11]);
		yMin = atof(argv[12]);
		xMax = atof(argv[13]);
		yMax = atof(argv[14]);
	
		nCol = ceil((xMax - xMin) / cellSize);
		nRow = ceil((yMax - yMin) / cellSize);

		xMax = xMin + cellSize * nCol;
		yMax = yMin + cellSize * nRow;

		nCol = nCol + 1;
		nRow = nRow + 1;

		printf("Window centers: %d rows and %d cols.\n", nRow, nCol);
	}

	double * x1;
	double * y1;
	double * x2;
	double * y2;
	int * nCass;
	int * nCons;

	int casCount;
	int conCount;
	int locCount;

	FILE * file;
	if(NULL == (file = fopen(argv[1], "r"))) {
		printf("ERROR: Can't open the input file.\n");
		exit(1);
	}

	locCount = getNumPoints(file);
	if(NULL == (x1 = (double *) malloc (locCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (y1 = (double *) malloc (locCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (x2 = (double *) malloc (locCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (y2 = (double *) malloc (locCount * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (nCass = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (nCons = (int *) malloc (locCount * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	readFile(file, x1, y1, x2, y2, nCass, nCons, casCount, conCount);
	
	printf("There is %d locations\n", locCount);
	printf("Total CASE count: \t%d\nTotal CONTROL count:\t%d\n", casCount, conCount);

	fclose(file);
	printf("Finish reading input files.\n");

	long long totalWindowsL;
	int totalCenters;
	int windowPerCen;
	if(windowShape == 2) {
		totalWindowsL = wCount * wCount;
		windowPerCen = wCount * wCount;
	}
	else {
		totalWindowsL = wCount;
		windowPerCen = wCount;
	}
	if(pAsCenter != 2) {
		totalWindowsL = totalWindowsL * locCount;
		totalCenters = locCount;
	}
	else {
		totalWindowsL = totalWindowsL * nRow * nCol * nRow * nCol;
		totalCenters = nRow * nCol * nRow * nCol;
	}

	printf("There are %ld total scanning windows\n", totalWindowsL);

	if(totalWindowsL > INT_MAX) {
		printf("The number of scanning windows is too much for this program\n");
		exit(1);
	}

	int totalWindows = totalWindowsL;

	int * casInW;
	int * conInW;
	double * ll;

	if(NULL == (casInW = (int *) malloc (totalWindows * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (conInW = (int *) malloc (totalWindows * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (ll = (double *) malloc (totalWindows * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	if(windowShape == 0) {
		if(pAsCenter == 2) {
			getCCCount4DSph(x1, y1, x2, y2, nCass, nCons, locCount, wSize, wCount, xMin, yMin, cellSize, nCol, nRow, casInW, conInW, elimIntersectOD);
		}
		else {
			getCCCount4DSph(x1, y1, x2, y2, nCass, nCons, locCount, wSize, wCount, casInW, conInW, elimIntersectOD);
		}
	}
	else if(windowShape == 1) {
		if(pAsCenter == 2) {
			getCCCount2DSame(x1, y1, x2, y2, nCass, nCons, locCount, wSize, wCount, xMin, yMin, cellSize, nCol, nRow, casInW, conInW, elimIntersectOD);
		}
		else {
			getCCCount2DSame(x1, y1, x2, y2, nCass, nCons, locCount, wSize, wCount, casInW, conInW, elimIntersectOD);
		}
	}
	else {
		if(pAsCenter == 2) {
			getCCCount2DDiff(x1, y1, x2, y2, nCass, nCons, locCount, wSize, wCount, xMin, yMin, cellSize, nCol, nRow, casInW, conInW, elimIntersectOD);
		}
		else {
			getCCCount2DDiff(x1, y1, x2, y2, nCass, nCons, locCount, wSize, wCount, casInW, conInW, elimIntersectOD);
		}
	
	}

//	printf("Checkpoint 1\n");
	
	loglikelihood(ll, casInW, conInW, totalWindows, casCount, conCount, highLow);
	
//	printf("Checkpoint 2\n");

	int * center;
	int * radiusO;
	int * radiusD;
	double * cLL;

	if(NULL == (center = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	if(NULL == (radiusO = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	if(NULL == (radiusD = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	if(NULL == (cLL = (double *) malloc (nClusters * sizeof(double)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}

	//printf("windowShape = %d, pAsCenter = %d\n", windowShape, pAsCenter);
	if(windowShape == 0) {
		if(pAsCenter == 2) {
			findTopNCluster4DSph(xMin, yMin, cellSize, nRow, nCol, ll, wSize, wCount, center, radiusO, cLL, nClusters);
		}
		else {
			findTopNCluster4DSph(x1, y1, x2, y2, locCount, ll, wSize, wCount, center, radiusO,  cLL, nClusters);
		}
	}
	else if(windowShape == 1) {
		if(pAsCenter == 2) {
			findTopNCluster2DSame(xMin, yMin, cellSize, nRow, nCol, ll, wSize, wCount, center, radiusO, cLL, nClusters);
		}
		else {
			findTopNCluster2DSame(x1, y1, x2, y2, locCount, ll, wSize, wCount, center, radiusO,  cLL, nClusters);
		}
	}
	else {
		if(pAsCenter == 2) {
			findTopNCluster2DDiff(xMin, yMin, cellSize, nRow, nCol, ll, wSize, wCount, center, radiusO, radiusD, cLL, nClusters);
		}
		else {
			findTopNCluster2DDiff(x1, y1, x2, y2, locCount, ll, wSize, wCount, center, radiusO, radiusD, cLL, nClusters);
		}
	}

//	printf("Checkpoint 3\n");

	int * clusterCas;
	int * clusterCon;
	double * cRadiusO;
	double * cRadiusD;
	int * highCluster;

	if(NULL == (clusterCas = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (clusterCon = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (cRadiusO = (double *) malloc (nClusters * sizeof(double)))) {	
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (cRadiusD = (double *) malloc (nClusters * sizeof(double)))) {	
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}
	if(NULL == (highCluster = (int *) malloc (nClusters * sizeof(int)))) {
		printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
		exit(1);
	}


	for(int i = 0; i < nClusters; i ++) {
		if(center[i] == -1) {
			nClusters = i;
			break;
		}
		if(windowShape == 2) {
			clusterCas[i] = casInW[center[i] * windowPerCen + radiusO[i] * wCount + radiusD[i]];
			clusterCon[i] = conInW[center[i] * windowPerCen + radiusO[i] * wCount + radiusD[i]];
			cRadiusO[i] = wSize * (radiusO[i] + 1);
			cRadiusD[i] = wSize * (radiusD[i] + 1);
			//printf("%d, %d, %d, %lf\n", i, radiusO[i], radiusD[i], wSize);
			//printf("%lf\n", wSize);
		}
		else {
			clusterCas[i] = casInW[center[i] * wCount + radiusO[i]];
			clusterCon[i] = conInW[center[i] * wCount + radiusO[i]];
			cRadiusO[i] = wSize * (radiusO[i] + 1);
		}

		double expCas = (double)casCount*(clusterCas[i] + clusterCon[i])/(casCount+conCount);
		if(clusterCas[i] > expCas)      //High clusters
			highCluster[i] = 1;
		else
			highCluster[i] = 0;
	}

//Here is the Monte Carlo Simulation
	int * nExtreme;

	if(nSim > 0) {
	
		int * locEnding;
		int accCount = 0;
		if(NULL == (locEnding = (int *) malloc (locCount * sizeof(int)))) {
			printf("ERROR: Out of memory at line %d in file %s\n", __LINE__, __FILE__);
			exit(1);
		}

		for(int i = 0; i < locCount; i ++) {
			accCount += nCass[i] + nCons[i];
			locEnding[i] = accCount;
		}

		nExtreme = monteCarlo(x1, y1, x2, y2, locEnding, locCount, casCount, casCount + conCount, windowShape, wSize, wCount, elimIntersectOD, highLow, cLL, nClusters, nSim, pAsCenter, xMin, yMin, cellSize, nRow, nCol);
		
		free(locEnding);
	}

	printf("############### Cluster Info ###############\n");
	printf("ID,HL,X1,Y1,X2,Y2,");
	if(windowShape == 2) {
		printf("RadiusO,RaiudsD,#Cas,#Con,Exp#Cas,Exp#Con");
	}
	else {
		printf("Radius,#Cas,#Con,Exp#Cas,Exp#Con");
	}
	if(nSim > 0)
		printf(",LL,P\n");
	else
		printf(",LL\n");
	for(int i = 0; i < nClusters; i ++) {

		double expCas = (double)casCount*(clusterCas[i] + clusterCon[i])/(casCount+conCount);
		double expCon = (double)conCount*(clusterCas[i] + clusterCon[i])/(casCount+conCount);

		printf("%d", i);
		if(highCluster[i]) 	//High clusters
			printf(",H");
		else			//Low clusters
			printf(",L");

		if(pAsCenter == 2) {
			int temp = center[i];
			int cX2 = temp % nCol;
			temp = temp / nCol;
			int cY2 = temp % nRow;
			temp = temp / nRow;
			int cX1 = temp % nCol;
			int cY1 = temp / nCol;
			
			printf(",%lf,%lf,%lf,%lf,", xMin + cellSize * cX1, yMin + cellSize * cY1, xMin + cellSize * cX2, yMin + cellSize * cY2);
		}
		else {
			printf(",%lf,%lf,%lf,%lf,", x1[center[i]], y1[center[i]], x2[center[i]], y2[center[i]]);
		}
		if(windowShape == 2) {
			printf("%lf,%lf", cRadiusO[i], cRadiusD[i]);
		}
		else {
			printf("%lf", cRadiusO[i]);
		}
		printf(",%d,%d,%lf,%lf", clusterCas[i], clusterCon[i], expCas, expCon);
		if(nSim > 0)
			printf(",%lf,%lf\n", cLL[i], (double)(nExtreme[i] + 1) / (nSim + 1));
		else
			printf(",%lf\n", cLL[i]);
		
	}

	printf("############ Cluster Membership ############\n");

	int inCluster;
	double distanceO, distanceD;		
	double radiusValue;

	for(int i = 0; i < locCount; i++) {
		inCluster = 0;

		for(int j = 0; j < nClusters; j++) {

			if(windowShape == 0) {
				distanceO = sqrt((x1[i] - x1[center[j]]) * (x1[i] - x1[center[j]]) + (y1[i] - y1[center[j]]) * (y1[i] - y1[center[j]]) + (x2[i] - x2[center[j]]) * (x2[i] - x2[center[j]]) + (y2[i] - y2[center[j]]) * (y2[i] - y2[center[j]]));
				if(distanceO <= cRadiusO[j]) {
					inCluster = 1;
				}
			}
			else if(windowShape == 1) {
				distanceO = sqrt((x1[i] - x1[center[j]]) * (x1[i] - x1[center[j]]) + (y1[i] - y1[center[j]]) * (y1[i] - y1[center[j]]));
				distanceD = sqrt((x2[i] - x2[center[j]]) * (x2[i] - x2[center[j]]) + (y2[i] - y2[center[j]]) * (y2[i] - y2[center[j]]));
				if(distanceO <= cRadiusO[j] && distanceD <= cRadiusO[j]) {
					inCluster = 1;
				}		
			}
			else {
				distanceO = sqrt((x1[i] - x1[center[j]]) * (x1[i] - x1[center[j]]) + (y1[i] - y1[center[j]]) * (y1[i] - y1[center[j]]));
				distanceD = sqrt((x2[i] - x2[center[j]]) * (x2[i] - x2[center[j]]) + (y2[i] - y2[center[j]]) * (y2[i] - y2[center[j]]));
				if(distanceO <= cRadiusO[j] && distanceD <= cRadiusD[j]) {
					inCluster = 1;
				}		
			}

			if(inCluster == 1) {
				printf("%lf,%lf,%lf,%lf,%d,%d,%d\n", x1[i], y1[i], x2[i], y2[i], nCass[i], nCons[i], j);
				break;
			}
		}

		if(inCluster == 0)
			printf("%lf,%lf,%lf,%lf,%d,%d,-1\n", x1[i], y1[i], x2[i], y2[i], nCass[i], nCons[i]);
	}

//	printf("Checkpoint 4\n");


	free(clusterCas);
	free(clusterCon);
	free(cRadiusO);
	free(cRadiusD);
	free(highCluster);

	free(center);
	free(radiusO);
	free(radiusD);
	free(cLL);


	free(casInW);
	free(conInW); 
	free(ll);

	free(x1);
	free(y1);
	free(x2);
	free(y2);
	free(nCass);
	free(nCons);

	if(nSim > 0)
		free(nExtreme);

	return 0;

}
