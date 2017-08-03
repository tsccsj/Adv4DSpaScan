/**
 * Adv4DScanBer.c
 * Authors: Yizhao Gao <yizhaotsccsj@gmail.com>
 * Date: {08/01/2017}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "io.h"
#include "scan.h"

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
	if(windowShape == 2) {
		totalWindowsL = wCount * wCount;
	}
	else {
		totalWindowsL = wCount;
	}
	if(pAsCenter != 2) {
		totalWindowsL = totalWindowsL * locCount;
	}
	else {
		totalWindowsL = totalWindowsL * nRow * nCol * nRow * nCol;
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
	if(NULL == (ll = (double *) malloc (locCount * wCount * sizeof(double)))) {
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
	



	free(casInW);
	free(conInW);
	free(ll);

	free(x1);
	free(y1);
	free(x2);
	free(y2);
	free(nCass);
	free(nCons);

	return 0;

}
