/**
 * scan.c
 * Authors: Yizhao Gao <yizhaotsccsj@gmail.com>
 * Date: {08/03/2017}
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void getCCCount4DSph(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, int * casInW, int * conInW, int elimIntersectOD) {
	double distance;
	int minWindow;

	for(int i = 0; i < locCount * wCount; i++) {
		casInW[i] = 0;
		conInW[i] = 0;
	}

#pragma omp parallel for private(distance, minWindow)
	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < locCount; j++) {
			distance = sqrt((x1[i] - x1[j]) * (x1[i] - x1[j]) + (y1[i] - y1[j]) * (y1[i] - y1[j]) + (x2[i] - x2[j]) * (x2[i] - x2[j]) + (y2[i] - y2[j]) * (y2[i] - y2[j]));
			minWindow = (int)(ceil(distance / wSize));
			if(minWindow > 0)
				minWindow --;
			for(int k = minWindow; k < wCount; k++) {
				casInW[i * wCount + k] += nCass[j];
				conInW[i * wCount + k] += nCons[j];
			}	
		}

		if(elimIntersectOD > 0) {
			double ODDistance = sqrt((x1[i] - x2[i]) * (x1[i] - x2[i]) + (y1[i] - y2[i]) * (y1[i] - y2[i])) / 2;
			int maxWindow = ODDistance / wSize;
			for(int k = maxWindow; k < wCount; k++) {
				casInW[i * wCount + k] = -1;
			}
		}
	}

	return;
}

void getCCCount4DSph(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, double xMin, double yMin, double cellSize, int nCol, int nRow, int * casInW, int * conInW, int elimIntersectOD) {

	double distance;
	int minWindow;

	int totalWindows = nCol * nCol * nRow * nRow * wCount;
	
	for(int i = 0; i < totalWindows; i++) {
		casInW[i] = 0;
		conInW[i] = 0;
	}

	int totalCenters = nRow * nRow * nCol * nCol;

#pragma omp parallel for private(distance, minWindow)
	for(int i = 0; i < totalCenters; i++) {

		int temp;	
		int cX2 = i % nCol;
		temp = i / nCol;
		int cY2 = temp % nRow;
		temp = temp / nRow;
		int cX1 = temp % nCol;
		int cY1 = temp / nCol;

		double centerX1 = xMin + cellSize * cX1;
		double centerY1 = yMin + cellSize * cY1;
		double centerX2 = xMin + cellSize * cX2;
		double centerY2 = yMin + cellSize * cY2;
		
		for(int j = 0; j < locCount; j++) {
		
			distance = sqrt((centerX1 - x1[j]) * (centerX1 - x1[j]) + (centerY1 - y1[j]) * (centerY1 - y1[j]) + (centerX2 - x2[j]) * (centerX2 - x2[j]) + (centerY2 - y2[j]) * (centerY2 - y2[j]));
			minWindow = (int)(ceil(distance / wSize));
			
			if(minWindow > 0)
				minWindow --;
			for(int k = minWindow; k < wCount; k++) {
				casInW[i * wCount + k] += nCass[j];
				conInW[i * wCount + k] += nCons[j];
			}	
		}

		if(elimIntersectOD > 0) {
			double ODDistance = sqrt((centerX1 - centerX2) * (centerX1 - centerX2) + (centerY1 - centerY2) * (centerY1 - centerY2)) / 2;
			int maxWindow = ODDistance / wSize;
			for(int k = maxWindow; k < wCount; k++) {
				casInW[i * wCount + k] = -1;
			}
		}
	}		
}

void getCCCount2DSame(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, int * casInW, int * conInW, int elimIntersectOD) {

	double distanceO, distanceD;
	int minWindow;

	for(int i = 0; i < locCount * wCount; i++) {
		casInW[i] = 0;
		conInW[i] = 0;
	}

#pragma omp parallel for private(distanceO, distanceD, minWindow)
	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < locCount; j++) {
			distanceO = sqrt((x1[i] - x1[j]) * (x1[i] - x1[j]) + (y1[i] - y1[j]) * (y1[i] - y1[j]));
			distanceD = sqrt((x2[i] - x2[j]) * (x2[i] - x2[j]) + (y2[i] - y2[j]) * (y2[i] - y2[j]));
			if(distanceO > distanceD) {
				minWindow = (int)(ceil(distanceO / wSize));
			}
			else {
				minWindow = (int)(ceil(distanceD / wSize));
			}

			if(minWindow > 0)
				minWindow --;
			for(int k = minWindow; k < wCount; k++) {
				casInW[i * wCount + k] += nCass[j];
				conInW[i * wCount + k] += nCons[j];
			}	
		}

		if(elimIntersectOD > 0) {
			double ODDistance = sqrt((x1[i] - x2[i]) * (x1[i] - x2[i]) + (y1[i] - y2[i]) * (y1[i] - y2[i])) / 2;
			int maxWindow = ODDistance / wSize;
			for(int k = maxWindow; k < wCount; k++) {
				casInW[i * wCount + k] = -1;
			}
		}
	}

	return;
}

void getCCCount2DSame(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, double xMin, double yMin, double cellSize, int nCol, int nRow, int * casInW, int * conInW, int elimIntersectOD) {

	double distanceO, distanceD;
	int minWindow;

	int totalWindows = nCol * nCol * nRow * nRow * wCount;
	
	for(int i = 0; i < totalWindows; i++) {
		casInW[i] = 0;
		conInW[i] = 0;
	}

	int totalCenters = nRow * nRow * nCol * nCol;

#pragma omp parallel for private(distanceO, distanceD, minWindow)
	for(int i = 0; i < totalCenters; i++) {

		int temp;	
		int cX2 = i % nCol;
		temp = i / nCol;
		int cY2 = temp % nRow;
		temp = temp / nRow;
		int cX1 = temp % nCol;
		int cY1 = temp / nCol;

		double centerX1 = xMin + cellSize * cX1;
		double centerY1 = yMin + cellSize * cY1;
		double centerX2 = xMin + cellSize * cX2;
		double centerY2 = yMin + cellSize * cY2;
		
		for(int j = 0; j < locCount; j++) {
		
			distanceO = sqrt((centerX1 - x1[j]) * (centerX1 - x1[j]) + (centerY1 - y1[j]) * (centerY1 - y1[j]));
			distanceD = sqrt((centerX2 - x2[j]) * (centerX2 - x2[j]) + (centerY2 - y2[j]) * (centerY2 - y2[j]));
			if(distanceO > distanceD) {
				minWindow = (int)(ceil(distanceO / wSize));
			}
			else {
				minWindow = (int)(ceil(distanceD / wSize));
			}
			
			if(minWindow > 0)
				minWindow --;
			for(int k = minWindow; k < wCount; k++) {
				casInW[i * wCount + k] += nCass[j];
				conInW[i * wCount + k] += nCons[j];
			}	
		}

		if(elimIntersectOD > 0) {
			double ODDistance = sqrt((centerX1 - centerX2) * (centerX1 - centerX2) + (centerY1 - centerY2) * (centerY1 - centerY2)) / 2;
			int maxWindow = ODDistance / wSize;
			for(int k = maxWindow; k < wCount; k++) {
				casInW[i * wCount + k] = -1;
			}
		}
	}		
}

void getCCCount2DDiff(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, int * casInW, int * conInW, int elimIntersectOD) {

	double distanceO, distanceD;
	int minWindowO, minWindowD;

	int windowPerCen = wCount * wCount;
	int totalWindows = windowPerCen * locCount; 

	for(int i = 0; i < totalWindows; i++) {
		casInW[i] = 0;
		conInW[i] = 0;
	}

#pragma omp parallel for private(distanceO, distanceD, minWindowO, minWindowD)
	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < locCount; j++) {
			distanceO = sqrt((x1[i] - x1[j]) * (x1[i] - x1[j]) + (y1[i] - y1[j]) * (y1[i] - y1[j]));
			distanceD = sqrt((x2[i] - x2[j]) * (x2[i] - x2[j]) + (y2[i] - y2[j]) * (y2[i] - y2[j]));
			minWindowO = (int)(ceil(distanceO / wSize));
			minWindowD = (int)(ceil(distanceD / wSize));

			if(minWindowO > 0)
				minWindowO --;
			if(minWindowD > 0)
				minWindowD --;

			for(int k = minWindowO; k < wCount; k++) {
				for(int l = minWindowD; l < wCount; l++) {

					casInW[i * windowPerCen + k * wCount + l] += nCass[j];
					conInW[i * windowPerCen + k * wCount + l] += nCons[j];
				}
			}	
		}

		if(elimIntersectOD > 0) {
			double ODDistance = sqrt((x1[i] - x2[i]) * (x1[i] - x2[i]) + (y1[i] - y2[i]) * (y1[i] - y2[i]));
			int sumWindow = (int)(ceil(ODDistance / wSize)) - 2;
			for(int k = 0; k < wCount; k++) {
				for(int j = 0; j < wCount; j++) {
					if(k + j > sumWindow) {
						casInW[i * windowPerCen + k * wCount + j] = -1;
					}
				}
			}
		}
	}

	return;
}

void getCCCount2DDiff(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, double xMin, double yMin, double cellSize, int nCol, int nRow, int * casInW, int * conInW, int elimIntersectOD) {

	double distanceO, distanceD;
	int minWindowO, minWindowD;

	int windowPerCen = wCount * wCount;
	int totalCenters = nRow * nRow * nCol * nCol;

	int totalWindows = totalCenters * windowPerCen;
	
	for(int i = 0; i < totalWindows; i++) {
		casInW[i] = 0;
		conInW[i] = 0;
	}

#pragma omp parallel for private(distanceO, distanceD, minWindowO, minWindowD)
	for(int i = 0; i < totalCenters; i++) {

		int temp;	
		int cX2 = i % nCol;
		temp = i / nCol;
		int cY2 = temp % nRow;
		temp = temp / nRow;
		int cX1 = temp % nCol;
		int cY1 = temp / nCol;

		double centerX1 = xMin + cellSize * cX1;
		double centerY1 = yMin + cellSize * cY1;
		double centerX2 = xMin + cellSize * cX2;
		double centerY2 = yMin + cellSize * cY2;
		
		for(int j = 0; j < locCount; j++) {
		
			distanceO = sqrt((centerX1 - x1[j]) * (centerX1 - x1[j]) + (centerY1 - y1[j]) * (centerY1 - y1[j]));
			distanceD = sqrt((centerX2 - x2[j]) * (centerX2 - x2[j]) + (centerY2 - y2[j]) * (centerY2 - y2[j]));
			minWindowO = (int)(ceil(distanceO / wSize));
			minWindowD = (int)(ceil(distanceD / wSize));

			if(minWindowO > 0)
				minWindowO --;
			if(minWindowD > 0)
				minWindowD --;

			for(int k = minWindowO; k < wCount; k++) {
				for(int l = minWindowD; l < wCount; l++) {

					casInW[i * windowPerCen + k * wCount + j] += nCass[j];
					conInW[i * windowPerCen + k * wCount + j] += nCons[j];
				}
			}
		}

		if(elimIntersectOD > 0) {
			double ODDistance = sqrt((centerX1 - centerX2) * (centerX1 - centerX2) + (centerY1 - centerY2) * (centerY1 - centerY2));
			int sumWindow = (int)(ceil(ODDistance / wSize)) - 2;
			for(int k = 0; k < wCount; k++) {
				for(int j = 0; j < wCount; j++) {
					if(k + j > sumWindow) {
						casInW[i * windowPerCen + k * wCount + j] = -1;
					}
				}
			}
		}
	}		
}



void loglikelihood(double * ll, int * casInW, int * conInW, int totalWindow, int casCount, int conCount, int highLow) {
	double cas, con, tot;
	double llTemp;
	int totCount = casCount + conCount;
	bool highCluster = true;
	bool lowCluster = true;
	if(highLow == 1)
		lowCluster = false;
	else if(highLow == -1)
		highCluster = false;

#pragma omp parallel for private(cas, con, tot, llTemp)
	for(int i = 0; i < totalWindow; i++) {
		cas = casInW[i];
		con = conInW[i];
		tot = cas + con;

		if(cas == -1) {
			ll[i] = 1;
		}
		else if(cas * conCount > con * casCount) { //High cluster of cases
			if(highCluster) {
				llTemp = cas * log(cas/tot);
				if(con > 0)
					llTemp += con * log(con/tot);
				if(casCount > cas)
					llTemp += (casCount - cas) * log((casCount - cas)/(totCount - tot));
				if(conCount > con)
					llTemp += (conCount - con) * log((conCount - con)/(totCount - tot));
				ll[i] = llTemp;
			}
			else
				ll[i] = 1;
		}
		else { //Low cluster of cases
			if(lowCluster) {
				llTemp = con * log(con/tot);
				if(cas > 0)
					llTemp += cas * log(cas/tot);
				if(casCount > cas)
					llTemp += (casCount - cas) * log((casCount - cas)/(totCount - tot));
				if(conCount > con)
					llTemp += (conCount - con) * log((conCount - con)/(totCount - tot));
				ll[i] = llTemp;
			}
			else
				ll[i] = 1;
		}
	}
	
}

void findTopNCluster4DSph(double * x1, double * y1, double * x2, double * y2, int locCount, double * ll, double wSize, int wCount, int * center, int * radius,  double * cLL, int nClusters) {
	if(nClusters < 1)
		return;

	int aCenter = -1;
	int aRadius = -1;

	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < wCount; j++) {
			if(ll[i * wCount + j] < 0) {
				if(aCenter < 0) {
					aCenter = i;
					aRadius = j;
				}
				else if(ll[i * wCount + j] > ll[aCenter * wCount + aRadius]) {
					aCenter = i;
					aRadius = j;
				}
			}
		}
	}

	center[0] = aCenter;
	radius[0] = aRadius;	
	cLL[0] = ll[aCenter * wCount + aRadius];

	double lastX1, lastY1, lastX2, lastY2, lastRad;
	lastX1 = x1[aCenter];
	lastY1 = y1[aCenter];
	lastX2 = x2[aCenter];
	lastY2 = y2[aCenter];

	lastRad = (aRadius + 1) * wSize;

	double distance;
	int maxWindow;

	for(int c = 1; c < nClusters; c ++) {
		//Remove intersecting clusters
		for(int i = 0; i < locCount; i++) {
			distance = sqrt((x1[i] - lastX1) * (x1[i] - lastX1) + (y1[i] - lastY1) * (y1[i] - lastY1) + (x2[i] - lastX2) * (x2[i] - lastX2) + (y2[i] - lastY2) * (y2[i] - lastY2)) - lastRad;
			maxWindow = ceil(distance / wSize) - 1;
			if(maxWindow < 0)
				maxWindow = 0;
			for(int j = maxWindow; j < wCount; j++) {
				ll[i * wCount + j] = 1;
			}			
		}
		

		//Find secoundary clusters
		aCenter = -1;
		aRadius = -1;

		for(int i = 0; i < locCount; i++) {
			for(int j = 0; j < wCount; j++) {
				if(ll[i * wCount + j] < 0) {
					if(aCenter < 0) {
						aCenter = i;
						aRadius = j;
					}
					else if(ll[i * wCount + j] > ll[aCenter * wCount + aRadius]) {
						aCenter = i;
						aRadius = j;
					}
				}
			}
		}
		center[c] = aCenter;
		radius[c] = aRadius;
		if(aCenter != -1) {
			cLL[c] = ll[aCenter * wCount + aRadius];
		}
		else {
			break;
		}

		lastX1 = x1[aCenter];
		lastY1 = y1[aCenter];
		lastX2 = x2[aCenter];
		lastY2 = y2[aCenter];

		lastRad = (aRadius + 1) * wSize;
	}

	return;	
}

void findTopNCluster4DSph(double xMin, double yMin, double cellSize, int nRow, int nCol, double * ll, double wSize, int wCount, int * center, int * radius,  double * cLL, int nClusters) {
	if(nClusters < 1)
		return;

	int totalCenters = nRow * nRow * nCol * nCol;
	
	int aCenter = -1;
	int aRadius = -1;

	for(int i = 0; i < totalCenters; i++) {
		for(int j = 0; j < wCount; j++) {
			if(ll[i * wCount + j] < 0) {
				if(aCenter < 0) {
					aCenter = i;
					aRadius = j;
				}
				else if(ll[i * wCount + j] > ll[aCenter * wCount + aRadius]) {
					aCenter = i;
					aRadius = j;
				}
			}
		}
	}

	center[0] = aCenter;
	radius[0] = aRadius;	
	cLL[0] = ll[aCenter * wCount + aRadius];

	int temp;	
	int cX2 = aCenter % nCol;
	temp = aCenter / nCol;
	int cY2 = temp % nRow;
	temp = temp / nRow;
	int cX1 = temp % nCol;
	int cY1 = temp / nCol;

	double lastX1 = xMin + cellSize * cX1;
	double lastY1 = yMin + cellSize * cY1;
	double lastX2 = xMin + cellSize * cX2;
	double lastY2 = yMin + cellSize * cY2;

	double lastRad = (aRadius + 1) * wSize;

	double distance;
	int maxWindow;

	double thisX1, thisY1, thisX2, thisY2;

	for(int c = 1; c < nClusters; c ++) {
		//Remove intersecting clusters
		for(int i = 0; i < totalCenters; i++) {
			cX2 = i % nCol;
			temp = i / nCol;
			cY2 = temp % nRow;
			temp = temp / nRow;
			cX1 = temp % nCol;
			cY1 = temp / nCol;

			thisX1 = xMin + cellSize * cX1;
			thisY1 = yMin + cellSize * cY1;
			thisX2 = xMin + cellSize * cX2;
			thisY2 = yMin + cellSize * cY2;	
		
			distance = sqrt((thisX1 - lastX1) * (thisX1 - lastX1) + (thisY1 - lastY1) * (thisY1 - lastY1) + (thisX2 - lastX2) * (thisX2 - lastX2) + (thisY2 - lastY2) * (thisY2 - lastY2)) - lastRad;
			maxWindow = ceil(distance / wSize) - 1;

			if(maxWindow < 0)
				maxWindow = 0;
			for(int j = maxWindow; j < wCount; j++) {
				ll[i * wCount + j] = 1;
			}			
		}
		

		//Find secoundary clusters
		aCenter = -1;
		aRadius = -1;

		for(int i = 0; i < totalCenters; i++) {
			for(int j = 0; j < wCount; j++) {
				if(ll[i * wCount + j] < 0) {
					if(aCenter < 0) {
						aCenter = i;
						aRadius = j;
					}
					else if(ll[i * wCount + j] > ll[aCenter * wCount + aRadius]) {
						aCenter = i;
						aRadius = j;
					}
				}
			}
		}
		center[c] = aCenter;
		radius[c] = aRadius;
		if(aCenter != -1) {
			cLL[c] = ll[aCenter * wCount + aRadius];
		}
		else {
			break;
		}	
		cX2 = aCenter % nCol;
		temp = aCenter / nCol;
		cY2 = temp % nRow;
		temp = temp / nRow;
		cX1 = temp % nCol;
		cY1 = temp / nCol;

		lastX1 = xMin + cellSize * cX1;
		lastY1 = yMin + cellSize * cY1;
		lastX2 = xMin + cellSize * cX2;
		lastY2 = yMin + cellSize * cY2;	
	
		lastRad = (aRadius + 1) * wSize;
	}

	return;	
}

void findTopNCluster2DSame(double * x1, double * y1, double * x2, double * y2, int locCount, double * ll, double wSize, int wCount, int * center, int * radius,  double * cLL, int nClusters) {
	if(nClusters < 1)
		return;

	int aCenter = -1;
	int aRadius = -1;

	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < wCount; j++) {
			if(ll[i * wCount + j] < 0) {
				if(aCenter < 0) {
					aCenter = i;
					aRadius = j;
				}
				else if(ll[i * wCount + j] > ll[aCenter * wCount + aRadius]) {
					aCenter = i;
					aRadius = j;
				}
			}
		}
	}

	center[0] = aCenter;
	radius[0] = aRadius;	
	cLL[0] = ll[aCenter * wCount + aRadius];

	double lastX1, lastY1, lastX2, lastY2, lastRad;
	lastX1 = x1[aCenter];
	lastY1 = y1[aCenter];
	lastX2 = x2[aCenter];
	lastY2 = y2[aCenter];

	lastRad = (aRadius + 1) * wSize;

	double distanceO, distanceD;
	int maxWindow;

	for(int c = 1; c < nClusters; c ++) {
		//Remove intersecting clusters
		for(int i = 0; i < locCount; i++) {
			distanceO = sqrt((x1[i] - lastX1) * (x1[i] - lastX1) + (y1[i] - lastY1) * (y1[i] - lastY1)) - lastRad;
			distanceD = sqrt((x2[i] - lastX2) * (x2[i] - lastX2) + (y2[i] - lastY2) * (y2[i] - lastY2)) - lastRad;

			if(distanceO > distanceD) {
				maxWindow = ceil(distanceO / wSize) - 1;
			}
			else {
				maxWindow = ceil(distanceD / wSize) - 1;
			}
			if(maxWindow < 0)
				maxWindow = 0;

			for(int j = maxWindow; j < wCount; j++) {
				ll[i * wCount + j] = 1;
			}			
		}
		

		//Find secoundary clusters
		aCenter = -1;
		aRadius = -1;

		for(int i = 0; i < locCount; i++) {
			for(int j = 0; j < wCount; j++) {
				if(ll[i * wCount + j] < 0) {
					if(aCenter < 0) {
						aCenter = i;
						aRadius = j;
					}
					else if(ll[i * wCount + j] > ll[aCenter * wCount + aRadius]) {
						aCenter = i;
						aRadius = j;
					}
				}
			}
		}
		center[c] = aCenter;
		radius[c] = aRadius;
		if(aCenter != -1) {
			cLL[c] = ll[aCenter * wCount + aRadius];
		}
		else {
			break;
		}

		lastX1 = x1[aCenter];
		lastY1 = y1[aCenter];
		lastX2 = x2[aCenter];
		lastY2 = y2[aCenter];

		lastRad = (aRadius + 1) * wSize;
	}

	return;	
}

void findTopNCluster2DSame(double xMin, double yMin, double cellSize, int nRow, int nCol, double * ll, double wSize, int wCount, int * center, int * radius,  double * cLL, int nClusters) {
	if(nClusters < 1)
		return;

	int totalCenters = nRow * nRow * nCol * nCol;
	
	int aCenter = -1;
	int aRadius = -1;

	for(int i = 0; i < totalCenters; i++) {
		for(int j = 0; j < wCount; j++) {
			if(ll[i * wCount + j] < 0) {
				if(aCenter < 0) {
					aCenter = i;
					aRadius = j;
				}
				else if(ll[i * wCount + j] > ll[aCenter * wCount + aRadius]) {
					aCenter = i;
					aRadius = j;
				}
			}
		}
	}

	center[0] = aCenter;
	radius[0] = aRadius;	
	cLL[0] = ll[aCenter * wCount + aRadius];

	int temp;	
	int cX2 = aCenter % nCol;
	temp = aCenter / nCol;
	int cY2 = temp % nRow;
	temp = temp / nRow;
	int cX1 = temp % nCol;
	int cY1 = temp / nCol;

	double lastX1 = xMin + cellSize * cX1;
	double lastY1 = yMin + cellSize * cY1;
	double lastX2 = xMin + cellSize * cX2;
	double lastY2 = yMin + cellSize * cY2;

	double lastRad = (aRadius + 1) * wSize;

	double distanceO, distanceD;
	int maxWindow;

	double thisX1, thisY1, thisX2, thisY2;

	for(int c = 1; c < nClusters; c ++) {
		//Remove intersecting clusters
		for(int i = 0; i < totalCenters; i++) {
			cX2 = i % nCol;
			temp = i / nCol;
			cY2 = temp % nRow;
			temp = temp / nRow;
			cX1 = temp % nCol;
			cY1 = temp / nCol;

			thisX1 = xMin + cellSize * cX1;
			thisY1 = yMin + cellSize * cY1;
			thisX2 = xMin + cellSize * cX2;
			thisY2 = yMin + cellSize * cY2;	
		
			distanceO = sqrt((thisX1 - lastX1) * (thisX1 - lastX1) + (thisY1 - lastY1) * (thisY1 - lastY1)) - lastRad;
			distanceD = sqrt((thisX2 - lastX2) * (thisX2 - lastX2) + (thisY2 - lastY2) * (thisY2 - lastY2)) - lastRad;

			if(distanceO > distanceD) {
				maxWindow = ceil(distanceO / wSize) - 1;
			}
			else {
				maxWindow = ceil(distanceD / wSize) - 1;
			}
			if(maxWindow < 0)
				maxWindow = 0;

			for(int j = maxWindow; j < wCount; j++) {
				ll[i * wCount + j] = 1;
			}			

		}
		

		//Find secoundary clusters
		aCenter = -1;
		aRadius = -1;

		for(int i = 0; i < totalCenters; i++) {
			for(int j = 0; j < wCount; j++) {
				if(ll[i * wCount + j] < 0) {
					if(aCenter < 0) {
						aCenter = i;
						aRadius = j;
					}
					else if(ll[i * wCount + j] > ll[aCenter * wCount + aRadius]) {
						aCenter = i;
						aRadius = j;
					}
				}
			}
		}
		center[c] = aCenter;
		radius[c] = aRadius;
		if(aCenter != -1) {
			cLL[c] = ll[aCenter * wCount + aRadius];
		}
		else {
			break;
		}	
		cX2 = aCenter % nCol;
		temp = aCenter / nCol;
		cY2 = temp % nRow;
		temp = temp / nRow;
		cX1 = temp % nCol;
		cY1 = temp / nCol;

		lastX1 = xMin + cellSize * cX1;
		lastY1 = yMin + cellSize * cY1;
		lastX2 = xMin + cellSize * cX2;
		lastY2 = yMin + cellSize * cY2;	
	
		lastRad = (aRadius + 1) * wSize;
	}

	return;	
}

void findTopNCluster2DDiff(double * x1, double * y1, double * x2, double * y2, int locCount, double * ll, double wSize, int wCount, int * center, int * radiusO, int * radiusD,  double * cLL, int nClusters) {
	if(nClusters < 1)
		return;

	int aCenter = -1;
	int aRadiusO = -1;
	int aRadiusD = -1;

	int windowPerCen = wCount * wCount;
	

	for(int i = 0; i < locCount; i++) {
		for(int j = 0; j < wCount; j++) {
			for(int k = 0; k < wCount; k++) {
				if(ll[i * windowPerCen + j * wCount + k] < 0) {
					if(aCenter < 0) {
						aCenter = i;
						aRadiusO = j;
						aRadiusD = k;
					}
					else if(ll[i * windowPerCen + j * wCount + k] > ll[aCenter * windowPerCen + aRadiusO * wCount + aRadiusD]) {
						aCenter = i;
						aRadiusO = j;
						aRadiusD = k;
					}
				}
			}	
		}
	}

	center[0] = aCenter;
	radiusO[0] = aRadiusO;
	radiusD[0] = aRadiusD;	
	cLL[0] = ll[aCenter * windowPerCen + aRadiusO * wCount + aRadiusD];

	double lastX1, lastY1, lastX2, lastY2, lastRadO, lastRadD;
	lastX1 = x1[aCenter];
	lastY1 = y1[aCenter];
	lastX2 = x2[aCenter];
	lastY2 = y2[aCenter];

	lastRadO = (aRadiusO + 1) * wSize;
	lastRadD = (aRadiusD + 1) * wSize;

	double distanceO, distanceD;
	int maxWindowO, maxWindowD;

	for(int c = 1; c < nClusters; c ++) {
		//Remove intersecting clusters
		for(int i = 0; i < locCount; i++) {
			distanceO = sqrt((x1[i] - lastX1) * (x1[i] - lastX1) + (y1[i] - lastY1) * (y1[i] - lastY1)) - lastRadO;
			distanceD = sqrt((x2[i] - lastX2) * (x2[i] - lastX2) + (y2[i] - lastY2) * (y2[i] - lastY2)) - lastRadD;

			maxWindowO = ceil(distanceO / wSize) - 1;
			maxWindowD = ceil(distanceD / wSize) - 1;
			
			if(maxWindowO < 0)
				maxWindowO = 0;
			if(maxWindowD < 0)
				maxWindowD = 0;

			for(int j = maxWindowO; j < wCount; j++) {
				for(int k = maxWindowD; k < wCount; k++) {
			
					ll[i * windowPerCen + j * wCount + k] = 1;
				}
			}			
		}
		

		//Find secoundary clusters
		aCenter = -1;
		aRadiusO = -1;
		aRadiusD = -1;

		for(int i = 0; i < locCount; i++) {
			for(int j = 0; j < wCount; j++) {
				for(int k = 0; k < wCount; k++) { 
					if(ll[i * windowPerCen + j * wCount + k] < 0) {
						if(aCenter < 0) {
							aCenter = i;
							aRadiusO = j;
							aRadiusD = k;
						}
						else if(ll[i * windowPerCen + j * wCount + k] > ll[aCenter * windowPerCen + aRadiusO * wCount + aRadiusD]) {
							aCenter = i;
							aRadiusO = j;
							aRadiusD = k;
						}
					}
				}
			}
		}

		center[c] = aCenter;
		radiusO[c] = aRadiusO;
		radiusD[c] = aRadiusD;	
		if(aCenter != -1) {
			cLL[c] = ll[aCenter * windowPerCen + aRadiusO * wCount + aRadiusD];
		}
		else {
			break;
		}

		lastX1 = x1[aCenter];
		lastY1 = y1[aCenter];
		lastX2 = x2[aCenter];
		lastY2 = y2[aCenter];

		lastRadO = (aRadiusO + 1) * wSize;
		lastRadD = (aRadiusD + 1) * wSize;
	}

	return;	
}

void findTopNCluster2DDiff(double xMin, double yMin, double cellSize, int nRow, int nCol, double * ll, double wSize, int wCount, int * center, int * radiusO, int * radiusD, double * cLL, int nClusters) {
	if(nClusters < 1)
		return;

	int totalCenters = nRow * nRow * nCol * nCol;
	int windowPerCen = wCount * wCount;
	
	int aCenter = -1;
	int aRadiusO = -1;
	int aRadiusD = -1;

	for(int i = 0; i < totalCenters; i++) {
		for(int j = 0; j < wCount; j++) {
			for(int k = 0; k < wCount; k++) {
				if(ll[i * windowPerCen + j * wCount + k] < 0) {
					if(aCenter < 0) {
						aCenter = i;
						aRadiusO = j;
						aRadiusD = k;
					}
					else if(ll[i * windowPerCen + j * wCount + k] > ll[aCenter * windowPerCen + aRadiusO * wCount + aRadiusD]) {
						aCenter = i;
						aRadiusO = j;
						aRadiusD = k;
					}
				}
			}
		}
	}

	center[0] = aCenter;
	radiusO[0] = aRadiusO;
	radiusD[0] = aRadiusD;	
	cLL[0] = ll[aCenter * windowPerCen + aRadiusO * wCount + aRadiusD];

	int temp;	
	int cX2 = aCenter % nCol;
	temp = aCenter / nCol;
	int cY2 = temp % nRow;
	temp = temp / nRow;
	int cX1 = temp % nCol;
	int cY1 = temp / nCol;

	double lastX1 = xMin + cellSize * cX1;
	double lastY1 = yMin + cellSize * cY1;
	double lastX2 = xMin + cellSize * cX2;
	double lastY2 = yMin + cellSize * cY2;

	double lastRadO = (aRadiusO + 1) * wSize;
	double lastRadD = (aRadiusD + 1) * wSize;

	double distanceO, distanceD;
	int maxWindowO, maxWindowD;;

	double thisX1, thisY1, thisX2, thisY2;

	for(int c = 1; c < nClusters; c ++) {
		//Remove intersecting clusters
		for(int i = 0; i < totalCenters; i++) {
			cX2 = i % nCol;
			temp = i / nCol;
			cY2 = temp % nRow;
			temp = temp / nRow;
			cX1 = temp % nCol;
			cY1 = temp / nCol;

			thisX1 = xMin + cellSize * cX1;
			thisY1 = yMin + cellSize * cY1;
			thisX2 = xMin + cellSize * cX2;
			thisY2 = yMin + cellSize * cY2;	
		
			distanceO = sqrt((thisX1 - lastX1) * (thisX1 - lastX1) + (thisY1 - lastY1) * (thisY1 - lastY1)) - lastRadO;
			distanceD = sqrt((thisX2 - lastX2) * (thisX2 - lastX2) + (thisY2 - lastY2) * (thisY2 - lastY2)) - lastRadD;
			maxWindowO = ceil(distanceO / wSize) - 1;
			maxWindowD = ceil(distanceD / wSize) - 1;
			
			if(maxWindowO < 0)
				maxWindowO = 0;
			if(maxWindowD < 0)
				maxWindowD = 0;

			for(int j = maxWindowO; j < wCount; j++) {
				for(int k = maxWindowD; k < wCount; k++) {
			
					ll[i * windowPerCen + j * wCount + k] = 1;
				}
			}	

		}
		

		//Find secoundary clusters
		aCenter = -1;
		aRadiusO = -1;
		aRadiusD = -1;

		for(int i = 0; i < totalCenters; i++) {
			for(int j = 0; j < wCount; j++) {
				for(int k = 0; k < wCount; k++) {
					if(ll[i * windowPerCen + j * wCount + k] < 0) {
						if(aCenter < 0) {
							aCenter = i;
							aRadiusO = j;
							aRadiusD = k;
						}
						else if(ll[i * windowPerCen + j * wCount + k] > ll[aCenter * windowPerCen + aRadiusO * wCount + aRadiusD]) {
							aCenter = i;
							aRadiusO = j;
							aRadiusD = k;
						}
					}
				}	
			}
		}
		center[c] = aCenter;
		radiusO[c] = aRadiusO;
		radiusD[c] = aRadiusD;	
		if(aCenter != -1) {
			cLL[c] = ll[aCenter * windowPerCen + aRadiusO * wCount + aRadiusD];
		}
		else {
			break;
		}

		cX2 = aCenter % nCol;
		temp = aCenter / nCol;
		cY2 = temp % nRow;
		temp = temp / nRow;
		cX1 = temp % nCol;
		cY1 = temp / nCol;

		lastX1 = xMin + cellSize * cX1;
		lastY1 = yMin + cellSize * cY1;
		lastX2 = xMin + cellSize * cX2;
		lastY2 = yMin + cellSize * cY2;	
	
		lastRadO = (aRadiusO + 1) * wSize;
		lastRadD = (aRadiusD + 1) * wSize;
	}

	return;	
}


