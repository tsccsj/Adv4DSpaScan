#ifndef MCH
#define MCH

int * monteCarlo(double * x1, double * y1, double * x2, double * y2, int * locEnding, int locCount, int casCount, int allCount, int windowShape, double wSize, int wCount, int elimIntersectOD, int highLow, double * clusterLL, int nClusters, int nSim, int pAsCenter, double xMin, double yMin, double cellSize, int nRow, int nCol);
int * monteCarlo4DSphe(double * x1, double * y1, double * x2, double * y2, int * locEnding, int locCount, int casCount, int allCount, double wSize, int wCount, int elimIntersectOD, int highLow, double * clusterLL, int nClusters, int nSim);

#endif
