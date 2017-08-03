#ifndef SCANH
#define SCANH

void getCCCount4DSph(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, int * casInW, int * conInW, int elimIntersectOD);
void getCCCount4DSph(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, double xMin, double yMin, double cellSize, int nCol, int nRow, int * casInW, int * conInW, int elimIntersectOD);

void getCCCount2DSame(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, int * casInW, int * conInW, int elimIntersectOD);
void getCCCount2DSame(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, double xMin, double yMin, double cellSize, int nCol, int nRow, int * casInW, int * conInW, int elimIntersectOD);

void getCCCount2DDiff(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, int * casInW, int * conInW, int elimIntersectOD);
void getCCCount2DDiff(double * x1, double * y1, double * x2, double * y2, int * nCass, int * nCons, int locCount, double wSize, int wCount, double xMin, double yMin, double cellSize, int nCol, int nRow, int * casInW, int * conInW, int elimIntersectOD);

void loglikelihood(double * ll, int * casInW, int * conInW, int totalWindow, int casCount, int conCount, int highLow);


void findTopNCluster(double * x1, double * y1, double * x2, double * y2, int locCount, double * ll, double wSize, int wCount, int * center, int * radius,  double * cLL, int nClusters);

#endif
