#ifndef CG_CRS_H
#define CG_CRS_H

#ifdef SPARSELIB_EXPORTS
#define SPARSELIB_API __declspec(dllexport) 
#else
#define SPARSELIB_API __declspec(dllimport) 
#endif

#include "Evaluator.h"
#include "MatrixData.h"
#include "TimeCounter.h"

class SPARSELIB_API CRSFormat : public Evaluator
{
public:
	const double EPSS = 1E-5;
	double **A;
	int **Ind;
	double *a, *x0;
	int n, nnzz;

	double *AA;
	int *JA;
	int *IA;

	MatrixData matrixData;

	CRSFormat(MatrixData matrixData, double * a, double * x0);
	~CRSFormat();

	virtual double fillMatrix();
	void arrayCpy(double * x, double * y, const int n);
	void intArrayCpy(int * x, int * y, const int n);
	virtual double minimal();
	virtual void generateNewYs();
	virtual void beforeMinimal();

	inline double vecProd(double *x, double *y, const int n);
	void MatrixByVectorCRS(double *aa, int *ja, int *ia, int k, double *x, double* res);
};

#endif