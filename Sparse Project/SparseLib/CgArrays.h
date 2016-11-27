#ifndef CG_ARRAY_H
#define CG_ARRAY_H

#ifdef SPARSELIB_EXPORTS
#define SPARSELIB_API __declspec(dllexport) 
#else
#define SPARSELIB_API __declspec(dllimport) 
#endif

#include "Evaluator.h"
#include "MatrixData.h"
#include "TimeCounter.h"

class SPARSELIB_API CgArrays : public Evaluator
{
public:
	const double EPS = 1E-5;
	MatrixData matrixData;

	double **A_ARR;
	double *a_ARR, *x0_ARR;

	CgArrays(MatrixData matrixData);
	~CgArrays();

	// Helper function
	inline	double vecProd(double *x, double *y, const int n);

	virtual double fillMatrix();
	virtual double minimal();
	virtual void generateNewYs();
	virtual void beforeMinimal();
};

#endif // CG_ARRAY_H