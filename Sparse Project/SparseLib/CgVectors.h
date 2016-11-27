#ifndef CG_VECTORS_H
#define CG_VECTORS_H

#ifdef SPARSELIB_EXPORTS
#define SPARSELIB_API __declspec(dllexport) 
#else
#define SPARSELIB_API __declspec(dllimport) 
#endif

#include <vector>
#include <numeric> // for inner_product

#include "Evaluator.h"
#include "MatrixData.h"

class SPARSELIB_API CgVectors : public Evaluator
{
public:
	const double EPS = 1E-5;
	MatrixData matrixData;

	std::vector<long double> x00;
	std::vector<long double> row1;
	std::vector<long double> row2;
	std::vector<std::vector<long double>> AA;
	std::vector<std::vector<long double>> AA0;
	std::vector<long double> aa;

	// Constructor
	CgVectors(MatrixData matrixData, double* a, double* x0);
	~CgVectors();

	virtual double fillMatrix();
	virtual double minimal();
	virtual void generateNewYs();
	virtual void beforeMinimal();
};

#endif // CG_VECTORS_H