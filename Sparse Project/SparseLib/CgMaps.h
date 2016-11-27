#ifndef CG_MAPS_H
#define CG_MAPS_H

#ifdef SPARSELIB_EXPORTS
#define SPARSELIB_API __declspec(dllexport) 
#else
#define SPARSELIB_API __declspec(dllimport) 
#endif

#include <map>

#include "Evaluator.h"
#include "MatrixData.h"
#include "TimeCounter.h"

#include <boost\numeric\ublas\matrix_sparse.hpp>
#include <boost\numeric\ublas\vector_sparse.hpp>
#include <boost\numeric\ublas\operation.hpp>

using namespace boost::numeric::ublas;

class SPARSELIB_API CgMaps : public Evaluator
{
public:
	const double EPS = 1E-5;
	MatrixData matrixData;

	mapped_vector<double> a_MAP, x0_MAP;

	//std::map<int, double> ARR_MAP;
	mapped_matrix<double> mappedMatrix;

	CgMaps(MatrixData matrixData, double* a, double* x0);
	~CgMaps();
	inline	double vecProd1(double *x, double *y, const int n);
	inline	void vecToMap(double *r, double *x0, std::map<int, double> map, const int n);

	virtual double fillMatrix();
	virtual double minimal();
	virtual void generateNewYs();
	virtual void beforeMinimal();
};

#endif // CG_MAPS_H