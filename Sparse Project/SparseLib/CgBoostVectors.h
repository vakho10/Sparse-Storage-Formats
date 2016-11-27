#ifndef CG_BOOST_VECTORS
#define CG_BOOST_VECTORS

#ifdef SPARSELIB_EXPORTS
#define SPARSELIB_API __declspec(dllexport) 
#else
#define SPARSELIB_API __declspec(dllimport) 
#endif

#include "Evaluator.h"
#include "MatrixData.h"
#include "TimeCounter.h"

#include "boost\numeric\ublas\vector_sparse.hpp"
#include "boost\numeric\ublas\matrix_vector.hpp"
#include "boost\numeric\ublas\matrix_sparse.hpp"

namespace ublas = boost::numeric::ublas;

class SPARSELIB_API CgBoostVectors : public Evaluator
{
public:
	const double EPS = 1E-5;
	MatrixData matrixData;

	ublas::vector<long double> x0;
	ublas::vector<long double> a;
	ublas::compressed_matrix<double, ublas::row_major, 0> A;

	CgBoostVectors(MatrixData matrixData, double* a, double* x0);
	~CgBoostVectors();

	virtual double fillMatrix();
	virtual double minimal();
	virtual void generateNewYs();
	virtual void beforeMinimal();
};

#endif // CG_BOOST_VECTORS