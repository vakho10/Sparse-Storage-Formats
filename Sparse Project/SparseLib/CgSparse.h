#ifndef CG_SPARSE_H
#define CG_SPARSE_H

#ifdef SPARSELIB_EXPORTS
#define SPARSELIB_API __declspec(dllexport) 
#else
#define SPARSELIB_API __declspec(dllimport) 
#endif

#include "Evaluator.h"
#include "MatrixData.h"
#include "TimeCounter.h"

// CG_Sparse
class SPARSELIB_API CgSparse : public Evaluator
{
public:
	const double EPSS = 1E-5;
	double **A;
	int **Ind;
	int n, nnzz;
	double *a, *x0;
	
	MatrixData matrixData;

	CgSparse();
	CgSparse(double* a, double* x0);
	CgSparse(double** A, int** Ind, int n, int nnzz);
	CgSparse(MatrixData matrixData, double* a, double* x0);
	~CgSparse();

	void arrayCpy(double *x, double *y, const int n);
	void intArrayCpy(int *x, int *y, const int n);

	// დამხმარე ფუნქცია მასივების სწრაფი გამრავლებისთვის
	inline double vecProd(double *x, double *y, const int n);

	// დამხმარე ფუნქცია მეჩხერი მატრიცის სტრიქონის ვექტორზე სწრაფი გამრავლებისთვის
	inline double matrVecProd(double *v, int * ind, double *y);

	// დამხმარე ფუნქცია სიმეტრიული მატრიცის და მასივის სწრაფი გამრავლებისთვის
	void MatrixByVector(double **m, int **index, double *x, double* res);

	virtual double fillMatrix(); 
	virtual double minimal();
	virtual double minimalCholesky();
	virtual void generateNewYs();
	virtual void beforeMinimal();
	int binarySearch(int * arr, int l, int r, int x);
	std::pair<double**, int**>* cholesky();
	double * forwardSubstitution(std::pair<double**, int**>* L);
	double * backwardSubstitution(std::pair<double**, int**>* U, double * b, bool isTransposed);
};

#endif // CG_SPARSE_H