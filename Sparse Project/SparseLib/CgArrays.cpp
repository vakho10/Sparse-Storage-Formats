#include "stdafx.h"

#include "CgArrays.h"

#include <random>

CgArrays::CgArrays(MatrixData matrixData)
{
	this->matrixData = matrixData;
}

CgArrays::~CgArrays() {
}

/**
 *	Helper function.
 */
inline double CgArrays::vecProd(double *x, double *y, const int n)
{
	int i, n5;
	double sum(0.0);
	if (n <= 0) return sum;
	n5 = n % 5;
	for (i = 0; i < n5; i++) {
		sum += x[i] * y[i];
	}
	for (; i < n; i += 5) {
		sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2] + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
	}
	return sum;
}

double CgArrays::fillMatrix()
{
	TimeCounter::StartCounter();
	std::ifstream ifs = MatrixData::ignoreComments(matrixData.path);
	int n, nnz;
	ifs >> n >> n >> nnz;

	// Get dimension
	int nARR = matrixData.n;

	// Memory allocation
	A_ARR = (double **)malloc(sizeof(double *)*nARR);
	for (int i = 0; i < nARR; i++) {

		A_ARR[i] = (double *)malloc(sizeof(double)*nARR);
	}
	a_ARR = (double *)malloc(sizeof(double)*nARR);
	x0_ARR = (double *)malloc(sizeof(double)*nARR);

	// Initialization
	for (int i = 0; i < nARR; ++i) {
		a_ARR[i] = 0;
		for (int j = 0; j < nARR; ++j) {
			A_ARR[i][j] = 0.;
		}
	}

	// Filling...
	int row, column;
	double value;
	for (int i = 0; i < nnz; i++) {
		ifs >> row >> column >> value;
		A_ARR[row - 1][column - 1] = value;
		A_ARR[column - 1][row - 1] = value;
	}

	ifs.close();

	// Input Ys
	std::ifstream inp("ys.txt");
	for (int i = 0; i < nARR; i++) {
		inp >> a_ARR[i];
		x0_ARR[i] = 0.;
	}
	inp.close();

	return TimeCounter::GetCounter();
}

double CgArrays::minimal()
{
	// Calculation time
	TimeCounter::StartCounter();

	// Get dimension
	int nARR = matrixData.n;

	double *r, *d, *q;
	r = (double *)malloc(sizeof(double)*nARR);
	d = (double *)malloc(sizeof(double)*nARR);
	q = (double *)malloc(sizeof(double)*nARR);

	for (int i = 0; i < nARR; i++)
		r[i] = d[i] = q[i] = 0.;

	for (int i = 0; i < nARR; i++)
		d[i] = r[i] = a_ARR[i] - vecProd(A_ARR[i], x0_ARR, nARR);

	double deltaNew(vecProd(r, r, nARR));

	while (deltaNew >= EPS*EPS) {
		for (int i = 0; i < nARR; i++)
			q[i] = vecProd(A_ARR[i], d, nARR);

		double alpha = deltaNew / vecProd(d, q, nARR);
		for (int i = 0; i < nARR; i++)
			x0_ARR[i] += alpha*d[i];

		for (int i = 0; i < nARR; i++)
			r[i] -= alpha*q[i];

		double deltaOld(deltaNew);
		deltaNew = vecProd(r, r, nARR);

		double beta(deltaNew / deltaOld);

		for (int i = 0; i < nARR; i++)
			d[i] = r[i] + (beta*d[i]);

		//cout << deltaNew << "\t\t" << EPS*EPS << endl;
	}

	// Cleaning
	for (int i = 0; i < nARR; i++)
		delete A_ARR[i];

	delete A_ARR;
	delete a_ARR;
	delete x0_ARR;

	return TimeCounter::GetCounter();
}

void CgArrays::generateNewYs()
{
	std::random_device rd;
	std::default_random_engine dre(rd());
	std::uniform_real_distribution<double> di(matrixData.min, matrixData.max);
	for (int i = 0; i < matrixData.n; i++) {
		//a_MAP(i) = di(dre);
	}
}

void CgArrays::beforeMinimal() {
	memset(x0_ARR, 0, sizeof(double)* matrixData.n);
}