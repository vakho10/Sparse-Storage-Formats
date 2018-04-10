#define USE_EXISTING_YS true
#define RESULT_JSON_PATH "results.json"
#define MIN_VALUE -100
#define MAX_VALUE +100

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>

#include <chrono>
#include <random>
#include <fstream>

#include "boost\filesystem.hpp"
#include "boost\regex.hpp"
#include "boost\algorithm\string\predicate.hpp"

#include "mmio.c"
#include "json.hpp"

#include "..\SparseLib\MatrixData.h"
#include "..\SparseLib\Evaluator.h"
#include "..\SparseLib\CgSparse.h"

int countNonZeroes(double* vec, int n) {
	int sum = 0;
	for (int i = 0; i < n; i++) {
		if (vec[i] != 0) {
			sum += 1;
		}
	}
	return sum;
}

double* cholesky(double *A, int n)
{
	double *L = (double*)calloc(n * n, sizeof(double));
	if (L == NULL) {
		exit(EXIT_FAILURE);
	}

	/*for (int i = 0; i < n; i++) {
		for (int j = 0; j < (i + 1); j++) {
			double s = 0;
			for (int k = 0; k < j; k++) {
				s += L[i * n + k] * L[j * n + k];
			}

			L[i * n + j] = (i == j) ? sqrt(A[i * n + i] - s) : (1.0 / L[j * n + j] * (A[i * n + j] - s));
		}
	}*/

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			double sum = 0;

			if (j == i) // summation for diagnols	
			{
				for (int k = 0; k < j; k++)
					sum += pow(L[j * n + k], 2);
				L[j * n + j] = sqrt(A[j * n + j] - sum);
			}
			else {

				// Evaluating L(i, j) using L(j, j)
				for (int k = 0; k < j; k++)
					sum += (L[i * n + k] * L[j * n + k]);
				L[i * n + j] = (A[i * n + j] - sum) / L[j * n + j];
			}
		}
	}
	return L;
}

void show_matrix(double *A, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			printf("%2.5f ", A[i * n + j]);
		printf("\n");
	}
}

void show_vector(double *v, int n) {
	for (int i = 0; i < n; i++)
		printf("%2.5f ", v[i]);

	printf("\n");
}

double* forwardSubstitution(double *L, double *b, int n) {

	double *y = (double*)calloc(n, sizeof(double));
	if (y == NULL)
		exit(EXIT_FAILURE);

	for (int i = 0; i < n; i++) {
		y[i] = b[i];

		for (int j = 0; j <= i - 1; j++)
			y[i] -= L[j + i * n] * y[j];

		y[i] /= L[i + i * n];
	}
	return y;
}

double* backwardSubstitution(double *U, double *b, int n, bool isTransposed = true) {

	double *x = (double*)calloc(n, sizeof(double));
	if (x == NULL)
		exit(EXIT_FAILURE);

	// Same codes but with i and j swapped for matrix U
	if (isTransposed) {
		for (int i = n - 1; i >= 0; i--) {
			x[i] = b[i];

			for (int j = i + 1; j < n; j++)
				x[i] -= U[j + i * n] * x[j];

			x[i] /= U[i + i * n];
		}
	}
	else {
		for (int i = n - 1; i >= 0; i--) {
			x[i] = b[i];

			for (int j = i + 1; j < n; j++)
				x[i] -= U[i + j * n] * x[j];

			x[i] /= U[i + i * n];
		}
	}
	return x;
}

int main()
{
	using namespace std;
	using namespace boost;
	using namespace boost::filesystem;
	using namespace boost::algorithm;
	using namespace std::chrono;
	using json = nlohmann::json;

	time_point<high_resolution_clock> start, finish;

	// reusable time variable
	long long time;

	// Random number generator for Ys
	random_device rd; // Random at each execution
	default_random_engine dre(rd());

	// Get the files with extension .mtx (specify full or relative path)
	path current_dir("..\\SparseMatrixProject\\matrices\\small");
	regex pattern("(.*\\.mtx)");

	json rootJson;

	// Iterate and read matrices
	for (recursive_directory_iterator iter(current_dir), end; iter != end; ++iter)
	{
		string name = iter->path().filename().string();

		if (regex_match(name, pattern))
		{
			// Or if it ends with 'b.mtx'!
			if (ends_with(name, "b.mtx") || ends_with(name, "coord.mtx")) {
				continue;
			}

			// Ignore this matrix
			if (starts_with(name, "bibd")) {
				continue;
			}

			FILE *f;
			int ret_code;
			MM_typecode matcode;
			int M, N, nz;
			int i, *I, *J;
			double *val, *inlineMat;

			string filePath = iter->path().string();

			// Open file
			if ((f = fopen(filePath.c_str(), "r")) == NULL)
				continue;

			if (mm_read_banner(f, &matcode) != 0)
			{
				printf("Could not process Matrix Market banner.\n");
				continue;
			}

			/*  This is how one can screen matrix types if their application */
			/*  only supports a subset of the Matrix Market data types.      */
			/*if (mm_is_complex(matcode) && mm_is_matrix(matcode) &&
				mm_is_sparse(matcode))
			{
				printf("Sorry, this application does not support ");
				printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
				exit(1);
			}*/

			/* find out size of sparse matrix .... */
			if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
				continue;

			// Ys
			double *b;
			b = (double *)malloc(sizeof(double)*N);
			//x = (double *)calloc(N, sizeof(double)); // allocates memory and fills with zeroes.

			string ysPath = filePath.substr(0, filePath.find(".mtx", 0)) + "_Ys.txt";
			if (USE_EXISTING_YS)
			{
				// Open Ys file
				FILE *ysFile;
				if ((ysFile = fopen(ysPath.c_str(), "r")) == NULL)
					continue;

				for (size_t i = 0; i < N; i++)
					fscanf(ysFile, "%lg\n", &b[i]);

				fclose(ysFile);
			}
			else {
				uniform_real_distribution<double> di(MIN_VALUE, MAX_VALUE);
				for (size_t i = 0; i < N; i++)
					b[i] = di(dre);
			}

			/* reseve memory for matrices */
			I = (int *)malloc(nz * sizeof(int));
			J = (int *)malloc(nz * sizeof(int));
			val = (double *)malloc(nz * sizeof(double));
			inlineMat = (double *)calloc(N * N, sizeof(double));

			/* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
			/*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
			/*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */
			for (i = 0; i < nz; i++)
			{
				fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
				I[i]--;  /* adjust from 1-based to 0-based */
				J[i]--;

				inlineMat[J[i] + (I[i] * N)] = val[i];
			}

			if (f != stdin) fclose(f);

			/************************/
			/* now write out matrix */
			/************************/
			mm_write_banner(stdout, matcode);
			mm_write_mtx_crd_size(stdout, M, N, nz);
			/*for (i = 0; i < nz; i++) {
				fprintf(stdout, "%d %d %20.19g\n", I[i] + 1, J[i] + 1, val[i]);
			}*/

			json matrixJson = { { "n", N }, { "nz", nz } },
				choleskyJson = { { "type", "cholesky" } },
				cgJson = { { "type", "cg" } };

			// [1] Solve with "Cholevsky Decomposition"
			start = high_resolution_clock::now();

			double* L = cholesky(inlineMat, N);
			double* y = forwardSubstitution(L, b, N);
			double* x1 = backwardSubstitution(L, y, N, false);

			finish = high_resolution_clock::now();
			time = duration_cast<nanoseconds>(finish - start).count();
			fprintf(stdout, "Cholevsky: %f milliseconds. \n", time / (double)1000000);
			choleskyJson["solve"] = time / (double)1000000;


			// [2] Solve with "Conjugate Gradient"
			double* x2 = (double *)calloc(N, sizeof(double));
			MatrixData matrixData;
			matrixData.n = N;
			matrixData.nnz = nz;
			matrixData.yPath = ysPath;
			matrixData.path = filePath;
			CgSparse* cgSparse = new CgSparse(matrixData, b, x2);
			cgSparse->fillMatrix();

			printf("New matrixo formatoo \n");
			std::pair<double**, int**>* res = cgSparse->cholesky();
			double* y_sparse = cgSparse->forwardSubstitution(res);
			double* x1_sparse = cgSparse->backwardSubstitution(res, y_sparse, false);

			for (int i = 0; i < N / 20; i++) {
				for (int j = 0; j <= i; j++) {
					fprintf(stdout, "%f ", res->first[i][j]);
				}
				fprintf(stdout, "%d nnz in row and in %d nnz real", countNonZeroes(res->first[i], res->second[0][i + 1]), cgSparse->Ind[0][i + 1]);
				printf("\n");
			}
			show_vector(L, 20);
			fprintf(stdout, "%f \n", L[0]);
			fprintf(stdout, "%f \n", L[N + 0]);
			fprintf(stdout, "%f \n", L[N + 1]);
			fprintf(stdout, "%f \n", L[N * 2 + 0]);
			fprintf(stdout, "%f \n", L[N * 2 + 1]);
			fprintf(stdout, "%f \n", L[N * 2 + 2]);
			fprintf(stdout, "%f \n", L[N * 3 + 0]);
			fprintf(stdout, "%f \n", L[N * 3 + 1]);
			fprintf(stdout, "%f \n", L[N * 3 + 2]);
			fprintf(stdout, "%f \n", L[N * 3 + 3]);

			double solveTime = cgSparse->getMinimal();
			fprintf(stdout, "CG: %f milliseconds.\n", solveTime);
			cgJson["solve"] = solveTime;

			

			matrixJson["results"].push_back(choleskyJson);
			matrixJson["results"].push_back(cgJson);

			rootJson.push_back(matrixJson);

			show_vector(y_sparse, 5);
			show_vector(y, 5);
			show_vector(x1_sparse, 5);
			show_vector(x1, 5);
			show_vector(x2, 5);

			free(I);
			free(J);
			free(val);
			free(inlineMat);

			free(L);
			free(y);
			free(x1);
			free(x2);

			free(b);

			printf("\n");
		}
	}

	// the setw manipulator was overloaded to set the indentation for pretty printing
	cout << setw(4) << rootJson << endl;

	// write prettified JSON to file
	FILE * pFile;
	pFile = fopen(RESULT_JSON_PATH, "w");
	if (pFile != NULL)
	{
		fputs(rootJson.dump().c_str(), pFile);
		fclose(pFile);
	}
}
