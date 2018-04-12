#include "stdafx.h"

#include <random>
#include <iostream>
#include <algorithm>
#include <utility>

#include "CgSparse.h"

CgSparse::CgSparse(MatrixData matrixData, double* a, double* x0) : CgSparse(a, x0)
{
	this->matrixData = matrixData;
}

CgSparse::CgSparse(double * a, double * x0)
{
	this->a = a;
	this->x0 = x0;
}

CgSparse::CgSparse(double ** A, int ** Ind, int n, int nnzz)
{
	this->A = A;
	this->Ind = Ind;
	this->n = n;
	this->nnzz = nnzz;
}

CgSparse::CgSparse()
{
	// Empty constructor
}

CgSparse::~CgSparse() {
}

void CgSparse::arrayCpy(double *x, double *y, const int n)
{
	int i, n5;
	if (n <= 0) return;
	n5 = n % 5;
	for (i = 0; i < n5; i++)  y[i] = x[i];
	for (; i < n; i += 5)
	{
		y[i] = x[i]; y[i + 1] = x[i + 1]; y[i + 2] = x[i + 2];  y[i + 3] = x[i + 3]; y[i + 4] = x[i + 4];
	}
}

void CgSparse::intArrayCpy(int *x, int *y, const int n)
{
	int i, n5;
	if (n <= 0) return;
	n5 = n % 5;
	for (i = 0; i < n5; i++)  y[i] = x[i];
	for (; i < n; i += 5)
	{
		y[i] = x[i]; y[i + 1] = x[i + 1]; y[i + 2] = x[i + 2];  y[i + 3] = x[i + 3]; y[i + 4] = x[i + 4];
	}
}

// დამხმარე ფუნქცია მასივების სწრაფი გამრავლებისთვის
inline double CgSparse::vecProd(double *x, double *y, const int n)
{
	int i, n5;
	double sum(0.0);
	if (n <= 0) return sum;
	n5 = n % 5;
	for (i = 0; i < n5; i++) sum += x[i] * y[i];
	for (; i < n; i += 5)
	{
		sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
			+ x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
	}
	return sum;
}

// დამხმარე ფუნქცია მეჩხერი მატრიცის სტრიქონის ვექტორზე სწრაფი გამრავლებისთვის
inline double CgSparse::matrVecProd(double *v, int * ind, double *y)
{
	double sum(0.0);
	int size = ind[0];
	for (int k = 1; k < size; k++)
		sum += v[k] * y[ind[k]];
	return sum;
}

// დამხმარე ფუნქცია სიმეტრიული მატრიცის და მასივის სწრაფი გამრავლებისთვის
void CgSparse::MatrixByVector(double **m, int **index, double *x, double* res)
{
	const int k(Ind[0][0]);

	int i, j;
	double *p;
	int *q;
	int size;

	for (i = 0; i < k; ++i)
		res[i] = 0.;

	for (i = 0; i < k; ++i)
	{
		p = m[i];
		q = index[i + 1];
		size = index[0][i + 1];
		res[i] += p[0] * x[q[0]];

		for (j = 1; j < size; ++j)
		{
			res[i] += p[j] * x[q[j]];
			res[q[j]] += p[j] * x[i];
		}
	}
}

double CgSparse::fillMatrix()
{
	TimeCounter::StartCounter();

	/* open an existing file for reading */
	FILE *infile = fopen(matrixData.path.c_str(), "r");

	/* declare a file pointer */
	char * buffer;
	long numbytes;

	/* if the file does not exist */
	if (infile == NULL)
		std::cout << "the file does not exist!" << std::endl;

	/* Get the number of bytes */
	fseek(infile, 0L, SEEK_END);
	numbytes = ftell(infile);

	/* reset the file position indicator to
	the beginning of the file */
	fseek(infile, 0L, SEEK_SET);

	/* grab sufficient memory for the
	buffer to hold the text */
	buffer = (char*)calloc(numbytes, sizeof(char));

	/* memory error */
	if (buffer == NULL)
		std::cout << "memory error!" << std::endl;

	/* copy all the text into the buffer */
	fread(buffer, sizeof(char), numbytes, infile);
	fclose(infile);

	// Ignore comment section
	size_t pos = 0;
	char *data = buffer;
	while (data[pos] == '%')
	{
		++pos;
		while (data[pos] != '\n')
			++pos;
		data += (pos + 1);
		pos = 0;
	}

	// რაოდენობები წავიკითხეთ	
	while (data[pos] != ' ')
		++pos;
	data[pos] = '\0';
	n = (int)atoi(data);
	++pos;
	data += pos;

	// There is second n in the matrix file
	pos = 0;
	while (data[pos] != ' ')
		++pos;
	data[pos] = '\0';
	n = (int)atoi(data);
	++pos;
	data += pos;

	pos = 0;
	while (data[pos] != '\n') ++pos;
	data[pos] = '\0';
	nnzz = (int)atoi(data);
	++pos;
	data += pos;
	pos = 0;
	// დავიწყოთ მეჩხერი აანულოვანი ქვემატრიცის შევსება
	// მეხსიერების გამოყოფა
	A = (double **)malloc(sizeof(double *)*n);
	Ind = (int**)malloc((n + 1) * sizeof(int*));
	Ind[0] = (int*)malloc((n + 1) * sizeof(int));
	Ind[0][0] = n;


	// მაქს. სიგრძის ორი ვექტორი, სტრიქონში ინდექსის და მნიშვნელობებისთვის
	int* ind = (int *)malloc(sizeof(int)*n);
	double*	val = (double *)malloc(sizeof(double)*n);

	int row = 1;
	int c = 0;		//მთვლელი
	// ფაილიდან სტრიქონის წამოსაღები 
	double v; 	int i, j;
	// ბოლო სტრიქონის გარდა
	for (int ii = 0; ii < nnzz - 1; ++ii)
	{
		// წავიკითხოთ ასეთი რიგით: j,i,v
		// j:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		j = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		// i:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		i = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		while (data[pos] != '\n') ++pos;
		data[pos] = '\0';
		v = (double)atof(data);
		++pos;
		data += pos;
		pos = 0;
		if (row != i)
		{
			--row;
			A[row] = (double*)malloc(c * sizeof(double));
			Ind[row + 1] = (int*)malloc(c * sizeof(int));

			arrayCpy(val, A[row], c);
			intArrayCpy(ind, Ind[row + 1], c);
			Ind[0][row + 1] = c;
			row = i;
			////////cout << c << endl;////////
			c = 0;
		}
		val[c] = v;
		ind[c] = j - 1;
		++c;
		/**/
	}

	{// ბოლო სტრიქონი
		// წავიკითხოთ ასეთი რიგით: j,i,v
		// j:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		j = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		// i:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		i = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		v = (double)atof(data);
		// მოვრჩით კითხვას..

		if (row != i)
		{
			--row;
			A[row] = (double*)malloc(c * sizeof(double));
			Ind[row + 1] = (int*)malloc(c * sizeof(int));
			arrayCpy(val, A[row], c);
			intArrayCpy(ind, Ind[row + 1], c);
			Ind[0][row + 1] = c;
			row = i;
			c = 0;
		}
		val[c] = v;
		ind[c] = j - 1;
		++c;
	}
	A[n - 1] = (double*)malloc(c * sizeof(double));
	Ind[n] = (int*)malloc(c * sizeof(int));
	arrayCpy(val, A[n - 1], c);
	intArrayCpy(ind, Ind[n], c);
	Ind[0][n] = c;

	delete[] ind;
	delete[] val;
	free(buffer);

	return TimeCounter::GetCounter();
}

double CgSparse::minimal()
{
	TimeCounter::StartCounter();

	double *r, *d, *q;
	r = (double *)malloc(sizeof(double)*n);
	d = (double *)malloc(sizeof(double)*n);
	q = (double *)malloc(sizeof(double)*n);

	for (int i = 0; i < n; i++)
		r[i] = d[i] = q[i] = 0.;

	for (int i = 0; i < n; i++)
		d[i] = r[i] = a[i];

	double deltaNew(vecProd(r, r, n));

	int count = 0;
	while (deltaNew >= EPSS*EPSS)
	{

		MatrixByVector(A, Ind, d, q);

		double alpha = deltaNew / vecProd(d, q, n);
		for (int i = 0; i < n; i++)
			x0[i] += alpha*d[i];

		for (int i = 0; i < n; i++)
			r[i] -= alpha*q[i];

		double deltaOld(deltaNew);
		deltaNew = vecProd(r, r, n);
		double beta(deltaNew / deltaOld);
		for (int i = 0; i < n; i++)
			d[i] = r[i] + (beta*d[i]);
		count++;
	}

	// Free memory
	free(r); r = NULL;
	free(d); d = NULL;
	free(q); q = NULL;

	return TimeCounter::GetCounter();
}

double CgSparse::minimalCholesky()
{
	TimeCounter::StartCounter();

	/*std::pair<double**, int**> L = cholesky();
	double* y = forwardSubstitution(L);
	double* x0 = backwardSubstitution(L, y, false);*/

	return TimeCounter::GetCounter();
}

void CgSparse::generateNewYs()
{
	std::random_device rd;
	std::default_random_engine dre(rd());
	std::uniform_real_distribution<double> di(matrixData.min, matrixData.max);
	for (int i = 0; i < matrixData.n; i++) {
		this->a[i] = di(dre);
	}
}

void CgSparse::beforeMinimal() {
	memset(x0, 0, sizeof(double)* matrixData.n);
}

// Iterative search
int CgSparse::iterativeSearch(int* arr, int l, int r, int x)
{
	for (int i = l; i < r; i++) {
		if (arr[i] == x) {
			return i;
		}
	}
	return -1;
}

// A recursive binary search function. It returns 
// location of x in given array arr[l..r] is present, 
// otherwise -1
int CgSparse::binarySearch(int* arr, int l, int r, int x)
{
	if (r >= l)
	{
		int mid = l + (r - l) / 2;

		// If the element is present at the middle 
		// itself
		if (arr[mid] == x)
			return mid;

		// If element is smaller than mid, then 
		// it can only be present in left subarray
		if (arr[mid] > x)
			return binarySearch(arr, l, mid - 1, x);

		// Else the element can only be present
		// in right subarray
		return binarySearch(arr, mid + 1, r, x);
	}

	// We reach here when element is not 
	// present in array
	return -1;
}


std::pair<double**, int**>* CgSparse::cholesky()
{
	// Construct L matrix's value and index arrays
	int **L_Ind = (int**)malloc((n + 1) * sizeof(int*));
	double** L_A = (double **)malloc(sizeof(double *) * n);
	L_Ind[0] = (int*)malloc((n + 1) * sizeof(int));
	L_Ind[0][0] = n;
	//for (int i = 1; i <= n; ++i) {
	//	L_Ind[i] = (int*)malloc(i * sizeof(int));
	//	L_Ind[0][i] = i;
	//	L_A[i - 1] = (double*)calloc(i, sizeof(double));
	//}

	for (int i = 0; i < n; i++)
	{
		// create temporary row 
		double* tmpRow = (double*)calloc(n, sizeof(double));

		int nnzsInTmpRow = 0;
		for (int j = 0; j <= i; j++)
		{
			double sum = 0;

			if (j == i) // summation for diagnols	
			{
				for (int k = 0; k < j; k++)
					sum += pow(tmpRow[k], 2); // TODO may be sped up
					//sum += pow(L_A[j][k], 2);


				// Find column index and element in A matrix
				int colIndex = iterativeSearch(Ind[j + 1], 0, Ind[0][j + 1], j);
				//L_A[j][j] = sqrt(((colIndex == -1) ? 0 : A[j][colIndex]) - sum);
				tmpRow[j] = sqrt(((colIndex == -1) ? 0 : A[j][colIndex]) - sum);
			}
			else {

				// Evaluating L(i, j) using L(j, j)
				for (int k = 0; k < j; k++) {
					int colIndex = binarySearch(L_Ind[j + 1], 0, L_Ind[0][j + 1], k); // [j, k]
					sum += (tmpRow[k] * (colIndex == -1 ? 0 : L_A[j][colIndex]));
					//sum += (L_A[i][k] * L_A[j][k]); // TODO same here as for diagonals
				}

				// Find column index and element in A matrix
				// We only save upper triangular matrix in sparse memory (so i is always less that j)!
				int findI = i <= j ? i : j;
				int findJ = i <= j ? j : i;
				int colIndexA = iterativeSearch(Ind[findI + 1], 0, Ind[0][findI + 1], findJ);
				//L_A[i][j] = ((colIndex == -1 ? 0 : A[findI][colIndex]) - sum) / L_A[j][j];
				int colIndexL_A = binarySearch(L_Ind[j + 1], 0, L_Ind[0][j + 1], j);
				tmpRow[j] = ((colIndexA == -1 ? 0 : A[findI][colIndexA]) - sum) / L_A[j][colIndexL_A]; // denominator must be present!
			}

			if (tmpRow[j] != 0) {
				nnzsInTmpRow++;
			}
		}

		// Store tmp row in real matrix 
		L_Ind[i + 1] = (int*)malloc(nnzsInTmpRow * sizeof(int));
		L_Ind[0][i + 1] = nnzsInTmpRow;
		L_A[i] = (double*)calloc(nnzsInTmpRow, sizeof(double));
		for (int j = 0, ind = 0; j < n; j++) {
			if (tmpRow[j] != 0) {
				L_A[i][ind] = tmpRow[j];
				L_Ind[i + 1][ind] = j;
				++ind;
			}
		}

		delete tmpRow;
	}
	std::pair<double**, int**>* L = new std::pair<double**, int**>();
	L->first = L_A;
	L->second = L_Ind;

	return L;

	//double** lower = new double*[n];


	/*for (int i = 0; i < n; i++) {
		for (int j = 0; j < (i + 1); j++) {
			double s = 0;
			for (int k = 0; k < j; k++) {
				s += L[i * n + k] * L[j * n + k];
			}

			L[i * n + j] = (i == j) ? sqrt(A[i * n + i] - s) : (1.0 / L[j * n + j] * (A[i * n + j] - s));
		}
	}*/

	// Decomposing a matrix into Lower Triangular
	//for (int i = 0; i < n; i++) 
	//{
	//	for (int j = 0; j <= i; j++) 
	//	{
	//		double sum = 0;

	//		if (j == i) // summation for diagnols	
	//		{
	//			for (int k = 0; k < j; k++)
	//				sum += pow(lower[j * n + k], 2);
	//			lower[j * n + j] = sqrt(matrix[j * n + j] -
	//				sum);
	//		}
	//		else {

	//			// Evaluating L(i, j) using L(j, j)
	//			for (int k = 0; k < j; k++)
	//				sum += (lower[i][k] * lower[j][k]);
	//			lower[i][j] = (matrix[i][j] - sum) /
	//				lower[j][j];
	//		}
	//	}
	//}
}

double* CgSparse::forwardSubstitution(std::pair<double**, int**> *L) {

	double *y = (double*)calloc(n, sizeof(double));
	if (y == NULL)
		exit(EXIT_FAILURE);

	for (int i = 0; i < n; i++) {
		y[i] = a[i];

		for (int j = 0, colIndex; j <= i - 1; j++) {
			colIndex = binarySearch(L->second[i + 1], 0, L->second[0][i + 1], j);
			y[i] -= (colIndex == -1 ? 0 : L->first[i][colIndex]) * y[j];
			//y[i] -= L->first[i][j] * y[j];
		}

		// Diagonal entries are real and positive!
		y[i] /= L->first[i][
			binarySearch(L->second[i + 1], 0, L->second[0][i + 1], i)
		];
	}
	return y;
}

double* CgSparse::backwardSubstitution(std::pair<double**, int**> *U, double *b, bool isTransposed = true) {

	double *x = (double*)calloc(n, sizeof(double));
	if (x == NULL)
		exit(EXIT_FAILURE);

	// Same codes but with i and j swapped for matrix U
	if (isTransposed) {
		for (int i = n - 1; i >= 0; i--) {
			x[i] = b[i];

			for (int j = i + 1, colIndex; j < n; j++) {
				colIndex = binarySearch(U->second[i + 1], 0, U->second[0][i + 1], j);
				x[i] -= (colIndex == -1 ? 0 : U->first[i][colIndex]) * x[j];
				//x[i] -= U->first[i][j] * x[j];
			}
			
			// Diagonal entries are real and positive!
			x[i] /= U->first[i][
				binarySearch(U->second[i + 1], 0, U->second[0][i + 1], i)
			];
			//x[i] /= U->first[i][i];
		}
	}
	else {
		for (int i = n - 1; i >= 0; i--) {
			x[i] = b[i];

			/*x[i] = b[i];

			for (int j = i + 1; j < n; j++)
				x[i] -= U[i + j * n] * x[j];

			x[i] /= U[i + i * n];*/

			for (int j = i + 1, colIndex; j < n; j++) {
				colIndex = binarySearch(U->second[j + 1], 0, U->second[0][j + 1], i);
				x[i] -= (colIndex == -1 ? 0 : U->first[j][colIndex]) * x[j];
				//x[i] -= U->first[i][j] * x[j];
			}

			// Diagonal entries are real and positive!
			x[i] /= U->first[i][
				binarySearch(U->second[i + 1], 0, U->second[0][i + 1], i)
			];
			//x[i] /= U->first[i][i];
		}
	}
	return x;
}