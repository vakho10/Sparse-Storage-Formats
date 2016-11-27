#include "stdafx.h"

#include <random>
#include <iostream>
#include <algorithm>

#include "CgSparse.h"

CgSparse::CgSparse(MatrixData matrixData, double* a, double* x0)
{
	this->matrixData = matrixData;
	this->a = a;
	this->x0 = x0;
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
	Ind[n] = (int*)malloc(c  * sizeof(int));
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