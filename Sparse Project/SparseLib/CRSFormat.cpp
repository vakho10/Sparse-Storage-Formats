#include "stdafx.h"

#include "CRSFormat.h"

#include <random>
#include <iostream>
#include <algorithm>

using namespace std;

CRSFormat::CRSFormat(MatrixData matrixData, double * a, double * x0)
{
	this->matrixData = matrixData;
	this->a = a;
	this->x0 = x0;
}

CRSFormat::~CRSFormat()
{
}

double CRSFormat::fillMatrix()
{
	TimeCounter::StartCounter();

	/* open an existing file for reading */
	FILE *infile = fopen(matrixData.path.c_str(), "r");

	/* declare a file pointer */
	char * buffer;
	long numbytes;

	/* if the file does not exist */
	if (infile == NULL)
		cout << "the file does not exist!" << endl;

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
		cout << "memory error!" << endl;

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

	//რაოდენობები წავიკითხეთ
	while (data[pos] != ' ') ++pos;
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

	//დავიწყოთ მეჩხერი არანულოვანი ქვემატრიცის შევსება
	//მეხსიერების გამოყოფა
	AA = (double *)malloc(sizeof(double)*nnzz);
	JA = (int*)malloc(sizeof(int)*nnzz);
	IA = (int*)malloc(sizeof(int)*(n + 1));
	IA[0] = 0;


	int row = 1;
	int iaInd = 1;
	int c = 0;		//მთვლელი სტრიქონების რაოდენობებისთვის

	double v; 	int i, j;
	//ბოლო სტრიქონის გარდა
	for (int ii = 0; ii < nnzz; ++ii)
	{
		//წავიკითხოთ ასეთი რიგით: j,i,v
		//j:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		j = (int)atoi(data);
		JA[ii] = j - 1;
		++pos;
		data += pos;
		pos = 0;
		//i:
		while (data[pos] != ' ') ++pos;
		data[pos] = '\0';
		i = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		while (data[pos] != '\n') ++pos;
		data[pos] = '\0';
		v = (double)atof(data);
		AA[ii] = v;
		++pos;
		data += pos;
		pos = 0;
		if (row != i)
		{
			IA[iaInd] = IA[iaInd - 1] + c;
			++iaInd;
			row = i;
			c = 0;
		}
		++c;
	}
	IA[n] = IA[n - 1] + c;

	free(buffer);

	//მარჯვენა მხარის შევსება, ფაზური ცვლადის ინიციალიზება
	a = (double *)malloc(sizeof(double)*n);
	x0 = (double *)malloc(sizeof(double)*n);
	default_random_engine dre;
	uniform_real_distribution<double> di(0, 20);
	for (int i = 0; i < n; i++)
	{
		a[i] = di(dre);			//ფაილში რომ წერია?
		x0[i] = 0.;
	}

	return TimeCounter::GetCounter();
}

void CRSFormat::arrayCpy(double *x, double *y, const int n)
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
void CRSFormat::intArrayCpy(int *x, int *y, const int n)
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

//დამხმარე ფუნქცია მასივების სწრაფი გამრავლებისთვის
inline double CRSFormat::vecProd(double *x, double *y, const int n)
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

//დამხმარე ფუნქცია სიმეტრიული მატრიცის და მასივის სწრაფი გამრავლებისთვის
void CRSFormat::MatrixByVectorCRS(double *aa, int *ja, int *ia, int k, double *x, double* res)
{
	int i, j;
	int left, right;

	for (i = 0; i < k; ++i)
		res[i] = 0.;

	for (i = 0; i < k; ++i)
	{
		left = ia[i];
		right = ia[i + 1];

		res[i] += aa[left] * x[ja[left]];

		for (j = left + 1; j < right; ++j)
		{
			res[i] += aa[j] * x[ja[j]];
			res[ja[j]] += aa[j] * x[i];
		}
	}
	/**/
}

double CRSFormat::minimal()
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
		MatrixByVectorCRS(AA, JA, IA, n, d, q);

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

void CRSFormat::generateNewYs()
{
	std::random_device rd;
	std::default_random_engine dre(rd());
	std::uniform_real_distribution<double> di(matrixData.min, matrixData.max);
	for (int i = 0; i < matrixData.n; i++)
	{
		this->a[i] = di(dre);
	}
}

void CRSFormat::beforeMinimal()
{
	memset(x0, 0, sizeof(double)* matrixData.n);
}
