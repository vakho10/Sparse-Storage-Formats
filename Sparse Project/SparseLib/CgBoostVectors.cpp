#include "stdafx.h"

#include "CgBoostVectors.h"

#include <random>

CgBoostVectors::CgBoostVectors(MatrixData matrixData, double* a, double* x0)
{
	this->matrixData = matrixData;

	int n = matrixData.n;
	this->x0 = ublas::vector<long double>(n);
	for (int i = 0; i < n; i++) {
		this->x0(i) = x0[i];
	}
	this->a = ublas::vector<long double>(n);
	for (int i = 0; i < n; i++) {
		this->a(i) = a[i];
	}
	this->A = ublas::compressed_matrix<double, ublas::row_major, 0>(n, n);
}

CgBoostVectors::~CgBoostVectors() {
}

double CgBoostVectors::fillMatrix()
{
	TimeCounter::StartCounter();
	int n, nnzz;

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

	double v; 	int i, j;
	// ბოლო სტრიქონის გარდა
	for (int ii = 0; ii < nnzz - 1; ++ii)
	{
		// წავიკითხოთ ასეთი რიგით: j,i,v
		// j:
		while (data[pos] != ' ') 
			++pos;
		data[pos] = '\0';
		j = (int)atoi(data);
		++pos;
		data += pos;
		pos = 0;
		// i:
		while (data[pos] != ' ') 
			++pos;
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
		
		// Filling up vectors
		A(i - 1, j - 1) = v;
		A(j - 1, i - 1) = v;

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

		// Filling up vectors
		A(i - 1, j - 1) = v;
		A(j - 1, i - 1) = v;
	}

	free(buffer);
	return TimeCounter::GetCounter();
}

double CgBoostVectors::minimal()
{
	TimeCounter::StartCounter();

	int dimension = matrixData.n;
	ublas::vector<long double> r(dimension);

	r = a - prod(A, x0);

	ublas::vector<long double> d(r);

	long double deltaNew = inner_prod(r, r);

	ublas::vector<long double> q(dimension);

	while (deltaNew >= EPS*EPS) {
		q = prod(A, d);
		long double alpha(deltaNew / inner_prod(d, q));
		x0 += (alpha*d);
		r -= (alpha * q);
		long double deltaOld(deltaNew);
		deltaNew = inner_prod(r, r);
		long double beta(deltaNew / deltaOld);
		d = r + (beta*d);
	}

	return TimeCounter::GetCounter();
}


void CgBoostVectors::generateNewYs()
{
	std::random_device rd;
	std::default_random_engine dre(rd());
	std::uniform_real_distribution<double> di(matrixData.min, matrixData.max);
	for (int i = 0; i < matrixData.n; i++) {
		a(i) = di(dre);
	}
}

void CgBoostVectors::beforeMinimal() {
	std::fill(x0.begin(), x0.end(), 0);
}