#include "stdafx.h"

#include "TimeCounter.h"

#include <iostream>
#include <random>

#include "CgVectors.h"

CgVectors::CgVectors(MatrixData matrixData, double* a, double* x0)
{
	this->matrixData = matrixData;

	int n	= matrixData.n;
	this->row1		= std::vector<long double>(n);
	this->row2		= std::vector<long double>(n);
	this->AA		= std::vector<std::vector<long double>>(n, row1);
	this->AA0		= std::vector<std::vector<long double>>(n, row2);

	this->aa = std::vector<long double>(n);
	for (int i = 0; i < n; i++) {
		this->aa[i] = a[i];
	}

	this->x00 = std::vector<long double>(n);
	for (int i = 0; i < n; i++) {
		this->x00[i] = x0[i];
	}
}

CgVectors::~CgVectors() {
}

double CgVectors::fillMatrix()
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
		AA[i - 1][j - 1] = v;
		AA[j - 1][i - 1] = v;

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
		AA[i - 1][j - 1] = v;
		AA[j - 1][i - 1] = v;
	}

	free(buffer);
	return TimeCounter::GetCounter();
}

double CgVectors::minimal()
{
	// Calculation time
	TimeCounter::StartCounter();

	// Get dimension
	int dimension = matrixData.n;

	std::vector<long double> r(dimension);
	for (int i = 0; i < dimension; i++)
		r[i] = aa[i] - inner_product(AA[i].begin(), AA[i].end(), x00.begin(), 0.);

	std::vector<long double> d(r);

	double deltaNew(inner_product(r.begin(), r.end(), r.begin(), 0.));

	std::vector<long double> q(dimension);

	while (deltaNew >= EPS*EPS) {
		for (int i = 0; i < dimension; i++)
			q[i] = inner_product(AA[i].begin(), AA[i].end(), d.begin(), 0.);

		double alpha(deltaNew / inner_product(d.begin(), d.end(), q.begin(), 0.));
		for (int i = 0; i < dimension; i++)
			x00[i] += alpha*d[i];

		for (int i = 0; i < dimension; i++)
			r[i] -= alpha*q[i];

		double deltaOld(deltaNew);
		deltaNew = inner_product(r.begin(), r.end(), r.begin(), 0.);

		double beta(deltaNew / deltaOld);

		for (int i = 0; i < dimension; i++)
			d[i] = r[i] + (beta*d[i]);

		//cout << deltaNew << "\t\t" << EPSS*EPSS << "\t\t" << countt << endl;
	}

	return TimeCounter::GetCounter();
}

void CgVectors::generateNewYs() 
{
	std::random_device rd; 
	std::default_random_engine dre(rd());
	std::uniform_real_distribution<double> di(matrixData.min, matrixData.max);
	for (int i = 0; i < matrixData.n; i++) {
		this->aa[i] = di(dre);
	}
}

void CgVectors::beforeMinimal() {
	std::fill(x00.begin(), x00.end(), 0);
}