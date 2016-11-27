#include "stdafx.h"

#include "CgMaps.h"

#include <random>
#include <iostream>

CgMaps::CgMaps(MatrixData matrixData, double* a, double* x0)
{
	this->matrixData = matrixData;
	this->mappedMatrix = mapped_matrix<double>(matrixData.n, matrixData.n, matrixData.nnz); // n, m, nnz

	this->a_MAP = mapped_vector<double>(matrixData.n);
	for (unsigned int i = 0; i < a_MAP.size(); i++) {
		a_MAP(i) = a[i];
	}

	this->x0_MAP = mapped_vector<double>(matrixData.n);
	for (unsigned int i = 0; i < x0_MAP.size(); i++) {
		x0_MAP(i) = x0[i];
	}
}

CgMaps::~CgMaps() {
}

inline double CgMaps::vecProd1(double *x, double *y, const int n)
{
	int i, n5;
	double sum(0.0);
	if (n <= 0) return sum;
	n5 = n % 5;
	for (i = 0; i < n5; i++) sum += x[i] * y[i];
	for (; i < n; i += 5) {
		sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2] + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
	}
	return sum;
}

inline void CgMaps::vecToMap(double *r, double *x0, std::map<int, double> map, const int n)
{
	for (auto it = map.begin(); it != map.end(); it++) {
		int i = it->first / n;

		int j = it->first%n;
		r[i] += it->second*x0[j];
	}
}

double CgMaps::fillMatrix()
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

		// Filling up maps
		int row = i;
		int column = j;

		row--;
		column--;

		mappedMatrix(row, column) = v;
		mappedMatrix(column, row) = v;
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

		// Filling up maps
		int row = i;
		int column = j;

		row--;
		column--;

		mappedMatrix(row, column) = v;
		mappedMatrix(column, row) = v;
	}

	free(buffer);
	return TimeCounter::GetCounter();
}

double CgMaps::minimal()
{
	// Calculation time
	TimeCounter::StartCounter();

	// Get dimension
	int nMAP = matrixData.n;

	/*double *r, *d, *q;
	r = (double *)malloc(sizeof(double)*nMAP);
	d = (double *)malloc(sizeof(double)*nMAP);
	q = (double *)malloc(sizeof(double)*nMAP);*/
	mapped_vector<double> r(nMAP);
	

	/*for (int i = 0; i < nMAP; i++)
		r[i] = d[i] = q[i] = 0.;*/

	r = a_MAP - prod(mappedMatrix, x0_MAP);

	mapped_vector<double> d(r);
	mapped_vector<double> q(nMAP);

	//vecToMap(r, x0_MAP, ARR_MAP, nMAP);
	// Multiply (Vector * Matrix) 

	/*for (int i = 0; i < nMAP; i++)
		d[i] = r[i] = a_MAP[i] - r[i];*/

	double deltaNew = inner_prod(r, r);

	int count = 0;
	while (deltaNew >= EPS*EPS) {
		/*for (int i = 0; i < nMAP; i++)
			q[i] = 0;*/

		q = prod(mappedMatrix, d); // q = M * d

		//vecToMap(q, d, ARR_MAP, nMAP);
		/*for (int i = 0; i < nMAP; i++) {
		cout << q[i] << "\t";
		}*/
		double alpha = deltaNew / inner_prod(d, q);

		for (int i = 0; i < nMAP; i++)
			x0_MAP[i] += alpha*d[i];

		for (int i = 0; i < nMAP; i++)
			r[i] -= alpha*q[i];

		double deltaOld(deltaNew);
		deltaNew = inner_prod(r, r);

		double beta(deltaNew / deltaOld);

		for (int i = 0; i < nMAP; i++)
			d[i] = r[i] + (beta*d[i]);
	}

	return TimeCounter::GetCounter();
}

void CgMaps::generateNewYs()
{
	std::random_device rd;
	std::default_random_engine dre(rd());
	std::uniform_real_distribution<double> di(matrixData.min, matrixData.max);
	for (int i = 0; i < matrixData.n; i++) {
		a_MAP(i) = di(dre);
	}
}

void CgMaps::beforeMinimal() {
	std::fill(x0_MAP.begin(), x0_MAP.end(), 0);
}