#include "stdafx.h"

#include "Helper.h"
#include <iostream>

void helper::setMatrixDataExtras(MatrixData& matrixData, bool& skipMatrix) {
	skipMatrix = false;

	double minElem = 10101;
	double maxElem = -10101;

	/* open an existing file for reading */
	FILE *infile = fopen(matrixData.path.c_str(), "r");

	/* declare a file pointer */
	char* buffer;
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

	int n1, n2;

	// რაოდენობები წავიკითხეთ	
	while (data[pos] != ' ')
		++pos;
	data[pos] = '\0';
	n1 = (int)atoi(data);
	++pos;
	data += pos;

	// There is second n in the matrix file
	pos = 0;
	while (data[pos] != ' ')
		++pos;
	data[pos] = '\0';
	n2 = (int)atoi(data);
	++pos;
	data += pos;

	if (n1 != n2) {
		skipMatrix = true;
		free(buffer);
		return;
	}
	matrixData.n = n1;

	pos = 0;
	while (data[pos] != '\n') ++pos;
	data[pos] = '\0';
	matrixData.nnz = (int)atoi(data);
	++pos;
	data += pos;
	pos = 0;

	double v; 	int i, j;
	// ბოლო სტრიქონის გარდა
	for (int ii = 0; ii < matrixData.nnz - 1; ++ii)
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

		if (v > maxElem) {
			maxElem = v;
		}
		if (v < minElem) {
			minElem = v;
		}
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

		if (v > maxElem) {
			maxElem = v;
		}
		if (v < minElem) {
			minElem = v;
		}
	}

	free(buffer);

	matrixData.min = minElem;
	matrixData.max = maxElem;
}