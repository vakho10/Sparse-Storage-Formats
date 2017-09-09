#include <Windows.h>

void arrayCpy(double *x, double *y, const int n)
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
void intArrayCpy(int *x, int *y, const int n)
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
void fillNnzData(string s)
{
	/* open an existing file for reading */
	FILE *infile = fopen(s.c_str(), "r");
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

	//გავიგოთ ელემენტების რაოდენობა
	size_t pos = 0;
	while (s[pos] != '_') ++pos;
	s[pos] = '\0';
	s = s.substr(s.rfind("\\") + 1);
	n = (int)atoi(s.c_str());
	//memory allocation for matries A and Ind
	A = (double**)malloc(n*sizeof(double*));
	Ind = (int**)malloc((n+1)*sizeof(int*));
	Ind[0] = (int*)malloc((n + 1)*sizeof(int));
	Ind[0][0] = n;
	//მარჯვენა მხარის შევსება, ფაზური ცვლადის ინიციალიზება
	a = (double *)malloc(sizeof(double)*n);
	x0 = (double *)malloc(sizeof(double)*n);
	default_random_engine dre1;
	uniform_real_distribution<double> di(0, 20);
	for (int i = 0; i < n; i++)
	{
		a[i] = di(dre1);
		x0[i] = 0.;
	}

	//მატრიცის  შევსება
	pos = 0;
	char *data = buffer;
	//მაქს. სიგრძის ორი ვექტორი, სტრიქონში ინდექსის და მნიშვნელობებისთვის
	int* ind = (int *)malloc(sizeof(int)*n);
	double*	val = (double *)malloc(sizeof(double)*n);
	int counter{ 0 };
	//შევსება
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			while (data[pos] != '\t') ++pos;
			data[pos] = '\0';
			val[counter] = (double)atof(data);
			if (val[counter])
			{
				ind[counter] = j;
				++counter;
			}
			++pos;
			data += pos;
			pos = 0;
		}
		++data;
		A[i] = (double*)malloc(counter*sizeof(double));
		arrayCpy(val, A[i], counter);
		Ind[0][i + 1] = counter;
		Ind[i+1] = (int*)malloc(counter*sizeof(int));
		intArrayCpy(ind, Ind[i + 1], counter);
		counter = 0;
	}
	free(buffer);
	delete[] ind;
	delete[] val;
}

void MatrixByVector(double **m, int **index, double *x, double* res)
{

	const int k(index[0][0]);
	int i, j, n5;
	double *p;
	int *q;
	int size;
	double result;

	for (i = 0; i < k; ++i)
	{
		result = 0.;
		p = m[i];
		q = index[i + 1];
		size = index[0][i + 1];
		n5 = size % 5;
		for (j = 0; j < n5; ++j)
			result += p[j] * x[q[j]];
		for (; j < size; j += 5)
		{
			result += p[j] * x[q[j]] + p[j + 1] * x[q[j + 1]]
				+ p[j + 2] * x[q[j + 2]] + p[j + 3] * x[q[j + 3]]
				+ p[j + 4] * x[q[j + 4]];
		}
		res[i] = result;
	}
}
long long getMinimal1(string s)
{
	fillNnzData(s);

	auto st_local = chrono::high_resolution_clock::now();

	double *r, *d, *q;
	r = (double *)malloc(sizeof(double)*n);
	d = (double *)malloc(sizeof(double)*n);
	q = (double *)malloc(sizeof(double)*n);

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
	for (int i = 0; i < n; ++i) {
		//free(A[i]);
		delete[] A[i];
	}
	//free(A);
	delete[] A;

	auto difference = chrono::high_resolution_clock::now() - st_local;
	auto local_time = chrono::duration_cast<chrono::milliseconds>(difference);
	return local_time.count();
}
