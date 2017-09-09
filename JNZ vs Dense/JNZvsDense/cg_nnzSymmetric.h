void fillNnzSymmetricData(string s)
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
	n = (int)atoi(s.c_str());
	//memory allocation for matries A and Ind
	A = (double**)malloc(n*sizeof(double*));
	Ind = (int**)malloc((n + 1)*sizeof(int*));
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
		for (int j = 0; j < i; j++)
		{
			while (data[pos] != '\t') ++pos;
			data[pos] = '\0';
			++pos;
			data += pos;
			pos = 0;
		}
		for (int j = i; j < n; j++)
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
		Ind[i + 1] = (int*)malloc(counter*sizeof(int));
		intArrayCpy(ind, Ind[i + 1], counter);
		counter = 0;
	}
	free(buffer);
/*		
	for (int i = 0; i < n; i++)
	{	
	for (int j = 0; j < Ind[0][i+1]; j++)
		cout << A[i][j] << '\t';
	cout << endl;
	}
	for (int i = -1; i < n; i++)
	{
		for (int j = 0; j < Ind[0][i + 1]; j++)
			cout << Ind[i+1][j] << '\t';
		cout << endl;
	}
*/
	delete[] ind;
	delete[] val;
}

void MatrixByVector2(double **m, int **index, double *x, double* res)
{
	const int k(Ind[0][0]);

	int i, j;
	double *p;
	int *q;
	int size;
	int t;

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
			t = q[j];
			res[i] += p[j] * x[t];
			res[t] += p[j] * x[i];
		}
	}
}
void getMinimal2(string s)
{
	fillNnzSymmetricData(s);
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
		MatrixByVector2(A, Ind, d, q);

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
	for (int i = 0; i < n; ++i)
		delete[] A[i];
	delete[] A;
}
