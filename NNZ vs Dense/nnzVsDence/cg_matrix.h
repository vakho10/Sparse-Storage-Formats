void fillData(string s)
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
	while (s[pos] != '_') ++ pos;
	s[pos] = '\0';
	n = (int)atoi(s.c_str());
	//memory allocation for matrix A
	A = (double**)malloc(n*sizeof(double*));
	for (int i = 0; i < n; ++i)
		A[i] = (double*)malloc(n*sizeof(double));
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
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			while (data[pos] != '\t') ++pos;
			data[pos] = '\0';
			A[i][j] = (double)atof(data);
			++pos;
			data += pos;
			pos = 0;
		}
		++data;
	}
	free(buffer);
}

void getMinimal(string s)
{
	fillData(s);
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
		for (int i = 0; i < n; i++)
			q[i] = vecProd(A[i], d, n);
	
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
