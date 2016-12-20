using namespace std;

const double EPSS = 1E-5;
double **A;
int **Ind;
double *a, *x0;
int n;

inline	double vecProd(double *x, double *y, const int n)
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

void make_dataFile
(
	int rows,
	int zerosPercent
)
{
	//memory allocation for matrix A
	double **A;
	A = (double**)malloc(rows*sizeof(double*));
	for(int i = 0; i < rows; ++i)
		A[i] = (double*)malloc(rows*sizeof(double));
	//fill diagonal
	for (int i = 0; i < rows; ++i)
		A[i][i] = 100;
	//filling upper part
	default_random_engine dre;
	uniform_int_distribution<short> uid(0,100);
	for (int i = 0; i < rows; ++i)
	for (int j = i + 1; j < rows; ++j)
		A[i][j] = uid(dre);
	//transform to reals
	uniform_real_distribution<float> urd;
	for (int i = 0; i < rows; ++i)
	for (int j = i ; j < rows; ++j)
		if (A[i][j] > zerosPercent) 
			A[i][j] = urd(dre);
		else A[i][j] = 0.;
	//filling lower part
	for (int i = 1; i < rows; ++i)
	for (int j = 0; j < i; ++j)
		A[i][j] = A[j][i];

	string name{ to_string(rows) };
	name += '_';
	name += to_string(zerosPercent);
	name += ".txt";
	
	freopen(name.c_str(), "w", stdout);
	//print
	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < rows; ++j)
			cout << A[i][j] << '\t';
		cout << endl;
	}
	for (int i = 0; i < rows; ++i)
		delete[] A[i];
	delete[] A;
}	

void make_rightHandSide
(
	int rows,
	int zerosPercent
)
{
	string name{ to_string(rows) };
	name += '_';
	name += to_string(zerosPercent);
	name += "_rightHandSide.txt";

	freopen(name.c_str(), "w", stdout);
	//print
	//fill 
	default_random_engine dre(zerosPercent);
	uniform_real_distribution<float> urd;
	for (int i = 0; i < rows; ++i)
		cout <<  urd(dre);
}
