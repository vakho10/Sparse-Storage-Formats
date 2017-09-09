#include <Windows.h>
#include <random>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>

#include "..\SparseLib\Helper.h"
#include "..\SparseLib\Evaluator.h"
#include "..\SparseLib\MatrixData.h"
#include "..\SparseLib\TimeCounter.h"

#include "..\SparseLib\CgSparse.h"
#include "..\SparseLib\CgVectors.h"
#include "..\SparseLib\CgBoostVectors.h"
#include "..\SparseLib\CgMaps.h"
#include "..\SparseLib\CRSFormat.h"

// test comment
//#include "boost\numeric\ublas\matrix_vector.hpp"
//#include "boost\numeric\ublas\vector_sparse.hpp"
//#include "boost\numeric\ublas\matrix_sparse.hpp"
//#include "boost\numeric\ublas\operation.hpp"

#include "boost\filesystem.hpp"
#include "boost\regex.hpp"
#include "boost\algorithm\string\predicate.hpp" // For string comparison

using namespace boost::filesystem;
using namespace std;

namespace ublas = boost::numeric::ublas;

// Initialize TimeCounter's static variables (once per compilation) [TimeCounter.h]
//double TimeCounter::PCFreq = 0.0;
//__int64 TimeCounter::CounterStart = 0;

pair<double, double> calculateAverageWorkingTime(Evaluator* abs, int typeIndex, int genType, double* a, double* x0, int n)
{
	pair<double, double> workingTimes;
	int t = 1;

	switch (typeIndex) {
	case 1:
		t = 10; // For small
		break;
	case 2:
		t = 5; // For medium
		break;
	default:
		t = 1; // For large
		break;
	}

	if (t > 1)
	{
		std::vector<double> fillVec;
		std::vector<double> solveVec;
		for (int g = t; g >= 1; g--)
		{
			if (genType == 2 && g < t) {
				// After next iterations generate new Ys
				abs->generateNewYs();
			}

			workingTimes.first = abs->fillMatrix();
			workingTimes.second = abs->getMinimal();

			fillVec.push_back(workingTimes.first);
			solveVec.push_back(workingTimes.second);
		}

		// If it worked with fixed Ys then calculate median. 
		// Otherwise, calculate average.
		if (genType == 1) {
			sort(fillVec.begin(), fillVec.end());
			sort(solveVec.begin(), solveVec.end());
			if (t % 2 == 0) {
				int ind1 = floor(t / 2) - 1;
				int ind2 = floor(t / 2);
				workingTimes.first = (fillVec[ind1] + fillVec[ind2]) / 2;
				workingTimes.second = (solveVec[ind1] + solveVec[ind2]) / 2;
			}
			else {
				int index = floor(t / 2);
				workingTimes.first = fillVec[index];
				workingTimes.second = solveVec[index];
			}
		}
		else {
			// Get average fill time
			double avgFill = 0.0;
			for (int i = 0; i < fillVec.size(); i++) {
				avgFill += fillVec[i];
			}
			workingTimes.first = (avgFill / fillVec.size());

			// Get average solve time
			double avgSolve = 0.0;
			for (int i = 0; i < solveVec.size(); i++) {
				avgSolve += solveVec[i];
			}
			workingTimes.second = (avgSolve / solveVec.size());
		}
	}
	else {
		// Worked only once (needs no median nor average)
		workingTimes.first = abs->fillMatrix();
		workingTimes.second = abs->getMinimal();
	}
	return workingTimes;
}

int main()
{
	// Random number generator for Ys
	random_device rd; // Random at each execution
	default_random_engine dre(rd());

	// Choose the size of matrices
	int typeIndex = -1;
	string types[3] = { "small", "medium", "large" };
	cout << "Choose matrix type:\n";
	cout << "1 - small\n2 - medium\n3 - large\n";
	while (true) {
		cin >> typeIndex;
		if (typeIndex > 3 || typeIndex < 1) {
			cout << "Wrong number! Try again." << endl;
			continue;
		}
		else {
			break;
		}
	}

	// Choose the generation type of Ys
	int genType = -1;
	cout << "1 - Existing Ys\n2 - Random Ys\n";
	while (true) {
		cin >> genType;
		if (genType > 2 || genType < 1) {
			cout << "Wrong number! Try again." << endl;
			continue;
		}
		else {
			break;
		}
	}

	// JSON RESULT
	std::ofstream jsonOfs("data.json");
	jsonOfs << "[";

	// Get the files with extension .mtx (specify full or relative path)
	string patternPath = "matrices\\";
	path current_dir(patternPath + types[typeIndex - 1]);
	boost::regex pattern("(.*\\.mtx)");

	/*
		Collect basic information about matrix files (path, name, n, nnz)
		*/
	std::vector<MatrixData> matrixDatas = std::vector<MatrixData>();

	for (recursive_directory_iterator iter(current_dir), end; iter != end; ++iter)
	{
		string name = iter->path().filename().string();

		if (regex_match(name, pattern))
		{
			// Or if it ends with 'b.mtx'!
			if (boost::algorithm::ends_with(name, "b.mtx") || boost::algorithm::ends_with(name, "coord.mtx")) {
				continue;
			}

			// Ignore this matrix
			if (boost::algorithm::starts_with(name, "bibd")) {
				continue;
			}

			string path = iter->path().string();		

			// Construct matrix data object
			MatrixData tmpData;
			tmpData.name = name;
			tmpData.path = path;
			
			bool skipMatrix = false;
			helper::setMatrixDataExtras(tmpData, skipMatrix);

			// n1 != n2
			if (skipMatrix) {
				continue;
			}

			uniform_real_distribution<double> di(tmpData.min, tmpData.max);

			// Create '*_Ys.mtx' file for matrix (if it doesn't exist).
			string partPath = path.substr(0, path.find(".mtx", 0));
			string yPath = partPath + "_Ys.txt";
			tmpData.yPath = yPath; // Save path in MatrixData object
			bool forceOverwrite = false; // Be careful with this parameter!
			if (!boost::filesystem::exists(yPath) || forceOverwrite) {
				ofstream ofs(yPath);
				for (int i = 0; i < tmpData.n; i++) {
					ofs << di(dre) << '\n';
				}
				ofs.close();
			}

			matrixDatas.push_back(tmpData);
		}
	}

	// Sort by N
	sort(matrixDatas.begin(), matrixDatas.end(), MatrixData::compareByN);

	cout << matrixDatas.size() << " matrices..." << endl;

	try {
		// Loop through matrices and do the calculations
		for (unsigned int i = 25; i < matrixDatas.size(); i++)
		{
			string name = matrixDatas[i].name;
			string path = matrixDatas[i].path;
			int n = matrixDatas[i].n;
			int nnz = matrixDatas[i].nnz;

			double *a, *x0;

			// მარჯვენა მხარის შევსება, ფაზური ცვლადის ინიციალიზება
			a = (double *)malloc(sizeof(double)*n);
			x0 = (double *)calloc(n, sizeof(double)); // allocates memory and fills with zeroes.

			if (genType == 1) {
				// FIXME Change to fast reading!!!!!
				std::ifstream ifs(matrixDatas[i].yPath);
				for (int i = 0; i < n; i++)
				{
					// Globally available static variables
					ifs >> a[i];
				}
				ifs.close();
			}
			else {
				uniform_real_distribution<double> di(matrixDatas[i].min, matrixDatas[i].max);
				for (unsigned int h = 0; h < n; h++) {
					a[h] = di(dre);
				}
			}

			cout << (i + 1) << endl;
			Evaluator* abs = nullptr;
			pair<double, double> workingTimes;

			jsonOfs << "{";
			jsonOfs << "\"name\":\"" << name << "\",";
			jsonOfs << "\"n\":" << n << ",";
			jsonOfs << "\"nnz\":" << nnz << ",";

			jsonOfs << "\"results\":{";

			// CG_Boost_Vectors
			cout << path << ", N: " << n << ", NNZ: " << nnz << endl;

			abs = new CgSparse(matrixDatas[i], a, x0);
			// Version 1			
			//TimeCounter::StartCounter(); // Fill
			//abs->fillMatrix();
			//workingTimes.first = TimeCounter::GetCounter();			
			//TimeCounter::StartCounter(); // Solve
			//abs->getMinimal();
			//workingTimes.second = TimeCounter::GetCounter();

			// Version 2
			/*workingTimes.first = abs->fillMatrix();
			workingTimes.second = abs->getMinimal();*/
			workingTimes = calculateAverageWorkingTime(abs, typeIndex, genType, a, x0, n);
			cout << "[CG_SPARSE] FILL: " << workingTimes.first << ", CALCULATE: " << workingTimes.second << endl;
			jsonOfs << "\"cg_sparse\":{\"fill\":" << workingTimes.first << ",\"solve\":" << workingTimes.second << "},";
			delete abs;

			// Clear Xs
			//memset(x0, 0, sizeof(double) * n);

			// CRS Format
			abs = new CRSFormat(matrixDatas[i], a, x0);			
			workingTimes = calculateAverageWorkingTime(abs, typeIndex, genType, a, x0, n);
			cout << "[CG_CRS_FORMAT] FILL: " << workingTimes.first << ", CALCULATE: " << workingTimes.second << endl;
			jsonOfs << "\"cg_crs_format\":{\"fill\":" << workingTimes.first << ",\"solve\":" << workingTimes.second << "},";
			delete abs;

			abs = new CgBoostVectors(matrixDatas[i], a, x0);
			// Version 1			
			//TimeCounter::StartCounter(); // Fill
			//abs->fillMatrix();
			//workingTimes.first = TimeCounter::GetCounter();
			//TimeCounter::StartCounter(); // Solve
			//abs->getMinimal();
			//workingTimes.second = TimeCounter::GetCounter();

			// Version 2
			/*workingTimes.first = abs->fillMatrix();
			workingTimes.second = abs->getMinimal();*/
			workingTimes = calculateAverageWorkingTime(abs, typeIndex, genType, a, x0, n);
			cout << "[CG_BOOST_VECTORS] FILL: " << workingTimes.first << ", CALCULATE: " << workingTimes.second << endl;
			jsonOfs << "\"cg_boost_vectors\":{\"fill\":" << workingTimes.first << ",\"solve\":" << workingTimes.second << "},";
			delete abs;

			/*abs = new CgVectors(matrixDatas[i], a, x0);
			workingTimes = calculateAverageWorkingTime(abs, typeIndex, x0, n);
			cout << "[CG_VECTORS] FILL: " << workingTimes.first << ", CALCULATE: " << workingTimes.second << endl;
			jsonOfs << "\"cg_vectors\":{\"fill\":" << workingTimes.first << ",\"solve\":" << workingTimes.second << "}";
			delete abs;*/

			/*abs = new CG_Arrays(matrixDatas[i]);
			workingTimes = abs->getMinimal();
			cout << "[CG_ARRAYS] FILL: " << workingTimes.first << ", CALCULATE: " << workingTimes.second << endl;
			jsonOfs << "\"cg_arrays\":{\"fill\":" << workingTimes.first << ",\"solve\":" << workingTimes.second << "},";
			delete abs;*/

			// Clear Xs
			//memset(x0, 0, sizeof(double) * n);

			abs = new CgMaps(matrixDatas[i], a, x0);
			// Version 1			
			//TimeCounter::StartCounter(); // Fill
			//abs->fillMatrix();
			//workingTimes.first = TimeCounter::GetCounter();
			//TimeCounter::StartCounter(); // Solve
			//abs->getMinimal();
			//workingTimes.second = TimeCounter::GetCounter();

			// Version 2
			/*workingTimes.first = abs->fillMatrix();
			workingTimes.second = abs->getMinimal();*/
			workingTimes = calculateAverageWorkingTime(abs, typeIndex, genType, a, x0, n);
			cout << "[CG_MAPS] FILL: " << workingTimes.first << ", CALCULATE: " << workingTimes.second << endl;
			jsonOfs << "\"cg_maps\":{\"fill\":" << workingTimes.first << ",\"solve\":" << workingTimes.second << "}";
			delete abs;

			cout << endl;

			jsonOfs << "}}";
			if (i != matrixDatas.size() - 1) {
				jsonOfs << ",";
			}

			jsonOfs.flush();

			// Free global variables for next iteration
			free(a);
			free(x0);
		}
		jsonOfs << "]";
		jsonOfs.close();

		/*TimeCounter::StartCounter();
		Sleep(1000);
		cout << "Testing one second! => " << GetCounter() << "\n";

		clock_t tStart = clock();
		TimeCounter::StartCounter();
		vector<double> vec(10000);
		for (unsigned int i = 0; i < vec.size(); i++) {
		vec[i] = fRand(1, 100);
		}
		sort(vec.begin(), vec.end());
		printf("Time taken: %.2fms\n --- %f", (double)(clock() - tStart) / CLOCKS_PER_SEC * 1000, GetCounter());*/

		// Solution used from http://stackoverflow.com/questions/1739259/how-to-use-queryperformancecounter
	}
	catch (std::overflow_error e) {
		cout << "ERROR!" << endl;
	}

	//using namespace boost::numeric::ublas;
	//mapped_matrix<double> m(3, 3, 3 * 3); // n, m, nnz
	//for (unsigned i = 0; i < m.size1(); ++i)
	//	for (unsigned j = 0; j < m.size2(); ++j)
	//		m(i, j) = 3 * i + j;


	//for (unsigned i = 0; i < m.size1(); ++i) {
	//	for (unsigned j = 0; j < m.size2(); ++j) {
	//		cout << m(i, j) << "  ";
	//	}
	//	cout << endl;
	//}

	return 0;
}