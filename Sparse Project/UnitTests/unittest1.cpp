#include "stdafx.h"
#include "CppUnitTest.h"

#include "..\SparseLib\TimeCounter.h"
#include "..\SparseLib\CgBoostVectors.h"

#include <string>
#include <algorithm>

// These defines are used to get path to the project directory
#define STRINGIFY(x) #x
#define EXPAND(x) STRINGIFY(x)

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTests
{
	TEST_CLASS(UnitTest1)
	{
	public:

		// Tests if TimeCounter gives the right time.
		TEST_METHOD(TestTimeCounter)
		{
			double expectedTime = 250;

			TimeCounter::StartCounter();
			Sleep(expectedTime);
			double time = TimeCounter::GetCounter();

			// Should be in between (num - 1) and (num + 1)
			Assert::IsTrue((expectedTime - 1) < time < (expectedTime + 1));

			// Double checking...
			expectedTime = 250;

			TimeCounter::StartCounter();
			Sleep(expectedTime);
			time = TimeCounter::GetCounter();

			// Should be in between...
			Assert::IsTrue((expectedTime - 1) < time < (expectedTime + 1));
		}

		// Test if solve works as expected.
		TEST_METHOD(TestSolve)
		{
			// Ax = b
			double a[3][3] = { { 25.5, 0.25, 0.1 }, { 0.75, 54.8, 0 }, { 0, 4.5, 7.8 } };
			double x[] = { 8.6, 1.0, 11.5 };
			double b[3];

			for (int i = 0; i < 3; i++) {
				b[i] = 0;
			}

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					b[i] += (a[i][j] * x[j]);
				}
			}

			double expB[3] = { 220.7, 61.25, 94.2 };

			// Check if calculation is correct...
			for (int i = 0; i < 3; i++) {
				Assert::AreEqual(expB[i], b[i]);
			}

			MatrixData matrixData = MatrixData();
			matrixData.n = 3;
			matrixData.min = 0;
			matrixData.max = 7;
			matrixData.nnz = 9;

			double x0[3] = { 0, 0, 0 };
			CgBoostVectors boostVecs = CgBoostVectors(matrixData, b, x0);

			// Fill by hand
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					boostVecs.A(i, j) = a[i][j];
				}
			}

			boostVecs.getMinimal();

			for (int i = 0; i < 3; i++) {
				x0[i] = boostVecs.x0(i);
			}

			// Check if the answer is in epsilon range
			for (int i = 0; i < 3; i++) {
				Assert::IsTrue(x[i] + boostVecs.EPS >= x0[i] && x[i] - boostVecs.EPS <= x0[i]);
			}
		}

		// Test whether fill works fine or not.
		TEST_METHOD(TestFill)
		{
			// Ax = b
			double a[3][3] = { { 14.0, 0, 8.0 }, { 0, 0, 7.0 }, { 8.0, 7.0, 0 } };
			double x[] = { 8.0, 1.0, 11.0 };
			double b[3];

			for (int i = 0; i < 3; i++) {
				b[i] = 0;
			}

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					b[i] += (a[i][j] * x[j]);
				}
			}

			double expB[3] = { 200, 77, 71 };

			// Check if calculation is correct...
			for (int i = 0; i < 3; i++) {
				Assert::AreEqual(expB[i], b[i]);
			}

			MatrixData matrixData = MatrixData();
			matrixData.n = 3;
			matrixData.min = 0;
			matrixData.max = 7;
			matrixData.nnz = 9;

			std::string s = EXPAND(UNITTESTPRJ);
			s.erase(0, 1); // erase the first quote
			s.erase(s.size() - 2); // erase the last quote and the dot
			matrixData.path = s + "dummy_matrix.mtx";

			double x0[3] = { 0, 0, 0 };
			CgBoostVectors boostVecs = CgBoostVectors(matrixData, b, x0);

			boostVecs.fillMatrix();

			typedef boost::numeric::ublas::compressed_matrix<double, ublas::row_major, 0>::iterator1 it1_t;
			typedef boost::numeric::ublas::compressed_matrix<double, ublas::row_major, 0>::iterator2 it2_t;

			std::vector<double> vec;
			for (it1_t it1 = boostVecs.A.begin1(); it1 != boostVecs.A.end1(); it1++)
			{
				for (it2_t it2 = it1.begin(); it2 != it1.end(); it2++)
				{
					vec.push_back(*it2);
				}
			}

			// And because it is symetric...
			double nnzs[5] = { 14, 7, 7, 8, 8 };

			// Sort them and check (FIXME bad test?!)
			std:sort(vec.begin(), vec.end());
			std::sort(nnzs, nnzs + 5);

			for (int i = 0; i < 3; i++) {
				Assert::AreEqual(nnzs[i], vec[i]);
			}
		}

	};
}