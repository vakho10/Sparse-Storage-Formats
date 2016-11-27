#ifndef MATRIX_DATA_H
#define MATRIX_DATA_H

#ifdef SPARSELIB_EXPORTS
#define SPARSELIB_API __declspec(dllexport) 
#else
#define SPARSELIB_API __declspec(dllimport) 
#endif

#include <string>
#include <fstream>

class SPARSELIB_API MatrixData
{
public:
	double min;
	double max;
	std::string name;
	std::string path;
	std::string yPath; // Path where Ys are stored
	int n;
	int nnz;

	/**
	 *	Compares matrices using their sizes.	
	 */
	static bool compareByN(MatrixData& a, MatrixData& b) 
	{
		return a.n < b.n;
	}

	/**
	 *	Compares matrices using number of their nonzero elements.
	 */
	static bool compareByNNZ(MatrixData& a, MatrixData& b) 
	{
		return a.nnz < b.nnz;
	}

	/**
	 *	Reads the file and ignores comments (which start with '%' symbol).
	 */
	static std::ifstream ignoreComments(std::string path) 
	{
		std::ifstream ifs(path);
		while (!ifs.eof()) {
			int c = ifs.peek(); // peek character
			std::string tmpStr;
			if (c != '%')
				break;
			getline(ifs, tmpStr);
		}
		return ifs;
	}
};

#endif // MATRIX_DATA_H