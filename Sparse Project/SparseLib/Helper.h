#ifndef HELPERS_H
#define HELPERS_H

#ifdef SPARSELIB_EXPORTS
#define SPARSELIB_API __declspec(dllexport) 
#else
#define SPARSELIB_API __declspec(dllimport) 
#endif

#include <random>
#include <string>

#include "MatrixData.h"

namespace helper {

	SPARSELIB_API void setMatrixDataExtras(MatrixData& matrixData, bool& skipMatrix);
}

#endif // HELPERS_H