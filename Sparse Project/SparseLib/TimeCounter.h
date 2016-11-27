#ifndef TIME_COUNTER_H
#define TIME_COUNTER_H

#ifdef SPARSELIB_EXPORTS
#define SPARSELIB_API __declspec(dllexport) 
#else
#define SPARSELIB_API __declspec(dllimport) 
#endif

#include <Windows.h>
#include <string>

/**
 *	This counter is not thread-safe!
 */
class TimeCounter
{
public:
	SPARSELIB_API static double PCFreq;
	SPARSELIB_API static __int64 CounterStart;

	SPARSELIB_API static void StartCounter();
	SPARSELIB_API static double GetCounter();
};

#endif // TIME_COUNTER_H