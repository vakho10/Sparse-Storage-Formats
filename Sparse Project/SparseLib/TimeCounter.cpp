#include "stdafx.h"

#include "TimeCounter.h"

double TimeCounter::PCFreq = 0.0;
__int64 TimeCounter::CounterStart = 0;

void TimeCounter::StartCounter() {
	LARGE_INTEGER li;
	if (!QueryPerformanceFrequency(&li))
		std::printf("QueryPerformanceFrequency failed!\n");
	PCFreq = double(li.QuadPart) / 1000.0; // Specifies milliseconds
	QueryPerformanceCounter(&li);
	CounterStart = li.QuadPart;
}

double TimeCounter::GetCounter() {
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return double(li.QuadPart - CounterStart) / PCFreq;
}

//LARGE_INTEGER frequency;        // ticks per second
//LARGE_INTEGER t1, t2;           // ticks
//double elapsedTime;
//
//void TimeCounter::StartCounter() {
//	// get ticks per second
//	QueryPerformanceFrequency(&frequency);
//}
//
//double TimeCounter::GetCounter() {
//	QueryPerformanceCounter(&t2);
//	return double((t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart);
//}