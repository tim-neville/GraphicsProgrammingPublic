#include "tims.h"

const float deg_to_rad = 0.01745329251;
clock_t timer_start;

void DebugAddress(void* pointer)
{
	printf("Address: %p\n", pointer);
}

void DebugText(const char* text)
{
	printf("%s\n", text);
}

void DebugInt(const char* text, int value)
{
	printf("%s%d\n", text, value);
}

void DebugFloat(char* text, float value)
{
	printf("%s%f\n", text, value);
}

void DebugError(const char* text)
{
	fprintf(stderr, "%s\n", text);
}

void StartTimer()
{
	timer_start = clock();
}

void StopTimer(const char* text)
{
	double time_spent = ((double)(clock() - timer_start) /  CLOCKS_PER_SEC) * 1000.0;
	printf("%s (ms): %.2f\n", text, time_spent);
}

float Remap_Float(float value, float low1, float high1, float low2, float high2)
{
	return low2 + (value - low1) * (high2 - low2) / (high1 - low1);
}

//uint32_t depth_to_greyscale(float depth)
//{
//	uint32_t grey = (uint32_t)(depth * UINT32_MAX);
//	grey = (grey << 8) | (grey >> 24);
//	return grey;
//}
