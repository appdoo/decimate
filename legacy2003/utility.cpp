// ============================================================================
// utils.cpp  
// 
// utility functions
//
// ============================================================================
#include "utility.h"

__int64 __start;
__int64 __stop;
__int64 __Frequency;

// retrieve the clock frequency
void setup_timer() 
{
   // setup frequency fot time measurements
  	QueryPerformanceFrequency((LARGE_INTEGER *)&__Frequency );
}

// read start time
void start_timer() 
{
   QueryPerformanceCounter((LARGE_INTEGER *) &__start);
}

// return time
double stop_timer() 
{
   QueryPerformanceCounter((LARGE_INTEGER *) &__stop);
   return ((double)(__stop - __start) / (double )__Frequency) * 1000.0;
}

char *progress = "|/-\\|/-\\";

// shows a progress information
void show_progress()
{
   // progress index
   static int index = 0;
   index = (index + 1) % 8;

   printf("\r%c", progress[index]);
}
