// ============================================================================
// utils.h  
// 
// utility functions
//
// ============================================================================
#ifndef UTILITY_H
#define UTILITY_H

#include <windows.h>
#include <stdio.h>

// retrieve the clock frequency
void setup_timer();
// read start time
void start_timer();
// return time
double stop_timer();
// shows a progress information
void show_progress();

#endif