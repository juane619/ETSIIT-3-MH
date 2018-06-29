#ifndef TIMER_H
#define TIMER_H

#include <time.h>

extern clock_t start_time;

extern double elapsed;

int time_expired();

void start_timers();

double elapsed_time();

typedef enum type_timer {REAL, VIRTUAL} TIMER_TYPE;

#endif