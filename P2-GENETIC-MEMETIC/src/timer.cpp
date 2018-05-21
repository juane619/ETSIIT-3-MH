
#include "timer.h"

clock_t start_time;

double elapsed;

void start_timers() {
    start_time = clock();
}

double elapsed_time() {
    elapsed = clock() - start_time;
    return elapsed / CLOCKS_PER_SEC;
}



