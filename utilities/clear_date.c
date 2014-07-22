#include "visibility.h"

void
clear_date (/* timetag) */
struct date *timetag)
    {
    timetag->year = 0;
    timetag->day = 0;
    timetag->hour = 0;
    timetag->minute = 0;
    timetag->second = 0.0;
    }
