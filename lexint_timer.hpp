#pragma once

#include <map>
#include <set>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <sys/time.h>

namespace LeXInt
{
    /// This timer class measures the elapsed time between two events. Timers can be
    /// started and stopped repeatedly. The total time as well as the average time
    /// between two events can be queried using the total() and average() methods,
    /// respectively.
    struct timer {
        timespec t_start;
        bool running;
        double elapsed = 0.0;
        unsigned counter;

        timer() {
            counter = 0;
            running = false;
        }

        void start() {
            clock_gettime(CLOCK_REALTIME, &t_start);
            running = true;
        }

        void restart() {
            elapsed = 0.0;
            counter = 0;
        }

        double stop() {
            if(running == false) {
                ::std::cout << "WARNING: timer::stop() has been called without calling timer::start() first." << ::std::endl;
                return 0.0;
            } else {
                timespec t_end;
                clock_gettime(CLOCK_REALTIME, &t_end);
                int sec  = t_end.tv_sec-t_start.tv_sec;
                double nsec = ((double)(t_end.tv_nsec-t_start.tv_nsec));
                if(nsec < 0.0) {
                    nsec += 1e9;
                    sec--;
                }
                double t = (double)sec + nsec/1e9;
                counter++;
                elapsed += t;
                return t;
            }
        }

        double total() {
            return elapsed;
        }

        double average() {
            return elapsed/double(counter);
        }

        unsigned count() {
            return counter;
        }
    };
}
