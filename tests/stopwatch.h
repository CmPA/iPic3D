#include <sys/time.h>
#include <assert.h>
#include <stdint.h>

#define myuint64_t int
//#define myuint64_t uint64_t

static inline myuint64_t tv_to_sec(struct timeval tv){
   return tv.tv_sec + tv.tv_usec/1000000;
}

static inline myuint64_t tv_to_ms(struct timeval tv){
   return tv.tv_sec*1000 + tv.tv_usec/1000;
}

static inline myuint64_t tv_to_us(struct timeval tv){
   return tv.tv_sec*1000000 + tv.tv_usec;
}


static inline struct timeval add_tv(const struct timeval a, const struct timeval b){
   const struct timeval res = { a.tv_sec+b.tv_sec, a.tv_usec+b.tv_usec };
   return res;
}

static inline struct timeval diff_tv(const struct timeval start, const struct timeval stop){
   const struct timeval diff = {stop.tv_sec - start.tv_sec, stop.tv_usec - start.tv_usec};
   return diff;
}

typedef enum {START, STOP, LAP, RESET} sw_action_t;
static inline int valid_sw_action(sw_action_t t){
   return t == START || t == STOP || t == LAP || RESET;
}


typedef enum {OFF, STARTED, STOPPED} sw_state_t;
static inline int valid_sw_state(sw_state_t s){
   return s == OFF || s == STOPPED || s == STARTED;
}


typedef struct {
   sw_state_t state;
   struct timeval total, now, last;
} stopwatch_t;

static inline struct timeval sw_start(stopwatch_t * const sw ){
   sw->state = STARTED;
   gettimeofday( &sw->now, 0 );
   sw->last = sw->now;
   return sw->total;
}

static inline struct timeval sw_stop(stopwatch_t * const sw ){
   sw->state = STOPPED;
   gettimeofday( &sw->now, 0 );
   sw->total = add_tv(sw->total, diff_tv(sw->last, sw->now));
   sw->last = sw->now;
   return sw->total;
}

static inline struct timeval sw_lap(stopwatch_t * const sw ){
   gettimeofday( &sw->now, 0 );
   const struct timeval elapsed = diff_tv(sw->last, sw->now);
   sw->total = add_tv(sw->total, elapsed);
   sw->last = sw->now;
   return elapsed;
}

static inline struct timeval sw_reset(stopwatch_t * const sw ){
   const static stopwatch_t sw_off = { OFF, {0, 0}, {0, 0}, {0, 0} };
   sw->state = OFF;
   *sw = sw_off;
   return sw->total;
}

static inline struct timeval stopwatch_mt(stopwatch_t * const sw, sw_action_t action){
   typedef struct timeval stopwatch_func_t( stopwatch_t *sw );

   static stopwatch_func_t * const  stopwatch_transitions[3][4] = { {sw_start, 0, 0, 0},
                                                                        {0, sw_stop, sw_lap, 0},
                                                                        {sw_start, 0, 0, sw_reset}};
   assert(sw != 0);
   assert(valid_sw_action(action));
   assert(valid_sw_state(sw->state));
   assert(stopwatch_transitions[sw->state][action] != 0);
   return stopwatch_transitions[sw->state][action](sw);
}


struct timeval stopwatch(sw_action_t action){
   static stopwatch_t sw = { OFF, {0, 0}, {0, 0}, {0, 0} };
   return stopwatch_mt(&sw, action);
}


