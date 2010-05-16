#include "timing.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif //HAVE_CONFIG_H

#include <cmath>
using namespace std;


#if defined(BLUEGENE_L)
// -------------------------------------------------------- //
// Timing code for BlueGene/L
// -------------------------------------------------------- //
#include <rts.h>

// this will return number of nanoseconds in a single BGL cycle
// use for converting from cycle units returned by rts_gettimebase
// to nanoseconds.
static double get_ns_per_cycle() {
  BGLPersonality personality;
  if (rts_get_personality(&personality, sizeof(personality)) != 0)
    return 0;
  return 1.0e9/((double) personality.clockHz);
}

// returns time in nanoseconds.
timing_t get_time_ns () {
  static double ns_per_cycle = get_ns_per_cycle();
  return (timing_t)(ns_per_cycle * rts_get_timebase());
}


#elif defined(BLUEGENE_P)
// -------------------------------------------------------- //
// Timing code for BlueGene/P
// -------------------------------------------------------- //
#define SPRN_TBRL         0x10C        // Time Base Read Lower Register (user & sup R/O)
#define SPRN_TBRU         0x10D        // Time Base Read Upper Register (user & sup R/O)
#define BGP_NS_PER_CYCLE  (1.0/0.85)   // Nanoseconds per cycle on BGP (850Mhz clock)

#define _bgp_mfspr(SPRN) ({ \
   unsigned int tmp; \
   do { \
      asm volatile ("mfspr %0,%1" : "=&r" (tmp) : "i" (SPRN) : "memory" ); \
   } while(0); \
   tmp; \
})

union bgp_time_reg {
  unsigned int ul[2];
  unsigned long long ull;
};

static inline unsigned long long timebase() {
  bgp_time_reg reg;
  unsigned int utmp;
  
  do {
    utmp      = _bgp_mfspr(SPRN_TBRU);
    reg.ul[1] = _bgp_mfspr(SPRN_TBRL);
    reg.ul[0] = _bgp_mfspr(SPRN_TBRU);
  }
  while( utmp != reg.ul[0] );
  
  return reg.ull;
}

timing_t get_time_ns() {
  return llround(BGP_NS_PER_CYCLE * timebase());
}



#elif (defined(HAVE_CLOCK_GETTIME) || defined(HAVE_LIBRT))
// -------------------------------------------------------- //
// Timing code using Linux hires timers.
// -------------------------------------------------------- //

#include <ctime>
#include <sys/time.h>

timing_t get_time_ns() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (ts.tv_sec * 1000000000ll + ts.tv_nsec);
}


#elif defined(HAVE_GETTIMEOFDAY)
// -------------------------------------------------------- //
// Generic timing code using gettimeofday.
// -------------------------------------------------------- //

#include <ctime>
#include <sys/time.h>

timing_t get_time_ns() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000000000ll + tv.tv_usec * 1000ll;
}


#else // if we get to here, we don't even have gettimeofday.
#error "NO SUPPORTED TIMING FUNCTIONS FOUND!"
#endif // types of timers
