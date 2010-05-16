#ifndef TIMING_H
#define TIMING_H

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

  /// Size for timings reported by get_time_ns.
  typedef unsigned long long timing_t;


  /// This is defined differently depending on compile-time timing options
  /// and on the platform we're using.  It returns a current time stamp in 
  /// nanoseconds using the most precise timing mechanism available on
  /// the host machine.
  timing_t get_time_ns();
  
#ifdef __cplusplus
}
#endif // __cplusplus

#endif // TIMING_H
