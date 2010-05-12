#
# LX_DETECT_BLUEGENE
#
# Defines some compile-time tests to detect BlueGene architectures.
# This will AC_DEFINE the following macros if detected:
#
#    BLUEGENE_L      If BlueGene/L architecture is detected
#    BLUEGENE_P      If BlueGene/P architecture is detected
#
AC_DEFUN([LX_DETECT_BLUEGENE],
[
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#ifndef __blrts__
#error "not BlueGene/L"
#endif /* __blrts__ */
]])],
  [
   AC_DEFINE([BLUEGENE_L],[1],[Define if we're compiling for a BlueGene/L system])
  ])

  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[#ifdef __linux__
#include <linux/config.h>
#endif /* __linux__ */
#if !(defined(__bgp__) || defined(CONFIG_BGP))
#error "not BlueGene/P"
#endif /* !(defined(__bgp__) || defined(CONFIG_BGP)) */
]])],
  [
   AC_DEFINE([BLUEGENE_P],[1],[Define if we're compiling for a BlueGene/P system])
  ])
])
