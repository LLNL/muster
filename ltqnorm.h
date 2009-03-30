#ifndef STAT_UTILS_H
#define STAT_UTILS_H

/**
 * Based on the code by Peter J. Acklam, and the adaptation by V. Natarajan,
 * available from http://home.online.no/~pjacklam/notes/invnorm/.
 * The code in this file is free to distribute or use for any purpose.
 *
 * Adapted for C++ and inclusion in AMPL by Todd Gamblin, 1/28/2005.  
 * tgamblin@cs.unc.edu
 *
 * Original header comments follow --------------------------------------
 *
 * Z = LTQNORM(P) returns the lower tail quantile for the standard normal
 * distribution function.  I.e., it returns the Z satisfying Pr{X < Z} = P,
 * where X has a standard normal distribution.
 *
 * LTQNORM(P) is the same as SQRT(2) * ERFINV(2*P-1), but the former returns a
 * more accurate value when P is close to zero.
 *
 * The algorithm uses a minimax approximation by rational functions and the
 * result has a relative error less than 1.15e-9.  A last refinement by
 * Halley's rational method is applied to achieve full machine precision.
 *
 * Author:      Peter J. Acklam
 * Time-stamp:  2003-04-23 08:26:51 +0200
 * E-mail:      pjacklam@online.no
 * URL:         http://home.online.no/~pjacklam
 */
long double ltqnorm(double probability);


/**
 * Gets the probability-P confidence interval for a normal 
 * distribution centered around the mean, in units of sigma.
 *
 * Convenience method for using ltqnorm().
 */
long double computeConfidenceInterval(double probability);


#endif // STAT_UTILS_H
