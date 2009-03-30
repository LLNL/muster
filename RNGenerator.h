/* 
 * Copyright (c) 2005  Renaissance Computing Institute. All rights reserved.
 * 
 * This software is open source. See the bottom of this file for the license.
 * 
 * Renaissance Computing Institute, 
 * (A Joint Institute between the University of North Carolina at Chapel Hill,
 * North Carolina State University, and Duke University)
 * http://www.renci.org
 * 
 * For questions, comments please contact software@renci.org
 * 
 */
#ifndef RN_GENERATOR_H
#define RN_GENERATOR_H

#include "cluster/MersenneTwister.h"

/**
 * Wrapper for random number generator that seeds the generator
 * from the system clock, with an optional perturbation, on
 * construction.
 *
 * The user can optionally provide some "salt" to change the seed. 
 * The current time is multiplied by one plus the salt value, and
 * the salt valueis subtracted from the result, and the resulting 
 * number is used as the seed.
 *
 * This could be used, for examlpe, to differentiate generators on 
 * different nodes by their MPI process ids.
 * 
 * This uses the Mersenne Twister generator internally.  For more 
 * info, see MersenneTwister.h
 */
class RNGenerator {
private:
    /** Wrapped mersenne twister generator. */
    MTRand generator;

public:

    /** 
     * Constructs a new random number generator and seeds from current time. 
     * <p>
     * @param salt The current time is multiplied by one plus the salt value,
     * and then the salt valueis subtracted from the result, and the resulting 
     * number is used as the seed.  If this is not provided, it defaults to 1.
     */
    RNGenerator(int salt=1);


    /** @return a random number in [0,1]. */
    double random() {
        return generator.rand();
    }


    /** @return a random integer in [0,2^32-1]. */
    unsigned randomInt() {
        return generator.randInt();
    }

        
    /** @return integer in [0,n] for n < 2^32. */
    unsigned randomInt( const int& n ) {
        return generator.randInt(n);
    }


    /** @return normally distributed random numbers. */
    double randNorm(const double& mean = 0.0, const double& variance = 0.0) {
        return generator.randNorm(mean, variance);
    }


    /**
     * Function call operator : allows this to model an STL RandomNumberGenerator,
     * and it can thus be passed to algorithms like random_sample or random_sample_n.
     *
     * @return an int in [0,N)
     */
    int operator()(int N) {
        return generator.randInt(N-1);
    }


    /**
     * This gets the system value and perturbs it using the value of salt
     * provided.
     * @param salt unique int value used to perturb the system time.
     */
    static int getSeed(int salt);
};

#endif //RN_GENERATOR_H

/* ***************************************************************************

RENCI Open Source Software License
The University of North Carolina at Chapel Hill

The University of North Carolina at Chapel Hill (the "Licensor") through 
its Renaissance Computing Institute (RENCI) is making an original work of 
authorship (the "Software") available through RENCI upon the terms set 
forth in this Open Source Software License (this "License").  This License 
applies to any Software that has placed the following notice immediately 
following the copyright notice for the Software:  Licensed under the RENCI 
Open Source Software License v. 1.0.

Licensor grants You, free of charge, a world-wide, royalty-free, 
non-exclusive, perpetual, sublicenseable license to do the following to 
deal in the Software without restriction, including without limitation the 
rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the 
Software is furnished to do so, subject to the following conditions:

. Redistributions of source code must retain the above copyright notice, 
this list of conditions and the following disclaimers.

. Redistributions in binary form must reproduce the above copyright 
notice, this list of conditions and the following disclaimers in the 
documentation and/or other materials provided with the distribution.

. Neither You nor any sublicensor of the Software may use the names of 
Licensor (or any derivative thereof), of RENCI, or of contributors to the 
Software without explicit prior written permission.  Nothing in this 
License shall be deemed to grant any rights to trademarks, copyrights, 
patents, trade secrets or any other intellectual property of Licensor 
except as expressly stated herein.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE CONTIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.

You may use the Software in all ways not otherwise restricted or 
conditioned by this License or by law, and Licensor promises not to 
interfere with or be responsible for such uses by You.  This Software may 
be subject to U.S. law dealing with export controls.  If you are in the 
U.S., please do not mirror this Software unless you fully understand the 
U.S. export regulations.  Licensees in other countries may face similar 
restrictions.  In all cases, it is licensee's responsibility to comply 
with any export regulations applicable in licensee's jurisdiction.

*************************************************************************** */
