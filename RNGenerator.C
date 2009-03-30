#include "RNGenerator.h"

#include <cstdlib>
#include <sys/time.h>


RNGenerator::RNGenerator(int salt) :
    //seed random number generator with perturbed system clock 
    generator(getSeed(salt))
{
    ; //nothing to do.
}


int RNGenerator::getSeed(int salt) {
    struct timeval time;
    gettimeofday(&time, 0);
    return (time.tv_sec * (salt+1) - salt);
}

