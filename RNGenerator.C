#include "RNGenerator.h"

#include <cstdlib>
#include <sys/time.h>


RNGenerator::RNGenerator(int salt) {
    //seed random number generator with perturbed system clock 
    seed(salt);
}


void RNGenerator::seed(int salt) {
    struct timeval time;
    gettimeofday(&time, 0);
    generator.seed(time.tv_sec * (salt+1) - salt);
}

