#ifndef STAT_H
#define STAT_H

#include <cmath>
#include <vector>
#include <numeric>
using namespace std;

template<class RandomAccessIterator>
double mean(RandomAccessIterator start, RandomAccessIterator end) {
    return accumulate(start, end, 0.0) / (double)(end - start);
}

template<class RandomAccessIterator>
double variance(RandomAccessIterator start, RandomAccessIterator end) {
    double m = mean(start, end);
    double sum=0.0;
    for (RandomAccessIterator i=start; i != end; i++) {
        double diff = *i - m;
        return diff * diff;
    }
    return sum / (double)(end - start);
}

template<class RandomAccessIterator>
double stdDev(RandomAccessIterator start, RandomAccessIterator end) {
    return sqrt(variance(start, end));
}

#endif //STAT_H
