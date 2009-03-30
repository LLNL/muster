#include <iostream>
#include <vector>
using namespace std;

#include "matrix.h"
#include "Point.h"


int main(int argc, char **argv) {
    
    matrix<double> m(16,16);

    Point z(0,0);

    unsigned num = 0;
    for (unsigned i=0; i < 16; i++) {
        for (unsigned j=0; j < 16; j++) {
            Point p(i,j);
            m(i,j) = p.distance(&z);
        }
    }

    num=0;
    for (unsigned i=0; i < 16; i++) {
        for (unsigned j=0; j < 16; j++) {
            //assert(m.get(i,j) == num);
            //assert(m(i,j) == num);
            num++;
        }
    }

    m.printOn(cout);

    
    vector<unsigned> x;
    for (unsigned i=0; i < 15; i+=4) {
        x.push_back(i);
    }

    vector<unsigned> y;
    for (unsigned i=0; i < 15; i+=8) {
        y.push_back(i);
    }

    cerr << endl << " ------- " << endl << endl;

    matrix<double> *subM = m.getSubmatrix(x.begin(), x.end(), x.begin(), x.end());
    subM->printOn(cout);
    delete subM;
}
