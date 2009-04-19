#include <algorithm>
#include <iterator>
#include <iostream>
#include <set>
#include <cstdlib>
#include <sstream>
using namespace std;

#include "kmedoids.h"
#include "counter.h"
using namespace cluster;

void makeClusterings(cluster_list& c1, cluster_list& c2) {
  const char *clusters1[] = {
    "1 3 5 7 9",
    "2 4",
    "6 8"
  };

  const char *clusters2[] = {
    "1 3 5 7 9",
    "2 4 11",
    "6 8"
  };
  
  c1.resize(sizeof(clusters1)/sizeof(char*));
  c2.resize(sizeof(clusters2)/sizeof(char*));
  
  for (size_t i=0; i < c1.size(); i++) {
    istringstream s1(clusters1[i]);
    istringstream s2(clusters2[i]);
    copy(istream_iterator<size_t>(s1), istream_iterator<size_t>(), inserter(c1[i], c1[i].begin()));
    copy(istream_iterator<size_t>(s2), istream_iterator<size_t>(), inserter(c2[i], c2[i].begin()));
  }
}


int main(int argc, char **argv) {
  cluster_list c1, c2;
  makeClusterings(c1, c2);

  cout << c1 << endl;
  cout << c2 << endl;
  cout << mirkin_distance(c1, c2) << endl;
}
