#include <algorithm>
#include <iterator>
#include <iostream>
#include <set>
#include <cstdlib>
#include <sstream>
using namespace std;

#include "KMedoids.h"
#include "counter.h"


void makeClusterings(KMedoids::clusterList& c1, KMedoids::clusterList& c2) {
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
  KMedoids::clusterList c1, c2;
  makeClusterings(c1, c2);

  printClustering(c1);
  cout << endl;

  printClustering(c2);
  cout << endl;
  
  cout << mirkin_distance(c1, c2) << endl;

  exit(1);

  set<size_t> s1;
  set<size_t> s2;

  set<size_t> *target = &s1;
  for (int i=1; i < argc; i++) {
    if (!strcmp(argv[i], ",")) {
      target = &s2;
    } else {
      target->insert(strtol(argv[i], NULL, 0));
    }
  }
  
  copy(s1.begin(), s1.end(), ostream_iterator<size_t>(cout));
  cout << endl;
  copy(s2.begin(), s2.end(), ostream_iterator<size_t>(cout));
  cout << endl;
  
  size_t count;
  set<size_t> real;
  set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), counter<size_t>(count));
  set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(real, real.begin()));
  cout << count << endl;
  cout << real.size() << endl;
}
