#include <iostream>
#include <iterator>
#include <sstream>
#include <algorithm>
#include "partition.h"

using namespace cluster;
using namespace std;

static bool passed = true;

template <class T, class Iterator2>
void validate(T* expected_begin, T* expected_end, Iterator2 actual_begin) {
  if (!equal(expected_begin, expected_end, actual_begin)) {
    passed = false;
    cerr << "expected: ";
    copy(expected_begin, expected_end, ostream_iterator<T>(cerr, " "));
    cerr << endl;
    cerr << "found:    ";
    int len = expected_end - expected_begin;
    copy(actual_begin, actual_begin + len, ostream_iterator<T>(cerr, " "));
    cerr << endl;
  }
}


int main(int argc, char **argv) {
  cluster::partition p;

  // indices         0  1  2  3  4  5  6  7  8  9 10 11 12 13
  medoid_id ids[] = {0, 0, 0, 1, 0, 1, 1, 1, 2, 1, 2, 2, 2, 0};
  const size_t num_ids = sizeof(ids) / sizeof(medoid_id);
  
  object_id medoids[] = {0, 3, 8};
  const size_t num_meds = sizeof(medoids) / sizeof(object_id);
  for (size_t i=0; i < num_meds; i++) {
    p.medoid_ids.push_back(medoids[i]); 
  }
  for (size_t i=0; i < num_ids; i++) {
    p.cluster_ids.push_back(ids[i]);
  }
  

  vector<medoid_id> out_ids;
  p.write_members(0, back_inserter(out_ids));
  medoid_id expected0[] = {0, 1, 2, 4, 13};
  validate(expected0, expected0 + 5, out_ids.begin());
  
  out_ids.clear();
  p.write_members(1, back_inserter(out_ids));
  medoid_id expected1[] = {3, 5, 6, 7, 9};
  validate(expected1, expected1 + 5, out_ids.begin());
  
  out_ids.clear();
  p.write_members(2, back_inserter(out_ids));
  medoid_id expected2[] = {8, 10, 11, 12};
  validate(expected2, expected2 + 4, out_ids.begin());


  ostringstream out;
  out << p.members(0);
  const char *expected = "0-2 4 13";
  validate(expected, expected + strlen(expected), out.str().c_str());
  
  out.str("");
  out << p.members(1);
  expected = "3 5-7 9";
  validate(expected, expected + strlen(expected), out.str().c_str());
  
  out.str("");
  out << p.members(2);
  expected = "8 10-12";
  validate(expected, expected + strlen(expected), out.str().c_str());
  

  if (passed) {
    cerr << "PASSED" << endl;
    exit(0);
  } else {
    cerr << "FAILED" << endl;
    exit(1);
  }
}
