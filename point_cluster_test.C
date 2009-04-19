#include <cmath>
#include <vector>
#include <iostream>
#include <cstdlib>
using namespace std;

#include <boost/numeric/ublas/matrix.hpp>
using boost::numeric::ublas::matrix;

#include "dissimilarity.h"
#include "kmedoids.h"
#include "color.h"
using namespace cluster;

const size_t num_colors = 15;
const char *colors[num_colors] = {
  Blue, Green, Cyan, Red, Purple, Brown, Light_Gray, 
  Light_Blue, Light_Green, Light_Cyan, Light_Red, Light_Purple, 
  White, Dark_Gray, Yellow
};


/// Simple 1 dimensional point class for testing medoids algorithms.
struct point {
public:
  int x, y;
    
  /// New point with position (x,y)
  point(int _x, int _y) : x(_x), y(_y) { }
  point(const point& other) : x(other.x), y(other.y) { }

  // Distance between this point and another. 
  double distance(const point& other) const {
    int dx = other.x - x;
    int dy = other.y - y;
    return sqrt((dx*dx) + (dy*dy));
  }
  
  point& operator+=(const point& other) {
    x += other.x;  y += other.y;
    return *this;
  }

  point operator+(const point& other) const {
    point result = *this;
    result += other;
    return result;
  }
  
  point& operator*=(const point& other) {
    x *= other.x;  y *= other.y;
    return *this;
  }

  point operator*(const point& other) const {
    point result = *this;
    result *= other;
    return result;
  }
  
  point& operator=(const point& other) { 
    x = other.x;  y = other.y;
    return *this;
  }
};


void draw(string label,  vector<point>& points, const cluster::partition& parts) {
  int max_x=0;
  int max_y=0;
  for (size_t i=0; i < points.size(); i++) {
    max_x = max(points[i].x, max_x);
    max_y = max(points[i].y, max_y);
  }
  
  matrix<int> pmat(max_x+1, max_y+1);
  for (size_t i=0; i < pmat.size1(); i++) {
    for (size_t j=0; j < pmat.size2(); j++) {
      pmat(i,j) = -1;
    }
  }

  for (size_t i=0; i < points.size(); i++) {
    pmat(points[i].x, points[i].y) = i;
  }
  
  const int min_hbar = 80;
  cout << label << " ";
  for (int x = 0; x <= max(max_x, min_hbar)-label.size()-1; x++) cout << "-";
  cout << endl;

  for (int y = max_y; y >= 0; y--) {
    for (int x = 0; x <= max_x; x++) {
      if (pmat(x,y) >= 0) {
        cout << colors[parts.cluster_ids[pmat(x,y)] % num_colors];
        cout << "o";
        cout << None;
      } else {
        cout << " ";
      }
    }
    cout << endl;
  }

  for (int x = 0; x <= max(max_x, min_hbar); x++) cout << "-";
  cout << endl;
}




/// Distance bt/w two points
struct point_distance {
  double operator()(const point& left, const point& right) {
    return left.distance(right);
  }
};


/// Test using points instead of QGrams to make sure clustering 
/// algorithms work.
int main(int argc, char **argv) { 
    //vector of test points
    vector<point> points;

    size_t clusters = 5;
    if (argc > 1) {
      clusters = strtol(argv[1], NULL, 0);    
    }

    //put 5-pt crosses inthe vector, offset by increasing distances
    point middle(1,1);
    point up(0,1), down(0,-1), left(-1,0), right(1,0);
    for (size_t i=0; i < clusters; i++) {
      points.push_back(middle);
      points.push_back(middle + up);
      points.push_back(middle + down);
      points.push_back(middle + left);
      points.push_back(middle + right);
      
      middle += (point(2*i+4, 0));
    }


    dissimilarity_matrix distance;
    build_dissimilarity_matrix(points, point_distance(), distance);

    kmedoids km;
    kmedoids clara;

    for (int k = 1; k <= clusters; k++) {
      km.pam(distance, k);
      clara.clara(points, point_distance(), k);

      cout << "k: " << k << ", Mirkin distance: " << mirkin_distance(km, clara) << endl;
      draw("PAM", points, km);
      draw("CLARA", points, clara);
      cout << endl;
    }
}
