#include <list>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/basic.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/algorithm.h>

using namespace CGAL;
using namespace std;

typedef Convex_hull_d_traits_3<Exact_predicates_inexact_constructions_kernel>  K;
typedef K::Point_d                                     P3;
typedef Convex_hull_d<K>::Vertex_handle                Vertex_handle;

int main() {
  Random_points_in_sphere_3<P3> gen(100.0);
  list<P3> points;

  // generate 250 points randomly on a sphere of radius 100.0
  // and insert them into the triangulation
  CGAL::copy_n(gen, 250, back_inserter(points));
  Convex_hull_d<K> hull(3);
  hull.insert(points.begin(), points.end());
  hull.is_valid(true);

  list<Vertex_handle> vertices(hull.all_vertices());
  cout << "This convex hull of the 250 points has "
       << hull.number_of_vertices()  << " vertices," << endl
       << hull.number_of_facets()    << " facets," << endl
       << hull.number_of_simplices() << " simplices." << endl;

  // start at origin and move out along z axis, printing 
  // whether we're still inside the hull
  for (size_t x=0; x < 200; x+= 10) {
    K::Point_d p(0,0,x);
    cout << p << "\t" << hull.bounded_side(p) << endl;
  }

  return 0;
}
