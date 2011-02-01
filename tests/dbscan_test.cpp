//////////////////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2010, Lawrence Livermore National Security, LLC.  
// Produced at the Lawrence Livermore National Laboratory  
// LLNL-CODE-433662
// All rights reserved.  
//
// This file is part of Muster. For details, see http://github.com/tgamblin/muster. 
// Please also read the LICENSE file for further information.
//
// Redistribution and use in source and binary forms, with or without modification, are
// permitted provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright notice, this list of
//    conditions and the disclaimer below.
//  * Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the disclaimer (as noted below) in the documentation and/or other materials
//    provided with the distribution.
//  * Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse
//    or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
// OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
// LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//////////////////////////////////////////////////////////////////////////////////////////////////

///
/// @file par_cluster_test.cpp
/// @author Juan Gonzalez juan.gonzalez@bsc.es
/// 
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>

#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/basic.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/Convex_hull_d_traits_3.h>
#include <CGAL/algorithm.h>

#include "Timer.h"
#include "point.h"
#include "point_set.h"
#include "spherical_clustering_generator.h"
#include "cdbw.h"
#include "density.h"

using namespace CGAL;
using namespace cluster;
using namespace std;

typedef Convex_hull_d_traits_3<Exact_predicates_inexact_constructions_kernel>  K;
typedef K::Point_d P3;

typedef boost::minstd_rand base_generator_type;

/// Distance bt/w two cgal points
struct cgal_2d_distance {
  double operator()(const P3& left, const P3& right) const {
    double dx = left[0] - right[0];
    double dy = left[1] - right[1];
    return ::sqrt(dx*dx + dy*dy);
  }
};


void usage() {
  cerr << "Usage: dbscan-test [-htvc] "
       << "[-r <number_of_points>] [-i <data_file>] [-o output_file] "
       << "-e <epsilon> -m <min_points>" 
       << endl;
  cerr << "  DBSCAN test case for clustering points" << endl;
  cerr << "Options:" << endl;
  cerr << "  -h         Show this message." << endl;
  cerr << "  -c         Use CDBW to evaluate the clustering." << endl;
  cerr << "  -t         Output timing info to file." << endl;
  cerr << "  -v         Verbose output." << endl;
  cerr << "  -r         Run using random points" << endl;
  cerr << "  -e         Epsilon parameter of DBSCAN algorithm" << endl;
  cerr << "  -m         MinPoints parameter of DBSCAN algorithm" << endl;
  cerr << "  -i         Load points from data file" << endl;
  cerr << "  -s         Use sampled density clustering." << endl;
  cerr << "  -o         Write out the resulting points and the cluster assignment" << endl;

  exit(1);
}

static bool   verbose = false;
static bool   timing = false;
static bool   random_run = false;
static size_t num_points = 0;
static bool   normalize_points = false;

static bool   epsilon_set = false;
static double epsilon;

static bool   min_points_set = false;
static size_t min_points;

static bool   sampling  = false;
static bool   cdbw = false;

static char*  input_file  = NULL;
static char*  output_file = NULL;

static point_set points;


/// Uses getopt to read in arguments.
void get_args(int *argc, char ***argv) {
  int c;
  char *err;

  if (*argc == 0) {
    cout << "Execution of 10,000 random points using epsilon = 0.01 and min_points = 4" << endl;
    cout << "Results written in output file 'dbscan_random_execution.csv'" << endl;

    random_run = true;
    num_points = 10000;

    epsilon_set = true;
    epsilon     = 0.01;

    min_points_set = true;
    min_points     = 4;

    output_file = strdup("dbscan_random_execution.csv");
  }
  
  while ((c = getopt(*argc, *argv, "htdncsr:e:m:i:o:")) != -1) {
    switch (c) {
      case 'h':
        usage();
        exit(1);
        break;
      case 't':
        timing = true;
        break;
      case 'c':
        cdbw = true;
        break;
      case 'v':
        verbose = true;
        break;
      case 's':
        sampling = true;
        break;
      case 'n':
        normalize_points = true;
        break;
      case 'r':
        random_run = true;
        num_points = strtol(optarg, &err, 0);
        if (*err)
          usage();
        break;
      case 'e':
        epsilon_set = true;
        epsilon     = strtod(optarg, &err);
        if (*err) usage();
        break;
      case 'm':
        min_points_set = true;
        min_points     = strtol(optarg, &err, 0);
        if (*err) usage();
        break;
      case 'i':
        input_file = strdup(optarg);
        break;
      case 'o':
        output_file = strdup(optarg);
        break;
      default:
        usage();
        exit(1);
        break;
    }
  }

  // adjust params
  *argc -= optind;
  *argv += optind;

  
  if (random_run && input_file != NULL) {
    random_run = false;
  }

  if (!epsilon_set || !min_points) {
    cerr << "You must provide -e and -m for epsilon and min_points." << endl;
    usage();
    exit(EXIT_FAILURE);
  }
}


void generate_random_points() {
  base_generator_type generator(42u);
  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

  spherical_clustering_generator cgen;

  /* Always 10 spherical clusters are generated */
  cgen.set_default_stddev(uni());
  for (size_t i = 0; i < 10; i++) {
    cgen.add_cluster(point(uni(),uni()*4));
  }

  for (size_t i = 0; i < num_points; i++) {
    point p;

    do {
      p = cgen.next_point();
    } while (p.x < 0 || p.y < 0);

    points.add_point(p);
  }
}


int main(int argc, char **argv) {
  Timer timer;
  density clustering;
  get_args(&argc, &argv);

  timer.record("Initialization");

  if (random_run) {
    cout << "Generating " << num_points << " random points...";
    generate_random_points();
    timer.record("Genarate Random Points");
    cout << " READY" << endl;

  } else {
    cout << "Loading input file...";
    ifstream csv_file(input_file, ifstream::in);
    points.load_csv_file(csv_file);
    timer.record("Load Input File");
    cout << " READY" << endl;
  }
  
  cout << "Converting point_set to CGAL points...";
  vector<P3> cgal_points;
  for (size_t i=0; i < points.size(); i++) {
    cgal_points.push_back(P3(points[i][0], points[i][1], 0));
  }

  // Points are available in 'cgal_points' vector
  cout << "Clustering Points...";

  if (!sampling) {
    clustering.dbscan(cgal_points, cgal_2d_distance(), epsilon, min_points);
  } else {
    clustering.sdbscan<K>(cgal_points, cgal_2d_distance(), epsilon, min_points);
  }
  timer.record("Clustering");

  cout << " READY (Clusters found = " << clustering.num_clusters() << ")" << endl;

  if (output_file != NULL) {
    ofstream csv_file(output_file, ios_base::trunc);
    if (!csv_file) {
      cerr << "Error opening output file!" << endl;
      exit (EXIT_FAILURE);
    }
    points.write_csv_file(csv_file, &clustering);
    timer.record("Write Output File");
  }

  if (cdbw) {
    // Checking the CDbw on everything but the noise.
    clustering.remove_cluster(density::NOISE);
    CDbw<P3, cgal_2d_distance> validation(clustering, cgal_points);
    double current_cdbw = validation.compute(10); // 10 representatives per cluster
    timer.record("CDBW");

    cout << "CDbw = " << current_cdbw << endl;
    cout << "  separtion:    " << validation.separation() << endl;
    cout << "  compactness:  " << validation.compactness() << endl;
    cout << "  cohesion:     " << validation.cohesion() << endl;
    cout << "  sep*compact:  " << validation.separation() * validation.compactness() << endl;
  }
  
  // Write timing
  if (timing) {
    cout << endl;
    cout << "**** EXECUTION TIMING ****" << endl;
    timer.write(cout);
  }

  exit(EXIT_SUCCESS);
}
