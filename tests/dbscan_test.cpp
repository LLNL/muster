//////////////////////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2010, Lawrence Livermore National Security, LLC.  
// Produced at the Lawrence Livermore National Laboratory  
// Written by Juan Gonzalez, juan.gonzalez@bsc.es

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


#include <iostream>
#include <fstream>
#include <cstring>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include "density_based.h"
#include "Timer.h"
#include "point.h"
#include "spherical_clustering_generator.h"

using namespace std;
using namespace cluster;

typedef boost::minstd_rand base_generator_type;

void usage() {
  cerr << "Usage: dbscan-test [-htdn] [-r number_of_points] [-e epsilon] [-m min_points] [-i <data_file> ] [-o output_file]" << endl;
  cerr << "  DBSCAN test case for clustering points" << endl;
  cerr << "Options:" << endl;
  cerr << "  -h         Show this message." << endl;
  cerr << "  -t         Output timing info to file." << endl;
  cerr << "  -d         Verbose debug output." << endl;
  cerr << "  -n         Normalize values of the dimensions" << endl;
  cerr << "  -r         Run using random points" << endl;
  cerr << "  -e         Epsilon parameter of DBSCAN algorithm" << endl;
  cerr << "  -m         MinPoints parameter of DBSCAN algorithm" << endl;
  cerr << "  -i         Load points from data file" << endl;
  cerr << "  -o         Write out the resulting points and the cluster assignment" << endl;

  exit(1);
}

bool debug = false;
bool timing = false;
bool validate = false;

bool   random_run = false;
size_t num_points = 0;

bool normalize_points = false;

const size_t dimensions = 2;

bool   epsilon_set = false;
double epsilon;

bool   min_points_set = false;
size_t min_points;

char* input_file  = NULL;
char* output_file = NULL;

vector<point> points;

/// Uses getopt to read in arguments.
void get_args(int *argc, char ***argv)
{
  int c;
  char *err;

  if (*argc == 0)
  {
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
  
  while ((c = getopt(*argc, *argv, "htdnr:e:m:i:o:")) != -1) {
    switch (c) {
      case 'h':
        usage();
        exit(1);
        break;
      case 't':
        timing = true;
        break;
      case 'd':
        debug = true;
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

  
  if (random_run && input_file != NULL)
  {
    random_run = false;
  }

  if (!epsilon_set || !min_points)
  {
    cerr << "You must set 'epsilon' (-e) and 'min_points' (-m)" << endl;
    exit (EXIT_FAILURE);
  }
}

void load_file()
{
  char line[200];
  
  ifstream input_stream (input_file, ifstream::in);

  if (!input_stream)
  {
    cerr << "Unable to open input file" << endl;
    exit (EXIT_FAILURE);
  }

  while (input_stream.good())
  {
    input_stream.getline(line, 200);
    parse_point_csv (line, points);
  }

  input_stream.close();
}

void random_generator()
{
  base_generator_type generator(42u);
  boost::uniform_real<> uni_dist(0,1);
  boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

  spherical_clustering_generator cgen;

  /* Always 10 spherical clusters are generated */
  cgen.set_default_stddev(uni());
  for (size_t i = 0; i < 10; i++)
  {
    cgen.add_cluster(point(uni(),uni()*4));
  }

  for (size_t i = 0; i < num_points; i++)
  {
    point p;

    do {
      p = cgen.next_point();
    } while (p.x < 0 || p.y < 0);

    points.push_back(p);
  }

  /*
  for (size_t i = 0; i < num_points; i++)
  {
    points.push_back(point(uni(), uni()));
  }
  */

}

void normalize()
{
  for (size_t i = 0; i < points.size(); i++)
  {
    points[i].normalize();
  }
}

void flush_points(density_based& clustering_results)
{
  ofstream output_stream (output_file, ios_base::trunc);

  if (!output_stream)
  {
    cerr << "Error opening output file!" << endl;
    exit (EXIT_FAILURE);
  }
  
  cout.setf(std::ios::fixed);
  for (size_t i = 0; i < points.size(); i++)
  {
    output_stream.setf(std::ios::fixed);

    output_stream << points[i].x << ", " << points[i].y << ", " << clustering_results.get_cluster(i) << endl;
    // output_stream << clustering_results;
  }
}



int main(int argc, char **argv)
{
  Timer timer;
  
  density_based clustering = density_based();
  
  get_args(&argc, &argv);

  if (random_run)
  {
    cout << "Generating " << num_points << " random points...";
    
    timer.record("Generation of random points starts");
    random_generator();
    timer.record("Genaration of random points ends");

    cout << " READY" << endl;
  }
  else
  {
    cout << "Loading input file...";
    
    timer.record("File load starts");
    load_file();
    timer.record("File load ends");

    cout << " READY" << endl;
  }

  // Points are available in 'points' vector
  
  cout << "Clustering Points...";
  
  timer.record("Clustering starts");
  clustering.dbscan(points, point_distance(), epsilon, min_points);
  timer.record("Clustering ends");

  cout << " READY (Clusters found = " << clustering.num_clusters() << ")" << endl;

  if (output_file != NULL)
  {
    timer.record("Output file write starts");
    flush_points(clustering);
    timer.record("Output file write ends");
  }

  // Write timing
  if (timing)
  {
    cout << endl;
    cout << "**** EXECUTION TIMING ****" << endl;
    timer.write(cout);
  }
  
  exit (EXIT_SUCCESS);
}

