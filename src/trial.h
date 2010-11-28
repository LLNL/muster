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
/// @file trial.h
/// @author Todd Gamblin tgamblin@llnl.gov
/// @brief Data structure representing a trial run of a partitioned clustering algorithm.
///
#ifndef TRIAL_H
#define TRIAL_H

#include <cstdlib>

namespace cluster {

  ///
  /// This struct represents parameters for a single trial run of kmedoids.
  /// We generate a bunch of these to farm all the trials out to processes when doing
  /// parallel clustering.
  /// 
  /// Each trial has a particular <i>k</i> (number of clusters), sample size, and 
  /// repetition number (for successive samples with the same k).
  ///
  struct trial {
    size_t k;
    size_t rep;
    size_t sample_size;

    trial() { }
    
    trial(size_t _k, size_t _rep, size_t _sample_size)
      : k(_k), rep(_rep), sample_size(_sample_size) { }
    
    trial(const trial& other)
      : k(other.k), rep(other.rep), sample_size(other.sample_size) { }
  };
  
  ///
  /// Class to generate a set of trials for clustering.  This packages up state for the main
  /// trial loop and allows work to be farmed out to worker processes.
  ///
  class trial_generator {    
  public:
    ///
    /// Constructor to generate trials from 1 to max_k.
    /// 
    trial_generator(size_t _max_k, size_t _max_reps, size_t _init_size, size_t _num_objects);

    ///
    /// Constructor to generate trials from min_k to max_k.
    ///
    trial_generator(size_t min_k, size_t _max_k, 
                    size_t _max_reps, size_t _init_size, size_t _num_objects);
    
    size_t count() const;       ///< return iterations so far.
    bool has_next() const;      ///< whether there are trials remaining.
    trial next();               ///< return parameters for next trial
    void reset();               ///< return to initial state
    size_t num_trials();        ///< Return total number of trials this will generate.

    const size_t max_k;         ///< maximum k to try
    const size_t max_reps;      ///< max number of repetitions per k
    const size_t init_size;     ///< initial size for samples before factoring in k, as per CLARA paper.
    const size_t num_objects;   ///< number of elements in the data set; determines maximum sample size.

  private:
    trial cur_trial;            ///< current state of the iterator.
    size_t number_of_trials;    ///< memoized total number of trials in this generator
    size_t iterations;          ///< number of iterations so far
    
    /// size of the sample to cluster for particular k
    size_t get_sample_size(size_t k);
  };
  
} // namespace cluster

#endif // TRIAL_H
