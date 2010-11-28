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
/// @file trial.cpp
/// @author Todd Gamblin tgamblin@llnl.gov
///
#include "trial.h"

#include <algorithm>
using namespace std;

namespace cluster {

  size_t trial_generator::get_sample_size(size_t k) {
    return min(init_size + 2 * k, num_objects);
  }

  static size_t count_all(const trial_generator& tg) {
    trial_generator dummy = tg;
    while (dummy.has_next()) dummy.next();
    return dummy.count();
  }

    
  trial_generator::trial_generator(size_t _max_k, size_t _max_reps, size_t _init_size, size_t _num_objects)
    : max_k(_max_k), 
      max_reps(_max_reps), 
      init_size(_init_size), 
      num_objects(_num_objects),
      cur_trial(1, 0, get_sample_size(1)),
      iterations(0)
  { 
    number_of_trials = count_all(*this);
  }
    
    
  trial_generator::trial_generator(size_t _min_k, size_t _max_k, size_t _max_reps, size_t _init_size, 
                                 size_t _num_objects)
    : max_k(_max_k), 
      max_reps(_max_reps), 
      init_size(_init_size), 
      num_objects(_num_objects),
      cur_trial(_min_k, 0, get_sample_size(_min_k)),
      iterations(0)
  { 
    number_of_trials = count_all(*this);
  }
    

  bool trial_generator::has_next() const {
    return cur_trial.k <= max_k;
  }


  trial trial_generator::next() {
    trial result = cur_trial;

    // iterate first through trials, then through k's, and only do one trial if the sample size is
    // equal to the maximum sample size.
    cur_trial.rep++;
    if (cur_trial.rep >= max_reps || cur_trial.sample_size == num_objects) {
      cur_trial.rep = 0;
      cur_trial.k++;
      cur_trial.sample_size = get_sample_size(cur_trial.k);
    }

    iterations++;
    return result;
  }


  size_t trial_generator::num_trials() {
    return number_of_trials;    
  }

  size_t trial_generator::count() const {
    return iterations;
  }



} // namespace cluster
