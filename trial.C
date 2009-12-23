#include "trial.h"

#include <algorithm>
using namespace std;

namespace cluster {

  size_t trial_generator::get_sample_size(size_t k) {
    return min(init_size + 2 * k, max_sample);
  }
    
  trial_generator::trial_generator(size_t _max_k, size_t _max_reps, size_t _init_size, size_t _max_sample)
    : max_k(_max_k), 
      max_reps(_max_reps), 
      init_size(_init_size), 
      max_sample(_max_sample),
      cur_trial(1, 0, get_sample_size(1)),
      iterations(0)
  { }
    
    
  trial_generator::trial_generator(size_t _min_k, size_t _max_k, size_t _max_reps, size_t _init_size, 
                                 size_t _max_sample)
    : max_k(_max_k), 
      max_reps(_max_reps), 
      init_size(_init_size), 
      max_sample(_max_sample),
      cur_trial(_min_k, 0, get_sample_size(_min_k)),
      iterations(0)
  { }
    

  bool trial_generator::has_next() const {
    return cur_trial.k <= max_k;
  }


  trial trial_generator::next() {
    trial result = cur_trial;

    // iterate first through trials, then through k's, and only do one trial if the sample size is
    // equal to the maximum sample size.
    cur_trial.rep++;
    if (cur_trial.rep >= max_reps || cur_trial.sample_size == max_sample) {
      cur_trial.rep = 0;
      cur_trial.k++;
      cur_trial.sample_size = get_sample_size(cur_trial.k);
    }

    iterations++;
    return result;
  }


  size_t trial_generator::count() const {
    return iterations;
  }



} // namespace cluster
