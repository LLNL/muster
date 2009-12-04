#include "par_kmedoids.h"

#include <cstdlib>
#include <stdint.h>
#include <sys/time.h>
using namespace std;

namespace cluster {

  size_t trial_iterator::get_sample_size(size_t k) {
    return min(init_size + 2 * k, max_sample);
  }
    
  trial_iterator::trial_iterator(size_t _max_k, size_t _max_trials, size_t _init_size, size_t _max_sample)
      : max_k(_max_k), max_trials(_max_trials), init_size(_init_size), max_sample(_max_sample)
  {
    cur_trial.k = 1;
    cur_trial.trial = 0;
    cur_trial.sample_size = get_sample_size(1);
  }
    
    
  trial_iterator::trial_iterator(size_t _min_k, size_t _max_k, size_t _max_trials, size_t _init_size, size_t _max_sample)
    : max_k(_max_k), max_trials(_max_trials), init_size(_init_size), max_sample(_max_sample)
  {
    cur_trial.k = _min_k;
    cur_trial.trial = 0;
    cur_trial.sample_size = get_sample_size(_min_k);
  }
    

  bool trial_iterator::has_next() {
    return cur_trial.k <= max_k;
  }


  trial trial_iterator::next() {
    trial result = cur_trial;

    // iterate first through trials, then through k's, and only do one trial if the sample size is
    // equal to the maximum sample size.
    cur_trial.trial++;
    if (cur_trial.trial >= max_trials || cur_trial.sample_size == max_sample) {
      cur_trial.trial = 0;
      cur_trial.k++;
      cur_trial.sample_size = get_sample_size(cur_trial.k);
    }
    
    return result;
  }


  void par_kmedoids::seed_random_uniform(MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    uint32_t seed = 0;
    if (rank == 0) {
      struct timeval time;
      gettimeofday(&time, 0);
    
      seed = time.tv_sec * time.tv_usec;
    }

    MPI_Bcast(&seed, 1, MPI_INT, 0, comm);
    random.seed(seed);
  }

} // namespace cluster
