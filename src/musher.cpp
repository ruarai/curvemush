#include <RcppThread.h>

#include <random>
#include "musher.h"


#define s_occupancy 0
#define s_transitions 1

static int g_seed = 0;

// Fast random number generator for picking samples
// from https://stackoverflow.com/a/3747462
inline int fastrand() { 
  g_seed = (214013*g_seed+2531011); 
  return (g_seed>>16)&0x7FFF; 
} 

// mush_curve
// Performs the discrete-time stochastic simulation of the compartmental model
std::vector<int> musher::mush_curve(mush_params params) {
  int n_types = 2;
  
  int n_steps = params.n_days * params.steps_per_day;
  
  int n_compartments = params.n_compartments;
  int n_array = n_steps * n_compartments * n_types;
  
  // Helper lambda function to index the data array
  auto ix = [n_compartments, n_types] (int t, int compartment, int type) {
    return t * n_compartments * n_types +
      compartment * n_types + 
      type;
  };
  
  
  std::vector<int> arr(n_array);
  
  int n_start = 1000;
  arr[ix(0, 0, 0)] = n_start;
  arr[ix(0, 0, 1)] = n_start;
  
  
  std::random_device rd;
  std::mt19937 rand(rd());
  
  std::vector<int> symptomatic_ward_delays = musher::make_delay_samples(
    params.n_delay_samples, 1.582, 4.589, rand
  );
  std::vector<int> ward_to_discharge_delays = musher::make_delay_samples(
    params.n_delay_samples, 0.907, 6.299, rand
  );
  std::vector<int> ward_to_ICU_delays = musher::make_delay_samples(
    params.n_delay_samples, 0.515, 4.560, rand
  );
  std::vector<int> ward_to_death_delays = musher::make_delay_samples(
    params.n_delay_samples, 1.931, 5.670, rand
  );
  
  
  int n_delay_samples = params.n_delay_samples;
  
  for(int t = 0; t < n_steps; t++) {
    // Update the current occupancy value to reflect the value at the last timestep
    // counted at arr[ix(t - 1, c, s_occupancy)]
    // plus/minus any transitions (counted by arr[ix(t, c, s_occupancy)])
    for(int c = 0; c < n_compartments; c++) {
      if(t > 0)
        arr[ix(t, c, s_occupancy)] = arr[ix(t, c, s_occupancy)] + arr[ix(t - 1, c, s_occupancy)];
    }
    
    
    
    musher::transition_symptomatic_ward(t, arr, ix, 
                                        symptomatic_ward_delays,
                                        n_steps);
    
    musher::transition_ward_next(t, arr, ix,
                                 ward_to_discharge_delays,
                                 ward_to_ICU_delays,
                                 ward_to_death_delays,
                                 0.1, 0.2,
                                 n_steps, rand);
  }       
  
  
  return arr;
}


template <typename F>
void musher::transition_symptomatic_ward(
    int t,
    std::vector<int> &arr,
    F ix,
    std::vector<int> &symptomatic_ward_delays,
    int n_steps
) {
  int n_to_transition = arr[ix(t, 0, s_transitions)];
  
  transition_generic_delayed(
    0, 1, n_to_transition,
    t, arr, ix,
    symptomatic_ward_delays,
    n_steps
  );
}


template <typename F>
void musher::transition_ward_next(
    int t,
    std::vector<int> &arr,
    F ix,
    
    std::vector<int> &ward_to_discharge_delays,
    std::vector<int> &ward_to_ICU_delays,
    std::vector<int> &ward_to_death_delays,
    float pr_ward_to_death,
    float pr_ward_to_ICU,
    
    int n_steps,
    
    std::mt19937 &rand
) {
  int n_to_transition = arr[ix(t, 1, s_transitions)];
  
  int n_to_discharge = 0;
  int n_to_ICU = 0;
  int n_to_death = 0;
  
  
  std::uniform_real_distribution<> dist(0, 1);
  
  // Need a faster way to do multinomial sampling, probably
  // but for now, this is faster than it looks
  for(int i = 0; i < n_to_transition; i++) {
    float rand_sample = dist(rand);
    
    if(rand_sample < pr_ward_to_death) {
      n_to_death++;
    } else if(rand_sample < pr_ward_to_death + pr_ward_to_ICU) {
      n_to_ICU++;
    } else{
      n_to_discharge++;
    }
  }
  
  
  transition_generic_delayed(
    1, 2, n_to_discharge,
    t, arr, ix,
    ward_to_discharge_delays,
    n_steps
  );
  transition_generic_delayed(
    1, 3, n_to_ICU,
    t, arr, ix,
    ward_to_ICU_delays,
    n_steps
  );
  transition_generic_delayed(
    1, 4, n_to_death,
    t, arr, ix,
    ward_to_death_delays,
    n_steps
  );
  
}

template <typename F>
void musher::transition_generic_delayed(
    int c_from, int c_to,
    int n_to_transition,
    int t,
    std::vector<int> &arr,
    F ix,
    std::vector<int> &delay_samples,
    int n_steps
) {
  
  int n_delay_samples = delay_samples.size();
  
  // Perform our transitions one at a time. This is surprisingly fast!
  for(int i = 0; i < n_to_transition; i++) {
    
    // How far forwards do we want to go?
    int delay_sample = delay_samples[fastrand() % n_delay_samples];
    
    // Add 1 so we never transition instantaneously (introducing a slight ~0.5 time step error)
    int t_set = delay_sample + t + 1;
    
    // Dropping transitions that occur outside our simulation
    if(t_set >= n_steps)
      continue;
    
    if(ix(t_set, c_to, s_transitions) >= 5 * 100 * 2) {
      RcppThread::Rcout << ix(t_set, c_to, s_transitions) << " ";
      RcppThread::Rcout << "(" << t_set << ", " << c_to << ", " << s_transitions << ") ";
    }
    
    // We need to handle this transition in the future
    arr[ix(t_set, c_to, s_transitions)]++;
    
    // Indicate that occupancy in the future for our new compartment increases
    arr[ix(t_set, c_to, s_occupancy)]++; 
    
    // Indicate that occupancy in the future for our current compartment decreases
    arr[ix(t_set, c_from, s_occupancy)]--;
  }
}


std::vector<int> musher::make_delay_samples(
    int n_samples, double shape, double scale, std::mt19937 &rand
  ) {
  std::vector<int> samples(n_samples);
  
  std::gamma_distribution<> delay_dist(shape, scale);
  
  for(int i = 0; i < n_samples; i++)
    samples[i] = std::floor(delay_dist(rand));
  
  return samples;
}



