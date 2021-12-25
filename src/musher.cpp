#include <random>
#include "musher.h"

#define s_occupancy 0
#define s_transitions 1

// Fast random number generator for picking samples
// from https://stackoverflow.com/a/3747462
static unsigned int g_seed = 1;
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
  std::mt19937 gen(rd());
  
  std::gamma_distribution<> delay_dist(2,1);
  
  int n_delay_samples = params.n_delay_samples;
  std::vector<int> delay_samples(n_delay_samples);
  for(int i = 0; i < n_delay_samples; i++)
    delay_samples[i] = delay_dist(gen) * params.steps_per_day;
  
  
  for(int t = 0; t < n_steps; t++) {
    for(int c = 0; c < n_compartments; c++) {
      // Update the current occupancy value to reflect the value at the last timestep
      // plus/minus any transitions (counted by arr[ix(t, c, s_occupancy)])
      if(t > 0)
        arr[ix(t, c, s_occupancy)] = arr[ix(t, c, s_occupancy)] + arr[ix(t - 1, c, s_occupancy)];
        
      
      // Do we apply delay to this compartment?
      if(c < n_compartments - 1) {
        int n_to_transition = arr[ix(t, c, s_transitions)];
        
        // Perform our transitions one at a time. This is surprisingly fast!
        for(int i = 0; i < n_to_transition; i++) {
          
          // How far forwards do we want to go?
          int delay_sample = delay_samples[fastrand() % n_delay_samples];
          
          int t_set = delay_sample + t;
          
          // Dropping transitions that occur outside our simulation
          if(t_set > n_steps)
            continue;
          
          // We need to handle this transition in the future
          arr[ix(t_set, c + 1, s_transitions)]++;
          
          // Indicate that occupancy in the future for our new compartment increases
          arr[ix(t_set, c + 1, s_occupancy)]++; 
          
          // Indicate that occupancy in the future for our current compartment decreases
          arr[ix(t_set, c, s_occupancy)]--;
        }
      }
    }
  }
  
  
  return arr;
}

