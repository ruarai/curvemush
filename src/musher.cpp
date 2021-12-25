#include <random>
#include "musher.h"

// Fast random number generator for picking samples
// from https://stackoverflow.com/a/3747462
static unsigned int g_seed = 1;
inline int fastrand() { 
  g_seed = (214013*g_seed+2531011); 
  return (g_seed>>16)&0x7FFF; 
} 


std::vector<int> musher::mush_curve(int n_start) {
  int n_compartments = 4;
  int n_types = 2;
  
  int n_days = 365;
  int steps_per_day = 4;
  
  int n_steps = n_days * steps_per_day;
  
  int n_array = n_steps * n_compartments * n_types;
  
  auto ix = [n_compartments, n_types] (int t, int compartment, int type) {
    return t * n_compartments * n_types +
      compartment * n_types + 
      type;
  };
  
  
  std::vector<int> arr(n_array);
  
  arr[ix(0, 0, 0)] = n_start;
  arr[ix(0, 0, 1)] = n_start;
  
  int n_delay_samples = 512;
  
  
  std::random_device rd;
  std::mt19937 gen(rd());
  
  std::gamma_distribution<> delay_dist(10,5);
  
  std::vector<int> delay_samples(n_delay_samples);
  for(int i = 0; i < n_delay_samples; i++)
    delay_samples[i] = delay_dist(gen) * steps_per_day;
  
  
  for(int t = 0; t < n_steps; t++) {
    for(int c = 0; c < n_compartments; c++) {
      if(t > 0) {
        arr[ix(t, c, 0)] = arr[ix(t, c, 0)] + arr[ix(t - 1, c, 0)];
      }
      
      // Do we apply delay to this compartment?
      if(c < n_compartments - 1) {
        int n_to_transition = arr[ix(t, c, 1)];
        
        // Surprisingly, performing transitions one-by-one is pretty fast
        // Maybe it gets compiled down to something better?
        for(int i = 0; i < n_to_transition; i++) {
          double delay_sample = delay_samples[fastrand() % n_delay_samples];
          
          int t_set = t + delay_sample;
          
          if(t_set >= n_steps)
            t_set = n_steps - 1;
          
          arr[ix(t_set, c + 1, 0)]++;
          arr[ix(t_set, c + 1, 1)]++;
          
          arr[ix(t_set, c, 0)]--;
        }
      }
    }
  }
  
  
  return arr;
}

