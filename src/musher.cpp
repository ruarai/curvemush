#include <RcppThread.h>

#include <random>
#include "musher.h"


#define def_n_slots 2

#define s_occupancy 0
#define s_transitions 1

// This ends up being shared across threads - but believe it won't cause issues
// Data race in this context doesn't really matter
static int g_seed = 0;

// Fast random number generator for picking samples
// from https://stackoverflow.com/a/3747462
inline int fastrand() { 
  g_seed = (214013*g_seed+2531011); 
  return (g_seed>>16)&0x7FFF; 
} 

// mush_curve
// Performs the discrete-time stochastic simulation of the compartmental model
mush_results musher::mush_curve(
    mush_params params,

    std::vector<int> case_curve,

    strat_data s_data,

    float scale_los,

    std::vector<float> pr_ICU_curve
  ) {
  int n_steps = params.n_days * params.steps_per_day;
  
  int n_array = n_steps * def_n_compartments * def_n_slots;
  
  // Helper lambda function to index the data array
  auto ix = [] (int t, int compartment, int slot) {
    return t * def_n_compartments * def_n_slots +
      compartment * def_n_slots + 
      slot;
  };
  
  
  std::vector<int> arr(n_array, 0);
  
  for(int d = 0; d < case_curve.size(); d++) {
    //arr[ix(d * params.steps_per_day, c_symptomatic, s_occupancy)] = case_curve[d];
    arr[ix(d * params.steps_per_day, c_symptomatic, s_transitions)] = case_curve[d];
  }
  
  std::random_device rd;
  std::mt19937 rand(rd());
  
  std::vector<int> symptomatic_ward_delays = musher::make_delay_samples(
    params.n_delay_samples, std::exp(std::log(s_data.d_shape_symptomatic_to_ward) + scale_los), s_data.d_scale_symptomatic_to_ward,
    params.steps_per_day, rand);
  
  std::vector<int> ward_to_discharge_delays = musher::make_delay_samples(
    params.n_delay_samples, std::exp(std::log(s_data.d_shape_ward_to_discharge) + scale_los), s_data.d_scale_ward_to_discharge,
    params.steps_per_day, rand);
  std::vector<int> ward_to_ICU_delays = musher::make_delay_samples(
    params.n_delay_samples, std::exp(std::log(s_data.d_shape_ward_to_ICU) + scale_los), s_data.d_scale_ward_to_ICU,
    params.steps_per_day, rand);
  std::vector<int> ward_to_death_delays = musher::make_delay_samples(
    params.n_delay_samples, std::exp(std::log(s_data.d_shape_ward_to_death) + scale_los), s_data.d_scale_ward_to_death,
    params.steps_per_day, rand);
  
  
  std::vector<int> ICU_to_discharge_delays = musher::make_delay_samples(
    params.n_delay_samples, s_data.d_shape_ICU_to_discharge, s_data.pr_ICU_to_discharge,
    params.steps_per_day, rand);
  std::vector<int> ICU_to_death_delays = musher::make_delay_samples(
    params.n_delay_samples, s_data.d_shape_ICU_to_death, s_data.d_scale_ICU_to_death,
    params.steps_per_day, rand);
  std::vector<int> ICU_to_postICU_delays = musher::make_delay_samples(
    params.n_delay_samples, s_data.d_shape_ICU_to_postICU, s_data.d_scale_ICU_to_postICU,
    params.steps_per_day, rand);
  
  std::vector<int> postICU_to_discharge_delays = musher::make_delay_samples(
    params.n_delay_samples, std::exp(std::log(s_data.d_shape_postICU_to_discharge) + scale_los), s_data.d_scale_postICU_to_discharge, 
    params.steps_per_day, rand);
  std::vector<int> postICU_to_death_delays = musher::make_delay_samples(
    params.n_delay_samples, std::exp(std::log(s_data.d_shape_postICU_to_death) + scale_los), s_data.d_scale_postICU_to_death,
    params.steps_per_day, rand);
  
  int n_delay_samples = params.n_delay_samples;
  
  // The primary simulation loop:

  // We move across t = [0 .. n_steps - 1], reading from n = arr[ix(t, c_current, s_transitions)] which indicates
  // how many individuals have /entered/ compartment c_current at timestep t. For each individual i = [0 .. n - 1],
  // we select their next compartment c_next and their delay d. 
  // We then set arr[ix(t + d + 1, c_next, s_transitions)]++, such that their next transition will be handled
  // at timestep t + d + 1.

  // To know the number of individuals in each compartment across time, we maintain a second value (in 'slot'
  // s_occupancy), which is a cumulative sum of the starting number of individuals in each compartment 
  // plus/minus the individuals entering that compartment. For the transition described above, we also set
  // arr[ix(t + d + 1, c_next, s_occupancy)]++ /and/ arr[ix(t + d + 1, c_current, s_occupancy)]--, where the second
  // subtraction distinguishes this from the s_transitions slot (which would not function correctly if transitions
  // that needed to be handled were lost)

  // By the end of the simulation loop, arr[ix(..., ..., s_occupancy)] contains the instantaneous occupancy count of
  // each compartment for every time t, and arr[ix(..., ..., s_transitions)] contains the number of transitions
  // into each compartment for every time t. Both of which are useful simulation outputs.

  // We utilise the arr[] array in such a way to maximise memory locality

  for(int t = 0; t < n_steps; t++) {
    for(int c = 0; c < def_n_compartments; c++) {
      
      // Update the current occupancy value to reflect the value at the last timestep
      // counted at arr[ix(t - 1, c, s_occupancy)]
      // plus/minus any transitions (counted by arr[ix(t, c, s_occupancy)])
      if(t > 0)
        arr[ix(t, c, s_occupancy)] = arr[ix(t, c, s_occupancy)] + arr[ix(t - 1, c, s_occupancy)];

      // Also increment occupancy by number of inwards transitions
      arr[ix(t, c, s_occupancy)] = arr[ix(t, c, s_occupancy)] + arr[ix(t, c, s_transitions)];
    }
    
    
    
    // Symptomatic -> ward
    musher::transition_generic_delayed(
      c_symptomatic, c_ward,
      arr[ix(t, c_symptomatic, s_transitions)],
      t, arr, ix,
      symptomatic_ward_delays,
      n_steps
    );

    float pr_ICU = pr_ICU_curve[std::max(t / params.steps_per_day - 5, 0)]; // Adjust for mean 5 day delay onset -> ward admission

    float pr_not_ICU = 1 - pr_ICU;
    float pr_discharge_given_not_ICU = s_data.pr_ward_to_discharge / (1 - s_data.pr_ward_to_ICU);
    float pr_discharge_adj = pr_discharge_given_not_ICU * pr_not_ICU;
    
    // Ward -> ICU, discharged, died
    musher::transition_ward_next(
      t, arr, ix,
      ward_to_discharge_delays, ward_to_ICU_delays, ward_to_death_delays,
      pr_discharge_adj, pr_ICU,
      n_steps, rand);
    
    // ICU -> discharged, died, postICU
    musher::transition_ICU_next(
      t, arr, ix,
      ICU_to_discharge_delays, ICU_to_death_delays, ICU_to_postICU_delays,
      s_data.pr_ICU_to_discharge, s_data.pr_ICU_to_postICU, s_data.pr_postICU_to_death,
      n_steps, rand);
    
    
    // PostICU -> discharged
    musher::transition_generic_delayed(
      c_postICU_to_discharge, c_discharged_postICU,
      arr[ix(t, c_postICU_to_discharge, s_transitions)],
      t, arr, ix,
      postICU_to_discharge_delays,
      n_steps
    );
    
    // PostICU -> died
    musher::transition_generic_delayed(
      c_postICU_to_death, c_died_postICU,
      arr[ix(t, c_postICU_to_death, s_transitions)],
      t, arr, ix,
      postICU_to_death_delays,
      n_steps
    );
  }

  // Pulling together grouped and ungrouped compartment counts for the results
  
  std::vector<int> grouped_occupancy_compartment_labels(params.n_days * def_n_compartment_groups, -1);
  std::vector<int> grouped_occupancy_counts(params.n_days * def_n_compartment_groups, 0);
  std::vector<int> grouped_transitions(params.n_days * def_n_compartment_groups, 0);

  std::vector<int> occupancy_compartment_labels(params.n_days * def_n_compartments, -1);
  std::vector<int> occupancy_counts(params.n_days * def_n_compartments, 0);
  std::vector<int> transitions(params.n_days * def_n_compartments, 0);
  
  for(int c = 0; c < def_n_compartments; c++) {
    
    int c_grp = group_compartment(c);
    for(int day = 0; day < params.n_days; day++) {
      // Sample occupancy counts from the steps_per_dayth slot
      int n_occupancy = arr[ix(day * params.steps_per_day, c, s_occupancy)];

      // But sum transitions across each slot within a day
      int n_transitions = 0;
      for(int i = 0; i < params.steps_per_day; i++) 
        n_transitions += arr[ix(day * params.steps_per_day + i, c, s_transitions)];
      
      
      occupancy_compartment_labels[c * params.n_days + day] = c;
      occupancy_counts[c * params.n_days + day] = n_occupancy;
      transitions[c * params.n_days + day] = n_transitions;
      
      grouped_occupancy_compartment_labels[c_grp * params.n_days + day] = c_grp;
      grouped_occupancy_counts[c_grp * params.n_days + day] += n_occupancy;
      grouped_transitions[c_grp * params.n_days + day] += n_transitions;
    }
  }
  
  mush_results results;
  
  results.grouped_occupancy_compartment_labels = grouped_occupancy_compartment_labels;
  results.grouped_occupancy_counts = grouped_occupancy_counts;
  results.grouped_transitions = grouped_transitions;

  results.occupancy_compartment_labels = occupancy_compartment_labels;
  results.occupancy_counts = occupancy_counts;
  results.transitions = transitions;
  
  return results;
}


template <typename F>
void musher::transition_ward_next(
    int t,
    std::vector<int> &arr,
    F ix,
    
    std::vector<int> &ward_to_discharge_delays,
    std::vector<int> &ward_to_ICU_delays,
    std::vector<int> &ward_to_death_delays,
    float pr_ward_to_discharge,
    float pr_ward_to_ICU,
    
    int n_steps,
    
    std::mt19937 &rand
) {
  int n_to_transition = arr[ix(t, c_ward, s_transitions)];
  
  int n_to_discharge = 0;
  int n_to_ICU = 0;
  int n_to_death = 0;
  
  
  std::uniform_real_distribution<> dist(0, 1);
  for(int i = 0; i < n_to_transition; i++) {
    float rand_sample = dist(rand);
    
    if(rand_sample < pr_ward_to_discharge) {
      n_to_discharge++;
    } else if(rand_sample < pr_ward_to_discharge + pr_ward_to_ICU) {
      n_to_ICU++;
    } else{
      n_to_death++;
    }
  }
  
  
  musher::transition_generic_delayed(
    c_ward, c_discharged_ward, n_to_discharge,
    t, arr, ix,
    ward_to_discharge_delays,
    n_steps
  );
  musher::transition_generic_delayed(
    c_ward, c_ICU, n_to_ICU,
    t, arr, ix,
    ward_to_ICU_delays,
    n_steps
  );
  musher::transition_generic_delayed(
    c_ward, c_died_ward, n_to_death,
    t, arr, ix,
    ward_to_death_delays,
    n_steps
  );
  
}



template <typename F>
void musher::transition_ICU_next(
    int t,
    std::vector<int> &arr,
    F ix,
    
    std::vector<int> &ICU_to_discharge_delays,
    std::vector<int> &ICU_to_death_delays,
    std::vector<int> &ICU_to_postICU_delays,
    float pr_ICU_to_discharge,
    float pr_ICU_to_postICU,
    float pr_postICU_to_death,
    
    int n_steps,
    
    std::mt19937 &rand
) {
  int n_to_transition = arr[ix(t, c_ICU, s_transitions)];
  
  int n_to_discharge = 0;
  int n_to_death = 0;
  int n_to_postICU_death = 0;
  int n_to_postICU_discharge = 0;
  
  
  std::uniform_real_distribution<> dist(0, 1);
  for(int i = 0; i < n_to_transition; i++) {
    float rand_sample = dist(rand);
    
    if(rand_sample < pr_ICU_to_discharge) {
      n_to_discharge++;
    } else if(rand_sample < pr_ICU_to_discharge + pr_ICU_to_postICU) {
      float postICU_sample = dist(rand);
      
      if(postICU_sample < pr_postICU_to_death)
        n_to_postICU_death++;
      else
        n_to_postICU_discharge++;
      
    } else {
      n_to_death++;
    }
  }
  
  
  musher::transition_generic_delayed(
    c_ICU, c_discharged_ICU, n_to_discharge,
    t, arr, ix,
    ICU_to_discharge_delays,
    n_steps
  );
  musher::transition_generic_delayed(
    c_ICU, c_died_ICU, n_to_death,
    t, arr, ix,
    ICU_to_death_delays,
    n_steps
  );
  musher::transition_generic_delayed(
    c_ICU, c_postICU_to_death, n_to_postICU_death,
    t, arr, ix,
    ICU_to_postICU_delays,
    n_steps
  );
  musher::transition_generic_delayed(
    c_ICU, c_postICU_to_discharge, n_to_postICU_discharge,
    t, arr, ix,
    ICU_to_postICU_delays,
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
    
    // We need to handle this transition in the future
    arr[ix(t_set, c_to, s_transitions)]++;
    
    // Indicate that occupancy in the future for our current compartment decreases
    arr[ix(t_set, c_from, s_occupancy)]--;
  }
}


std::vector<int> musher::make_delay_samples(
    int n_samples, double shape, double scale, int steps_per_day, std::mt19937 &rand
  ) {
  std::vector<int> samples(n_samples);
  
  std::gamma_distribution<> delay_dist(shape, scale);
  
  for(int i = 0; i < n_samples; i++)
    samples[i] = std::floor(delay_dist(rand) * steps_per_day);
  
  return samples;
}


int musher::group_compartment(int compartment_id) {
  switch(compartment_id) {
  case c_symptomatic:
    return cgrp_symptomatic;
  case c_ward:
  case c_postICU_to_death:
  case c_postICU_to_discharge:
    return cgrp_ward;
  case c_ICU:
    return cgrp_ICU;
  case c_discharged_ward:
  case c_discharged_ICU:
  case c_discharged_postICU:
    return cgrp_discharged;
  case c_died_ward:
  case c_died_ICU:
  case c_died_postICU:
    return cgrp_died;
  default:
    return -1;
    
  }
}
