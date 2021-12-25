
#pragma once
#include "mush_params.h"
#include <random>

class musher {
public:
  static std::vector<int> mush_curve(mush_params params);
  
  template <typename F>
  static void transition_symptomatic_ward(
      int t,
      std::vector<int> &arr,
      F ix,
      std::vector<int> &symptomatic_ward_delays,
      int n_steps
    );
  
  template <typename F>
  static void transition_ward_next(
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
    );
  
  
  
  template <typename F>
  static void transition_generic_delayed(
      int c_from, int c_to,
      int n_to_transition,
      
      int t,
      std::vector<int> &arr,
      F ix,
      std::vector<int> &delay_samples,
      int n_steps
    );
  
private:
  static std::vector<int> make_delay_samples(
      int n_samples, double shape, double scale, std::mt19937 &rand
  );
};
