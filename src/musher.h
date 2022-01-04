
#pragma once
#include "mush_params.h"
#include "mush_compartments.h"
#include "group_data.h"
#include <random>


struct mush_results {
  std::vector<int> grouped_occupancy_counts;
  std::vector<int> grouped_occupancy_compartment_labels;
};

class musher {
public:
  static mush_results mush_curve(
      mush_params params,
      
      std::vector<int> case_curve,

      group_data g_data
  );
  
  template <typename F>
  static void transition_ward_next(
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
    );
  
  template <typename F>
  static void transition_ICU_next(
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
      int n_samples, double shape, double scale, int steps_per_day, std::mt19937 &rand
  );
  
  static int group_compartment(int compartment_id);
  
};
