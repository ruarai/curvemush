

struct group_data {
  float pr_ward_to_death;
  std::vector<float> pr_ward_to_ICU;
  
  float pr_ICU_to_death;
  float pr_ICU_to_discharge;
  
  float pr_postICU_to_death;
  
  
  float d_shape_symptomatic_to_ward, d_scale_symptomatic_to_ward;
  
  float d_shape_ward_to_discharge, d_scale_ward_to_discharge;
  float d_shape_ward_to_ICU, d_scale_ward_to_ICU;
  float d_shape_ward_to_death, d_scale_ward_to_death;
  
  
  float d_shape_ICU_to_discharge, d_scale_ICU_to_discharge;
  float d_shape_ICU_to_postICU, d_scale_ICU_to_postICU;
  float d_shape_ICU_to_death, d_scale_ICU_to_death;
  
  
  float d_shape_postICU_to_discharge, d_scale_postICU_to_discharge;
  float d_shape_postICU_to_death, d_scale_postICU_to_death;
};