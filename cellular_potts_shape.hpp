#ifndef __CPM_SHAPE__
#define __CPM_SHAPE__
#include <vector>
#include <cmath>
#include "toolbox.hpp"
#include "cellular_potts_definition.hpp"
#include "cellular_potts_io.hpp"
#include "cellular_potts_cell.hpp"
#include "cellular_potts_type.hpp"
#include "cellular_potts_site.hpp"
#include "cellular_potts_state.hpp"
class shape_class{
public: std::vector<double> polar_moments;
  //
public: void initialize(
			long int & number_of_cells
			);
  //
public: shape_class();
  //
};
//
class shape_tensor_class{
public: std::vector<double> shape_moment_tensors;
public: void initialize(
			long int & number_of_cells,
			int & space_dimension
			);
  //
public: shape_tensor_class();
};
//
class shape_observable_for_type_class{
public: std::vector<double> shape_observables;
public: long int start_pointer;
public: long int end_pointer;
  //
public: shape_observable_for_type_class();
};
//
class shape_system_class{
private: std::vector<shape_class> shapes;
private: shape_tensor_class shape_tensors;
private: std::vector<shape_observable_for_type_class> polar_observables;
private: std::vector<shape_observable_for_type_class> polar_observable_variances;
private: std::vector<shape_observable_for_type_class> squared_polar_observables;
private: std::vector<shape_observable_for_type_class> shape_traces;
private: std::vector<shape_observable_for_type_class> shape_trace_variances;
private: std::vector<shape_observable_for_type_class> squared_shape_traces;
private: std::vector<shape_observable_for_type_class> shape_traceless_determinants;
private: std::vector<shape_observable_for_type_class> shape_traceless_determinant_variances;
private: std::vector<shape_observable_for_type_class> squared_shape_traceless_determinants;
private: std::vector<double> tmp_matrix;
private: std::vector<double> tmp_vector;
private: std::vector<double> traces;
private: std::vector<double> traceless_determinants;
private: int number_of_polarities;
private: int number_of_components;
private: long int number_of_cells;
private: long int number_of_sites;
private: int space_dimension;
private: int number_of_cell_types;
private: long int counter;
private: long int buffer_cell;
private: io_cellular_potts io_method;
private: std::string filename;
private: int sweep_step;
private: toolbox tool;
  // work_memory
private: std::vector<long int> work_long_int; // dim=space_dimension
private: std::vector<double> relative_coordinate; // dim=space_dimension
public: bool polar_gen_flag;
  //
public:shape_system_class();
public:void initialize(
		       const model_parameters_cellular_potts_class & model,
		       const cell_system_class & cell_system,
		       const type_system_class & cell_type_system,
		       site_system_class & site_system,
		       int & input_sweep_step
		       );
public:void calculate_shape_tensor(
				   const state_system_class & state,
				   const model_parameters_cellular_potts_class & model,
				   site_system_class & site_system,
				   const cell_system_class & cell_system
				   );
private:void accumlate_tensor(
			      std::vector<double> & characters,
			      std::vector<shape_observable_for_type_class> & observables,
			      std::vector<shape_observable_for_type_class> & observable_variances,
			      std::vector<shape_observable_for_type_class> & squared_observables
			      );
private:void accumlate_shape_tensor();
public:void calculate_polar_correlation(
					const state_system_class & state,
					const cell_system_class & cell_system
					);
private:void calculate_characters(
				  shape_tensor_class & tensors
				  );
private:void accumlate_polar_correlation();
public:void output_polar_correlation(
				     const std::vector<std::string> & parameter_titles,
				     const std::vector<double> & model_parameters
				     );
public:void output_shape_characters(
				     const std::vector<std::string> & parameter_titles,
				     const std::vector<double> & model_parameters
				     );
				    
private:void output_labels(
			   const std::vector<std::string> & parameter_titles,
			   const std::string & filename_head,
			   const int & type_index
			   );
  //
private:inline void calculate_relative_coordinate(
						  const state_system_class & state,
						  const model_parameters_cellular_potts_class & model,
						  site_system_class & site_system,
						  const long long int & site_index,
						  const long int cell_index
						  )
  {    
    site_system.origin_shift(
			     model,
			     state.cell_origins[cell_index],
			     site_index,
			     work_long_int
			     );
    for(int direction_index=0;
	direction_index<space_dimension;
	direction_index++)
      {
	relative_coordinate[direction_index]
	  =(double)work_long_int[direction_index]
	  -(double)state.cell_total_positions[cell_index*space_dimension+direction_index]
	  /(double)state.cell_volumes[cell_index];
      };
  };
};
#endif // __CPM_SHAPE___
