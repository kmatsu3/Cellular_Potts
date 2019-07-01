#ifndef __OBSERVATION_TYPE__
#define __OBSERVATION_TYPE__
#include "cellular_potts_definition.hpp"
#include "cellular_potts_io.hpp"
#include "cellular_potts_cell.hpp"
#include "cellular_potts_type.hpp"
#include "cellular_potts_state.hpp"
#include <vector>
#include <cmath>
class observables_type_class{
public: int type;
public: int space_dimension;
public: long int start_pointer;
public: long int end_pointer;
public: std::vector<long int> number_of_living_cells;
public: std::vector<double> net_polarity;
public: std::vector<double> variance_of_polarity;
public: std::vector<double> net_displacement;
public: std::vector<double> variance_of_displacement;
public: std::vector<double> total_displacement;
public: std::vector<double> represent_position;
public: std::vector<double> averaged_position;
  //
public: observables_type_class();
  //
};
//
class observables_type_system_class{
private: int number_of_cell_types;
private: long int number_of_cells;
private: int space_dimension;
private: int sweep_step;
private: long long int track_steps;
private: long long int track_index;
private: std::vector<observables_type_class> observables_types;
private: std::vector<long long int> cell_volumes;
private: std::vector<double> cell_polarities;
private: std::vector<double> cell_displacements;
private: std::vector<double> total_cell_displacements;
private: std::vector<long int> steps;
  // memory for output
private: std::string filename;
private: std::string message;
private: std::vector<double> averaged_net_polarity;
private: std::vector<double> averaged_variance_of_polarity;
private: std::vector<double> averaged_net_displacement;
private: std::vector<double> averaged_variance_of_displacement;
private: std::vector<double> work_position;
private: double averaged_number_of_cells;
  //
private: long long int iteration_counter;
private: long long int number_of_observations;
private: long long int matrix_dimension;
private: std::vector<double> work_polarity_vector;
private: std::vector<double> work_polarity_variance;
private: std::vector<double> work_displacement_vector;
private: std::vector<double> work_displacement_variance;
private: std::vector<double> work_norm;
private: double work_average_norm;
private: io_cellular_potts io_method;
  //
public: void initialize(
			const long long int & input_number_of_observations,
			const int & input_sweep_step,
			const model_parameters_cellular_potts_class & model,
			const cell_system_class & cell_system,	
			const type_system_class & cell_type_system
			);
  //
private: void load_data(
		       const state_system_class & state
		       );
  //
public: void calculation(
			 const state_system_class & state
			 );
public: void track(
		   const state_system_class & state
		   );
public: void track_push();
public: void output(
		    const std::vector<std::string> & parameter_titles,
		    const std::vector<double> & model_parameters
		    );
private: void output_initialize();
private: void output_labels(
			    const std::vector<std::string> & parameter_titles,
			    const std::string & filename_head,
			    const int & type_index
			    );
private: void calculation_average(const int & type_index);
private: void store_living_cell_number();
private: void value_set(const int & type_index);
private: void calculate_net_value(
				  const int & type_index,
				  const std::vector<double> & input_vectors,
				  std::vector<double> & output_vectors
				  );
private: void calculate_variance(
				 const int & type_index,
				 const std::vector<double> & input_vectors,
				 std::vector<double> & output_vectors
				 );
private: double calculate_norm_average(
				       const int & type_index,
				       const std::vector<double> & input_vectors,
				       const double & diviser
				       );
private: inline double norm(
			    const std::vector<double> & input_vector,
			    const long int & index,
			    const long int & vector_size
			    ) const;
private: inline double vector_sum(
				  const std::vector<double> & input_vector,
				  const long int & start_pointer,
				  const long int & end_pointer
				  ) const;
private: inline void sorted_vector_sum(
				       std::vector<double> & sum_vector,
				       const std::vector<double> & sorted_vectors,
				       const long int & start_pointer,
				       const long int & end_pointer,
				       const long int & number_of_vectors,
				       const int & vector_size
				       ) const;
private: inline void sorted_vector_inner_product(
						 std::vector<double> & sum_vector,
						 const std::vector<double> & sorted_vectors_1,
						 const std::vector<double> & sorted_vectors_2,
						 const long int & start_pointer,
						 const long int & end_pointer,
						 const long int & number_of_vectors,
						 const int & vector_size
						 ) const;
  //
public: observables_type_system_class();
};
#endif // __OBSERVATION_TYPE__
