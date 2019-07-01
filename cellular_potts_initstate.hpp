#ifndef __INIT_STATE__
#define __INIT_STATE__
#include "cellular_potts_definition.hpp"
#include "cellular_potts_io.hpp"
#include "cellular_potts_state.hpp"
#include <vector>
class initial_region_class{
  /*======================
    Members
   =======================*/
  // Configuration
private: std::vector<long int> cell_indicators;
private: std::vector<long int> generation_origin;
private: std::vector<long int> generation_array;
private: std::vector<long int> initial_cell_dimensions;
private: std::vector<long int> initial_cell_separations;
private: std::vector<double> initial_polarity;
private: std::string generation_method_for_configuration;
private: std::string generation_type_for_configuration;
private: std::string default_type_for_configuration;
private: std::string generation_method_for_polarity;
private: std::string generation_type_for_polarity;
private: std::string default_type_for_polarity;
  // Constants
private: int space_dimension;
  /*======================
    Methods
   =======================*/
public: void set_longint_vector_value(
				      const std::string & value_identifier,
				      const std::vector<long int> & input_vector
				      );
  //
public: void set_string_value(
			      const std::string & value_identifier,
			      const std::string & input_value
			      );
  //
public: void set_double_vector_value(
				     const std::string & value_identifier,
				     const std::vector<double> & input_vector
				     );
  //
public: void get_longint_vector_value(
				      const std::string & value_identifier,
				      std::vector<long int> & output_vector
				      );
  //
public: void get_string_value(
			      const std::string & value_identifier,
			      std::string & output_value
			      );
  //
public: void get_double_vector_value(
				     const std::string & value_identifier,
				     std::vector<double> & output_vector
				     );
  /*======================
    Constructor
   =======================*/
public: initial_region_class(
			     const model_parameters_cellular_potts_class & model
			     );
};
//
class initial_state_system_class{
  /*======================
    Members
   =======================*/
private: std::vector<initial_region_class> initial_region;
  // dimensions
private: int space_dimension;
private: long int number_of_cells;
private: long long int number_of_sites;
private: int number_of_regions;
private: long int buffer_cell;
  // temporary values
private: std::string structure_item;
  // temporary classes
private: io_cellular_potts io_method;
private: std::vector<long int> system_dimensions;
private: std::vector<long int> cell_region_end;
private: std::vector<long int> cell_indicators;
private: std::vector<long int> generation_origin;
private: std::vector<long int> generation_array;
private: std::vector<long int> initial_cell_dimensions;
private: std::vector<long int> initial_cell_separations;
  /*======================
    Methods
   =======================*/
private: void initialize_initial_region(
					const model_parameters_cellular_potts_class & model
					);
  //
private: void set_structure_region_item(
					const int & region_index,
					std::string child
					);
  //
private: void generate_configuration(
				     std::vector<long int> & configuration,
				     const model_parameters_cellular_potts_class & model,
				     const cell_system_class & cell_system,
				     const type_system_class & cell_type_system,
				     const site_system_class & site_system,
				     simulation_system_class & simulation
				     );
  //
private: void fill_buffer(
			  std::vector<long int> & configuration
			  );
  //
private: void generate_for_each_region(
				       std::vector<long int> & configuration,
				       const model_parameters_cellular_potts_class & model,
				       const cell_system_class & cell_system,
				       const type_system_class & cell_type_system,
				       const site_system_class & site_system,
				       simulation_system_class & simulation,
				       const int & region_index
				       );
  
  /*-------------------*/
  //
private: void random_checkerboard(
				  std::vector<long int> & configuration,
				  const model_parameters_cellular_potts_class & model,
				  const cell_system_class & cell_system,
				  const type_system_class & cell_type_system,
				  const site_system_class & site_system,
				  simulation_system_class & simulation,
				  const int & region_index
				  );
  //
private: long int site_to_cell(
			       const model_parameters_cellular_potts_class & model,
			       const cell_system_class & cell_system,
			       const type_system_class & cell_type_system,
			       const site_system_class & site_system,
			       const long long int & site_index,
			       const std::string & pattern
			       );
  //
private: long long int cell_coordinate_to_index(
						const model_parameters_cellular_potts_class & model,
						const cell_system_class & cell_system,
						const std::vector<long int> & cell_coordinates
						) const;
  //
private: void random_polar(
			   std::vector<double> & cell_polarities,
			   const model_parameters_cellular_potts_class & model,
			   const cell_system_class       & cell_system,
			   const type_system_class       & cell_type_system,
			   const site_system_class       & site_system,
			   simulation_system_class & simulation
			   );
  //
private: void random_shuffle_type_of_cell(
					  const model_parameters_cellular_potts_class & model,
					  cell_system_class & cell_system,
					  const type_system_class & cell_type_system,
					  simulation_system_class & simulatiom
					  )const;
  //
private: void random_shuffle_for_cells(
				       std::vector<long int> & configuration,
				       simulation_system_class & simulation
				       );
  //
private: void shuffle_table_generator_long_integer(
						   simulation_system_class & simulation,
						   std::vector<long int> & shuffle_table
						   );
  //
private: void shuffle_table_generator_for_type(
					       const model_parameters_cellular_potts_class & model,
					       const type_system_class & cell_type_system,
					       simulation_system_class & simulation,
					       std::vector<int> & shuffle_table
					       ) const;
  //
  /*======================
    Constructor
   =======================*/
public: initial_state_system_class(
				   const model_parameters_cellular_potts_class & model,
				   const cell_system_class & cell_system
				   );
};
#endif //__INIT_STATE__

