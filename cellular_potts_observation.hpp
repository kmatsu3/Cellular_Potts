#ifndef __CPM_OBSERVATION__
#define __CPM_OBSERVATION__
#include <vector>
#include <boost/random.hpp>
#include <cmath>
#include "cellular_potts_definition.hpp"
#include "cellular_potts_io.hpp"
#include "cellular_potts_cell.hpp"
#include "cellular_potts_type.hpp"
#include "cellular_potts_site.hpp"
#include "cellular_potts_state.hpp"
class observables_class{
  /*======================
    Members
   =======================*/
private: int space_dimension;
private: double number_of_living_cells;
private: std::vector<double> net_polarity;
private: std::vector<double> variance_of_polarity;
private: double absolute_value_of_polarity;
private: std::vector<double> net_displacements;
private: std::vector<double> variance_of_displacements;
private: double absolute_value_of_displacement;
private: double cell_squared_displacement;
private: double total_adhesion_energy_for_cell;
private: double total_volume_energy_for_cell;
private: double total_intercell_contact_for_cell;
private: double total_intercell_head_contact_for_cell;
private: double total_intercell_tail_contact_for_cell;
private: double total_intercell_vartical_contact_for_cell;
private: double total_intercell_lateral_contact_for_cell;
private: double total_intercell_ordered_contact_for_cell;
private: double total_cellmedium_contact_for_cell;
private: double total_square_polar_product;
private: double volume_fraction;
private: double beta;
  /*======================
    Methods
   =======================*/
public: void initialize(
			const model_parameters_cellular_potts_class & model
			);
  //
public: observables_class(
			  const model_parameters_cellular_potts_class & model
			  );
  //
public: void set_number_of_living_cells( 
					const double & number_of_living_cells
					);
  //
public: double get_number_of_living_cells() const ;
  //
public: void set_net_polarity( 
			      const std::vector<double> & net_polarity,
			      const std::vector<double> & valiance_of_polarity
			       );
  //
public: void get_net_polarity( 
			      std::vector<double> & net_polarity,
			      std::vector<double> & valiance_of_polarity,
			      double & absolute_value_of_polarity
			       ) const ;
  //
public: void set_net_displacement( 
				  const std::vector<double> & net_displacement,
				  const std::vector<double> & valiance_of_displacement
				   );
  //
public: void get_net_displacement( 
				  std::vector<double> & net_displacement,
				  std::vector<double> & valiance_of_displacement
				   ) const;
  //
public: void set_cell_squared_displacement( 
					   const double & input_average
					   );
  //
public: double get_cell_squared_displacement() const ;
  //
public: void set_volume( 
			const double & input_total_volume_energy,
			const double & input_volume_fraction
			 );
  //
public: double get_component( 
			     const std::string & data_identifier,
			     const int & component_index
			      ) const ;
  //
public: void set_total_adhesion(
				const double & total_adhesion_energy_for_cell,
				const double & total_intercell_contact_for_cell,
				const double & total_intercell_head_contact_for_cell,
				const double & total_intercell_tail_contact_for_cell,
				const double & total_intercell_vartical_contact_for_cell,
				const double & total_intercell_lateral_contact_for_cell,
				const double & total_intercell_ordered_contact_for_cell,
				const double & total_cellmedium_contact_for_cell,
				const double & total_square_polar_product
				); 
  //
};
//
class observation_system_class{
  /*======================
    Members
   =======================*/
private: std::vector<observables_class> observables;
private: std::vector<std::string> observable_settings;
private: toolbox tool;
private: io_cellular_potts io_method;
  // Model parameters
private: int space_dimension;
private: long int number_of_cells;
private: int number_of_cell_types;
private: int number_of_observables;
private: int number_of_control_parameters;
  // Obsrvation parameters
private: long long int steps_for_observation;
private: long long int number_of_observations;
private: long long int period_of_observation;
  // Work Memory
  // Work flags
private: int load_volume_flag;
private: int load_polarity_flag;
private: int load_displacement_flag;
private: int load_value;
private: int unload_value;
  // Flags
private: std::string observation_output_format;
  // Work Vectors
private: std::vector<long long int> cell_volumes;
private: std::vector<long long int> natural_volumes;
private: std::vector<double> balk_modulus;
private: std::vector<double> cell_polarities;
private: std::vector<double> cell_displacements;
private: std::vector<double> work_average;
private: std::vector<double> work_variance;
private: std::vector<double> work_producted_vectors; // number_of_cells*space_dimension
private: long int buffer_cell;
private: int buffer_type;
  // Instance
  /*======================
    Methods
   =======================*/
public: void set_observable_settings();
  //
private:void initialize_load_flags();
  //
public:void observation(
			const model_parameters_cellular_potts_class & model,
			const cell_system_class & cell_system,	
			const type_system_class & cell_type_system,
			site_system_class & site_system,
			state_system_class & state,
			const long long int & monte_carlo_step
			);
  //
public:void output_observation(
			       std::vector<std::string> & parameter_titles,
			       std::vector<double> & model_parameters
			       );
  //
public:void output_average_observation(
				       std::vector<std::string> & parameter_titles,
				       std::vector<double> & model_parameters
				       );
  //
private: void load_volumes(
			   const state_system_class & state
			   );
  //
private: double calculate_number_of_living_cells(
						 const model_parameters_cellular_potts_class & model,
						 const state_system_class & state  
						 ) const;
  //
private: void load_polarities(
			      const state_system_class & state
			      );
  //
private: void calculate_net_polarity(
				     const model_parameters_cellular_potts_class & model,
				     const state_system_class & state,
				     std::vector<double> & net_polarity
				     ) const;
  //
private: inline void sorted_vector_sum(
				       std::vector<double> & sum_vector,
				       const std::vector<double> & vectors,
				       const long int & number_of_vectors,
				       const int & vector_size
				       ) const;
  //
public: void calculate_cell_variance_of_polarity(
						 const model_parameters_cellular_potts_class & model,
						 const state_system_class & state,
						 const std::vector<double> & net_polarity,
						 std::vector<double> & variance_of_polarity
						 ) const;
  //
private: inline void squared_vector(
		 		    const std::vector<double> & input_vector,
				    std::vector<double> & output_vector,
				    const long int & vector_size
				    );
  //
  /*
private: inline void sorted_vector_producted_sum(
						 std::vector<double> & sum_vector,
						 const std::vector<double> & vectors,
						 const long int & number_of_vectors,
						 const int & vector_size
						 ) const;
  */
  //
public: void calculate_net_displacement(
					const model_parameters_cellular_potts_class & model,
					const state_system_class & state,
					std::vector<double> & net_displacement
					) const;
  //
public: void calculate_cell_variance_of_displacement(
						     const model_parameters_cellular_potts_class & model,
						     const state_system_class & state,
						     const std::vector<double> & net_displacement,
						     std::vector<double> & variance_of_displacement
						     ) const;
public: void calculate_cell_average_squared_displacement(
							 const model_parameters_cellular_potts_class & model,
							 const state_system_class & state,
							 double & average_squared_displacement
							 ) const;
  //
private: void load_cell_displacements(
				      const state_system_class & state
				      );
  //
public: double calculate_adhesion_hamiltonian(
					      const cell_system_class & cell_system,
					      const type_system_class & cell_type_system,
					      const state_system_class & state  
					      ) const;
public: double calculate_volume_hamiltonian(
					    const cell_system_class & cell_system,
					    const type_system_class & cell_type_system,
					    const state_system_class & state  
					    ) const;
public: double calculate_contact_number(
					const std::string anisotropy
					) const;
  //
public: void calculate_volume(
			      const type_system_class & cell_type_system,
			      const cell_system_class & cell_system,
			      double & volume_energy,
			      double & volume_fraction
			      );
  //
public: void get_sorted_observables(
				    const long long int & observation_index,
				    std::vector<double> & sorted_observables,
				    std::vector<std::string> & observable_titles,
				    const std::string & job_flag
				    ) const;
  //
  /*======================
    Constructor
   =======================*/
public:observation_system_class(
				const model_parameters_cellular_potts_class & model,
				const cell_system_class & cell_system,
				const type_system_class & cell_type_system,
				const simulation_system_class & simulation
				);
};
#endif // __CPM_OBSERVATION___
