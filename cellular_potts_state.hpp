#ifndef __STATE__
#define __STATE__
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <cmath>
#include <numeric>
#include <limits>
#include <ctime>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include "toolbox.hpp"
#include "cellular_potts_definition.hpp"
#include "cellular_potts_simulation.hpp"
#include "cellular_potts_cell.hpp"
#include "cellular_potts_adhesion.hpp"
#include "cellular_potts_type.hpp"
#include "cellular_potts_site.hpp"
#include "cellular_potts_region.hpp"
#include "cellular_potts_gnuplot.hpp"
class state_system_class{
  /*======================
    Members
   =======================*/
public: std::vector<long int> configuration;
public: std::vector<long long int> cell_volumes;
private: std::vector<long long int> original_cell_volumes;
private: std::vector<long long int> cell_perimaters;
public: std::vector<double> cell_polarities;
private: std::vector<double> adhesion_polarites;
public: std::vector<long long int> cell_origins;
public: std::vector<long long int> cell_total_positions;
private: double cell_scale;
private: int space_dimension;
  /*
    cell center of mass 
    = (double)(cell_total_positions[]/cell_volumes[]+origion_position[])
   */
  /*======================
    Constants
   =======================*/
private: long int buffer_cell;
private: int buffer_type;
  /*=====================
    2 Temoirary states
   =====================*/
  // Configuration
  /* Initialization setting */
private: std::vector<long int> generation_origin;
private: std::vector<long int> generation_array;
private: std::vector<long int> initial_cell_dimensions;
private: std::vector<long int> initial_cell_separations;
private: std::string generation_method_for_configuration;
private: std::string generation_type_for_configuration;
private: std::string default_type_for_configuration;
private: std::vector<long int> load_configuration_origin;
  //Polarity
private: std::string generation_method_for_polarity;
private: std::string generation_type_for_polarity;
private: std::string default_type_for_polarity;
  // energy_difference
private: double present_adhesion_energy;
private: double candidate_adhesion_energy;
private: double present_volume_energy_difference;
private: double candidate_volume_energy_difference;
private: double polarity_driving_difference;
private: double field_driving_difference;
  // Work calculation
private: double polarity_work;
private: double field_work;
private: double dummy_work;
  /* I/O setting*/
private: std::vector<std::string> configuration_output_type;
private: std::vector<std::string> polarity_output_type;
private: long int state_finilization_step;
  /*Local instance*/
private: io_cellular_potts io_method;
private: toolbox tool;
  /* random number */
private: std::vector<double> random_for_trial;
private: std::vector<long long int> random_for_site_choice;
private: std::vector<int> random_for_neighbor_choice;
  /* update_flag */
private:int update_flag;
  /* simulation parameter*/
private:long int number_of_cells;
private:long long int number_of_sites;
private:int number_of_neighbor_sites;
private:int number_of_cell_types;
private:double beta;
  /* work memory*/
private: std::vector<double> cell_displacements;
private: std::vector<double> plot_memory_position;
private: std::vector<double> plot_displacements;
public: std::vector<double> polar_product; // size = number_of_sites
private: std::vector<double> work_relative_coordinate; // size = space_dimension
private: std::vector<double> work_cell_polarity; // size = space_dimension
private: std::vector<long int> work_long_int; // size = space_dimension
private: std::vector<int> cell_types;
private: std::vector< std::vector <double> > isotropic_adhesion_couplings;
private: std::vector< std::vector <double> > dipolar_adhesion_couplings;
private: std::vector<double> dipolar_adhesion_basals;
private: std::vector< std::vector <double> > quadrapolar_adhesion_couplings;
private: std::vector<double> quadrapolar_adhesion_basals;
private: std::vector<double> work_vector_coupling;
private: std::vector<double> external_field;
private: std::vector<int> finite_field_flag;
private: std::vector<double> work_external_field;
  /* adhesion field*/
private: std::vector<double> adhesion_field;
private: std::vector<double> work_adhesion_field;
  /*Cummulative time caluculater*/
private: clock_t cumulative_time;
  /*======================
    error_return_code
   =======================*/
public: const static int error_code_int=-1;
public: const static long int error_code_longint=-1;
public: const static long long int error_code_longlongint=-1;
public: const static int success_return=1;
public: const static int fail_return=0;
private: int false_value;
private: int true_value;
  /*======================
    Methods
   =======================*/
  //input methods
public: void initialize_state( 
			      const model_parameters_cellular_potts_class & model,
			      cell_system_class & cell_system,
			      const type_system_class & cell_type_system,
			      site_system_class & site_system,
			      simulation_system_class & simulation
			      );
public: void advance_time(
			  const model_parameters_cellular_potts_class & model,
			  cell_system_class & cell_system,
			  const type_system_class & cell_type_system,
			  const adhesion_system_class & adhesion_system,
			  site_system_class & site_system,
			  region_system_class & region_system,
			  simulation_system_class & simulation
			  );
private: void input_initial_configuration_setting();
  //configuration initialiation
private: void generate_configuration(
				     const model_parameters_cellular_potts_class & model,
				     const cell_system_class & cell_system,
				     const type_system_class & cell_type_system,
				     site_system_class & site_system,
				     simulation_system_class & simulation
				     );
private: void load_configuration(
				 const model_parameters_cellular_potts_class & model,
				 const site_system_class & site_system
				 );
private: void set_structure_site_item(
				      const long long int & site_index,
				      const std::string & child,
				      std::string & structure_item
				      );
private: void normal_checkerboard(
				  const model_parameters_cellular_potts_class & model,
				  const cell_system_class & cell_system,
				  const type_system_class & cell_type_system,
				  const site_system_class & site_system,
				  simulation_system_class & simulation
				  );
private: void random_checkerboard(
				  const model_parameters_cellular_potts_class & model,
				  const cell_system_class & cell_system,
				  const type_system_class & cell_type_system,
				  const site_system_class & site_system,
				  simulation_system_class & simulation
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
  //polarity initialization
private: void generate_polarity(
				const model_parameters_cellular_potts_class & model,
				const cell_system_class       & cell_system,
				const type_system_class       & cell_type_system,
				site_system_class       & site_system,
				simulation_system_class & simulation
				) ;
  //
private: void load_polarity(
			    const model_parameters_cellular_potts_class & model,
			    const site_system_class & site_system
			    );
  //
private: void set_structure_polarity_item(
					  const long int & cell_index,
					  const std::string & child,
					  std::string & structure_item
					  );
  //
private: void random_polar(
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
  //output methods
public: void plot_mid_state(
			    const model_parameters_cellular_potts_class & model,
			    const cell_system_class & cell_system,
			    const type_system_class & cell_type_system,
			    const site_system_class & site_system,
			    const long int & time_step,
			    const int & sweep_step
			    ) ;
  //
public: void finalize_state(
			    const model_parameters_cellular_potts_class & model,
			    const cell_system_class & cell_system,
			    const type_system_class & cell_type_system,
			    const site_system_class & site_system
			    );
public: void input_configuration_output_setting();
  //
private: void output_configuration(
				   const model_parameters_cellular_potts_class & model,
				   const cell_system_class & cell_system,
				   const site_system_class & site_system,
				   const long int & time_index,
				   const int & sweep_step
				   ) const;
private: void output_polarity(
			      const model_parameters_cellular_potts_class & model,
			      const site_system_class & site_system,
			      const cell_system_class & cell_system,
			      const long int & time_index,
			      const int & sweep_step
			      ) const;
private: void output_displacement(
				  const model_parameters_cellular_potts_class & model,
				  const site_system_class & site_system,
				  const cell_system_class & cell_system,
				  const long int & time_index,
				  const int & sweep_step
				  ) ;
private: void output_cell_position(
				   const std::string & file_header,
				   const long int & time_index,
				   const int & sweep_step
				   );
  // random number generation
private: void generate_random_number(
				     const model_parameters_cellular_potts_class & model,
				     simulation_system_class & simulation_system
				     );
  //
public: void allocate_random_number_memory(
					   const model_parameters_cellular_potts_class & model,
					   simulation_system_class & simulation
					   );
  // for debug
private: std::string state_consistency_checker(
					       const model_parameters_cellular_potts_class & model,
					       const cell_system_class & cell_system,
					       const type_system_class & cell_type_system,
					       const site_system_class & site_system,
					       const std::string & job_flag
					       ) const;
  //	
private: void init_work();
  //
private: void add_work();
  //
private: double get_energy_difference() const;
  //
private: int Metropholis_check_update(
				      const double energy_difference,
				      const double random_number
				      ) const;
  //
private: void update(
		     const model_parameters_cellular_potts_class & model,
		     const site_system_class & site_system,
		     const long long int & candidate_site,
		     const long int & candidate_cell
		     );
  //
private: void configuration_light_read(
				       const model_parameters_cellular_potts_class & model,
				       site_system_class & site_system
				       );
private: void polarity_light_read(
				  const model_parameters_cellular_potts_class & model,
				  site_system_class & site_system
				  );
  //public: const std::vector<long long int> get_cell_volumes() const; 
public: void get_cell_volumes(
			      std::vector<long long int> & volumes
			      ) const; 
  //public: const std::vector<double> get_cell_polarities() const;
public: void get_cell_polarities(
				 std::vector<double> & cell_polarities
				 ) const;
public: void get_sorted_cell_polarities(
					std::vector<double> & cell_polarities
					) const;
public: void get_cell_displacements(
				    std::vector<double> & output_cell_displacements
				    ) const ;
  //
public: void get_sorted_cell_displacements(
					   std::vector<double> & output_sorted_cell_displacementsF
					   ) const ;
  //
public: void set_model_parameter(
				 const double & parameter,
				 const std::string & parameter_identifier
				 );
  //
public: void set_model_parameter_array(
				       const std::vector<double> & parameters,
				       const std::string & parameter_identifier
				       );
  //
private: inline double external_product_sign(
					     const std::vector<double> & vector_1,
					     const std::vector<double> & vector_2
					     );
  //
public: void get_total_adhesion(
				const model_parameters_cellular_potts_class & model,
				const cell_system_class & cell_system,
				site_system_class & site_system,
				double & total_adhesion_energy_for_cell,
				double & total_intercell_contact_for_cell,
				double & total_intercell_head_contact_for_cell,
				double & total_intercell_tail_contact_for_cell,
				double & total_intercell_vartical_contact_for_cell,
				double & total_intercell_lateral_contact_for_cell,
				double & total_intercell_ordered_contact_for_cell,
				double & total_cellmedium_contact_for_cell,
				double & total_square_polar_product
				);
  //
public: void get_adhesion_field(
				const model_parameters_cellular_potts_class & model,
				const cell_system_class & cell_system,
				site_system_class & site_system,
				const adhesion_system_class & adhesion_system
				);
  //
private: inline double norm(
			    const std::vector<double> & input_vector
			    );
  //
private: inline void gen_polar_product(
				       const model_parameters_cellular_potts_class & model,
				       site_system_class & site_system,
				       std::vector<double> & polar_product
				       )
{
  int site_index;
  for(site_index=0;
      site_index<number_of_sites;
      site_index++)
    {
      if(configuration[site_index]!=buffer_cell)
	{
	  calculate_relative_coordinate(
					model,
					site_system,
					site_index,
					work_relative_coordinate
					);
	  //
	  if(iszero_norm(work_relative_coordinate)==false_value)
	    {
	      get_cell_polarity(work_cell_polarity,site_index);
	      //
	      polar_product[site_index]
		=normalized_product(
				    work_relative_coordinate,
				    work_cell_polarity
				    );
	      //
	    }
	  else
	    {
	      polar_product[site_index]=0.0;
	    };
	}
      else
	{
	};
    };
  
};
  //
private: inline void calculate_relative_coordinate(
						   const model_parameters_cellular_potts_class & model,
						   site_system_class & site_system,
						   const long long int & site_index,
						   std::vector<double> & relative_coordinate
						   ) 
{
  int direction_index;
  long int cell_index=configuration[site_index];
  site_system.origin_shift(
			   model,
			   cell_origins[cell_index],
			   site_index,
			   work_long_int
			   );
  for(direction_index=0;
      direction_index<space_dimension;
      direction_index++)
    {
     relative_coordinate[direction_index]
	=(double)work_long_int[direction_index]
	-(double)cell_total_positions[cell_index*space_dimension+direction_index]
	/(double)cell_volumes[cell_index];
    };
};
  //
private: inline void get_cell_polarity(
				       std::vector<double> & cell_polarity,
				       const long long int & site_index
				       ) 
  const {
  int direction_index;
  long int cell_index=configuration[site_index];
  for(direction_index=0;
      direction_index<space_dimension;
      direction_index++)
    {
      cell_polarity[direction_index]
	=cell_polarities[cell_index*space_dimension+direction_index];
    };
};
  //
private: inline int iszero_norm(
				const std::vector<double> & vector_data
				)
  const {
  int return_value=false_value;
  double norm = std::inner_product(vector_data.begin(),vector_data.end(),vector_data.begin(),0.0);
  if(norm<0.00000000001)
    {
      return_value=true_value;
    };
  return return_value;
};
private: inline double normalized_product(
					  const std::vector<double> & vector_1,
					  const std::vector<double> & vector_2
					  )
  const {
  //
  double norm_1     =std::inner_product(vector_1.begin(),vector_1.end(),vector_1.begin(),0.0);
  double norm_2     =std::inner_product(vector_2.begin(),vector_2.end(),vector_2.begin(),0.0);
  double product_12 =std::inner_product(vector_1.begin(),vector_1.end(),vector_2.begin(),0.0);
  return product_12/pow(norm_1*norm_2,0.5);
  //
};
//
private: inline double isotropic_adhesion_energy(
						 const int & cell_type_1,
						 const int & cell_type_2
						 ) 
  const {
  return isotropic_adhesion_couplings[cell_type_1][cell_type_2];
    };
  //
private: inline double dipolar_adhesion_energy(
					       const int & cell_type_1,
					       const double & polar_product_1,
					       const int & cell_type_2,
					       const double & polar_product_2
					       ) 
  const {
  return dipolar_adhesion_couplings[cell_type_1][cell_type_2]
    *polar_product_1*polar_product_2;
};
  //
public:void set_adhesion_coupling(
				  const type_system_class & cell_type_system
				  );
  //
public:void set_external_field(
			       const site_system_class & site_system,
			       region_system_class & region_system
			       );
private:void copy_volume(
						 const std::vector<long long int> & origin_volume,
						 std::vector<long long int> & copy_volume
						);
  //
  /*======================
    Constructor
   =======================*/
public: state_system_class(
			   const model_parameters_cellular_potts_class & model,
			   const site_system_class & site_system,
			   const cell_system_class & cell_system,
			   const type_system_class & type_system
			   );
};
//
  /*======================
    Members
   =======================*/
class local_state_class{
  /* Variables */
private:long long int local_site;
private:long int cell;
private:int cell_type;
private:long int volume;
private:std::vector<double> polarities;
private:std::vector<long long int> total_coordinates;
private:std::vector<double> relative_coordinates;
private:std::vector<double> external_field;
private:double product_polarity;
private:double product_external_field;
private:int field_flag;
  // for neighbors
private:std::vector<long long int> neighbor_sites;
private:std::vector<long int> neighbor_cells;
private:std::vector<long int> neighbor_types; 
private:std::vector< std::vector <long long int> > neighbor_total_coordinates;
private:std::vector< std::vector <double> > neighbor_polarities;
private:std::vector< std::vector <double> > neighbor_relative_coordinates;
  /* Constant members */
private:int space_dimension;
private:int number_of_neighbor_sites;
private:long int number_of_cells;
private:int number_of_regions;
private:int buffer_type;
private:long int buffer_cell;
  // Isotropic constant
private:std::vector< std::vector <double> > isotropic_adhesion_couplings;
  //private:std::vector< std::vector <int> > isotropic_adhesion_coupling_zeros;
  // Anisotropic constant
private:std::vector< std::vector <double> > dipolar_adhesion_couplings;
private:std::vector< double > dipolar_adhesion_basals;
  //private:std::vector< std::vector <int> > dipolar_adhesion_coupling_zeros;
private:std::vector< std::vector <double> > quadrapolar_adhesion_couplings;
private:std::vector< double > quadrapolar_adhesion_basals;
  //private:std::vector< std::vector <int> > quadrapolar_adhesion_coupling_zeros;
  // Polarity_sensitivity
private:std::vector< double > polarity_sensitivity;
private:std::vector< double > quadrapolarity_sensitivity;
private:std::vector< double > field_sensitivity;
  // Volume parameters
private:std::vector <long int> natural_volumes;
private:std::vector <double> balk_moduluses;
  // Geometry 
private:std::vector< std::vector <long int> > neighbor_shift;
  /*Work memory members*/
private:std::vector<double> work_vector_double_for_neighbors;
private:std::vector<long int> point_vector;
private:std::vector<long long int> work_vector;
private:std::vector<double> work_vector_double;
private:std::vector<long int> original_vector;
private:std::vector<double> work_vector_coupling;
  /*Local instance*/
private: toolbox tool;
private: io_cellular_potts io_method;
  /*Return values*/
private: int false_value;
private: int true_value;
  /*Cumulative time caluculater*/
private: clock_t sub_start_time, sub_end_time, sub_total_time;
private: clock_t cumulative_time_local, cumulative_time;
public: void get_cpu_time(
			  clock_t & cumulative_time,
			  const std::string & job_flag
			  );
  /*======================
    Methods
   =======================*/
public:void initialize_site(
			    const long long int & site_index,
			    const std::vector<long int> & configuration,
			    const std::vector<double> & cell_polarity,
			    const std::vector<double> & field,
			    const std::vector<int> & finite_field_flag,
			    const model_parameters_cellular_potts_class & model,
			    const cell_system_class & cell_system,
			    const type_system_class & cell_type_system,
			    const site_system_class & site_system
			    );
  //
public:void initialize_local_state(
				   const long int & cell_index,
				   const std::vector<long int> & configuration,
				   const std::vector<long long int> & cell_volumes,
				   const std::vector<long long int> & original_cell_volumes,
				   const std::vector<double> & cell_polarities,
				   const std::vector<long long int> & cell_total_positions,
				   const std::vector<long long int> & cell_origins,
				   const model_parameters_cellular_potts_class & model,
				   const cell_system_class & cell_system,
				   const type_system_class & cell_type_system,
				   site_system_class & site_system
				   );
  //
public: void debug_output_local_state(
				   const std::string & data_identifier
				) const;
  //
public: void show_local_state(const std::vector<long int> & configuation) const;
  //
public: void show_local_site_list(
				  const model_parameters_cellular_potts_class &model,
				  const site_system_class & site_system,
				  const std::vector<long int> & configuration
				  ) const;
  //
public: void recalculate_cell_center(
				     const model_parameters_cellular_potts_class &model,
				     const site_system_class & site_system,
				     const long int & cell_index,
				     const std::vector<long int> & configuration,
				     std::vector<double> & center_coordinates
				     ) const;
  //
public: double get_local_isotropic_adhesion_energy();
  //
public: double get_local_dipolar_adhesion_energy();
  //
public: double get_local_polarity_driving_energy();
  //
public: double get_local_adhesion_energy(
					 const adhesion_system_class & adhesion_system
					 );
  //
public: double get_local_field_driving_energy();
  //
public: void set_product_polarity();
  //
public: void set_product_external_field();
  //
private: inline int iszero_norm(
				const std::vector<double> & vector_data
				) const;
  //
private: inline double normalized_product(
					  const std::vector<double> & vector_1,
					  const std::vector<double> & vector_2
					  );
  //
private: inline int iszero_norm_in_adhesive_hamiltonian(
							const std::vector<double> & vector_data
							) const;
  //
private: inline double normalized_product_in_adhesive_hamiltonian(
								  const std::vector<double> & vector_1,
								  const std::vector<double> & vector_2
								  );
  //
private: class add_square_double{
public: inline double operator()(const double summation, const double value)
  {
    return summation + value * value;
  };
};
  // Caution! :: this function calculates volume deformation energy difference for cell on local_site.
  // It does not calculate the volume energy for neither the present state nor the candidate one. 
public: double get_volume_energy_difference(const std::vector<long int> & configuration) const; 
  //
  /*======================
    Constructor
   =======================*/
public:local_state_class(
			 const model_parameters_cellular_potts_class & model,
			 const site_system_class & site_system,
			 const cell_system_class & cell_system,
			 const type_system_class & cell_type_system
			 );
};
//
class polarity_motion_class{
  /*======================
    Members
   =======================*/
  /* Variables */
private: std::vector<double> memorized_positions;
  /* Constants */
private: int space_dimension;
private: long int number_of_cells;
private: int number_of_cell_types;
private: std::vector<double> persistent_time;
private: std::vector<double> adhesion_sensitivity;
private: std::vector<int> cell_types;
private: toolbox tool;
  /* Work_memory */
private: std::vector<double> difference_polarities;
private: std::vector<double> work_vector_double;
private: double projection_to_parpendicular_direction;
private: io_cellular_potts io_method;
  // for update_polarities
  // for get_displacements
private: std::vector<double> work_vector_double_1; //dimension=space_dimension
private: std::vector<double> work_vector_double_2; //dimension=space_dimension
  // 
  //private: std::vector<double> rotation_angles;
  //private: std::vector<double> rotation_axes;
private: double angle;
  /*Work class */
  /*======================
    Methods
   =======================*/
public:void initialize_polarities(
				  const model_parameters_cellular_potts_class & model,
				  const cell_system_class & cell_system,
				  const site_system_class & site_system,
				  const std::vector<double> & input_polarities,
				  const std::vector<long long int> & input_cell_total_positions,
				  const std::vector<long long int> & input_cell_volumes,
				  const std::vector<long long int> & input_cell_origins
				  );
  //
public:void update_polarities(
			      const model_parameters_cellular_potts_class & model,
			      const cell_system_class & cell_system,
			      const site_system_class & site_system,
			      const std::vector<long long int> & input_cell_total_positions,
			      const std::vector<long long int> & input_cell_volumes,
			      const std::vector<long long int> & input_cell_origins,
			      const std::vector<double> & adhesion_field,
			      std::vector<double> & polarities
			      );
  //
public:void get_displacements(
			      const model_parameters_cellular_potts_class & model,
			      const cell_system_class & cell_system,
			      const site_system_class & site_system,
			      const std::vector<long long int> & input_cell_total_positions,
			      const std::vector<long long int> & input_cell_volumes,
			      const std::vector<long long int> & input_cell_origins,
			      std::vector<double> & cell_displacements
			      ) ;
  /*
private: const std::vector<double> polarity_motion_class:: rotation_2d(
								       const std::vector<double> vectors,
								       const std::vector<double> difference_vectors,
								       const long int array_size
								       ) const;
  */
  //
private: inline double external_product_sign(
					     const std::vector<double> & vector_1,
					     const std::vector<double> & vector_2
					     );
  //
  /*======================
    Constructor
   =======================*/
public:polarity_motion_class(
			     const model_parameters_cellular_potts_class & model,
			     const cell_system_class & cell_system,
			     const type_system_class & cell_type_syste
			     );
};
#endif // __STATE__
