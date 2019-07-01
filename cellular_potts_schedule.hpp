#ifndef __CPM_SCHEDULE__
#define __CPM_SCHEDULE__
#include <vector>
#include <boost/random.hpp>
#include <cmath>
#include <ctime>
#include "cellular_potts_definition.hpp"
#include "cellular_potts_io.hpp"
#include "cellular_potts_cell.hpp"
#include "cellular_potts_type.hpp"
#include "cellular_potts_site.hpp"
#include "cellular_potts_state.hpp"
#include "cellular_potts_region.hpp"
#include "cellular_potts_adhesion.hpp"
#include "cellular_potts_observation.hpp"
#include "cellular_potts_observation_type.hpp"
#include "cellular_potts_motion_trace.hpp"
#include "cellular_potts_simulation.hpp"
class model_parameter_schedule_class {
  /*======================
    Members
    =======================*/
private: std::string class_identifier;
private: long int type_identifier;
private: std::string parameter_identifier;
private: int component_identifier;
private: double start_value;
private: double incriment;
  /*======================
    Methods
   =======================*/
public: void set_class_identifier(const std::string & identifier);
public: std::string get_class_identifier() const;
public: void set_type_identifier(const long int & identifier);
public: long int get_type_identifier() const;
public: void set_component_identifier(const long int & identifier);
public: long int get_component_identifier() const;
public: void set_parameter_identifier(const std::string & identifier);
public: std::string get_parameter_identifier() const;
public: void set_schedule_value(
				const std::string & identifier,
				const double & value
				);
public: double get_schedule_value(
				  const std::string & identifier
				  ) const;
public: void output_schedule() const;
  /*======================
    Constructor
   =======================*/
public: model_parameter_schedule_class();
};
class schedule_system_class {
  /*======================
    Members
   =======================*/
  //parameter schedule class
private: std::vector<model_parameter_schedule_class> parameter_schedule;
  //simulation setting parameter
private: long long int monte_carlo_steps_for_observation;
private: long long int monte_carlo_steps_for_relaxzation;
private: long long int number_of_observations;
private: double beta;
private: double temperature;
  //simulation plotting schedule
private: long int number_of_midplot_time;
private: long int period_of_midplots;
  //
private: long int period_of_observation;
private: int number_of_control_parameters;
private: int number_of_sweep_steps;
  //
  // private class
private: io_cellular_potts io_method;
  // work memory
private:std::string structure_item;
private:std::vector<double> work_type_vector;
private:std::vector<double> work_field_vector;
private:std::vector<double> model_parameters;
private:std::vector<std::string> parameter_titles;
private:int number_of_types;
private:int space_dimension; 
  /*======================
    Methods
   =======================*/
public: void monte_carlo(
			 const model_parameters_cellular_potts_class & model,
			 const cell_system_class & cell_system,
			 type_system_class & cell_type_system,
			 site_system_class & site_system,
			 region_system_class & region_system,
			 simulation_system_class & simulation,
			 state_system_class & state,
			 cell_displacement_system & cell_track,
			 observation_system_class & macro_data,
			 observables_type_system_class & observables,
			 adhesion_system_class & adhesion_system
			 );
  //    
public: void input_parameters(
			      type_system_class & cell_type_system,
			      region_system_class & region_system,
			      adhesion_system_class & adhesion_system,
			      const int & sweep_step
			      );
  //
public: void input_mid_plot_schedule();
  //
public: void set_schedule();
public: void set_structure_type_item(
				     int parameter_index,
				     std::string child
				     );
public: void output_schedule() ;
private: void make_control_parameter_titles(
					     const type_system_class & cell_type_system
					    );
private: inline double incrimentation(const int & parameter_index, const int & sweep_step);
  /*======================
    Constructor
   =======================*/
public: schedule_system_class(
			      const model_parameters_cellular_potts_class & model,
			      const type_system_class & cell_type_system,
			      const simulation_system_class & simulation
			      );
};
#endif // __CPM_SCHEDULE___
