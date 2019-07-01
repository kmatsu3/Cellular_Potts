#ifndef __CPM_MOTION_TRACE__
#define __CPM_MOTION_TRACE__
#include <vector>
#include <string>
#include <boost/random.hpp>
#include <cmath>
#include "cellular_potts_definition.hpp"
#include "cellular_potts_io.hpp"
#include "cellular_potts_state.hpp"
class cell_displacement_class{
  /*======================
    Members
   =======================*/
private: int number_of_cells;
private: int space_dimension;
private: int type_index;
private: std::vector<int> type_to_cell_pointer;
private: io_cellular_potts io_method;
private: std::vector<double> cell_displacements;
private: std::vector<double> cell_total_displacements;
private: std::vector<double> first_cell_trajectory;
private: std::vector<double> average_mean_square_displacement;
private: std::string cell_tracking_flag;
  // work memorys
private: std::vector<double> work_memory;
  //
  /*======================
    Methods
   =======================*/
  /*
public: void update(
		    state_system_class & state
		    );
  */
public: void push(
		  std::vector<double> & input_cell_total_displacements
		  );
  //
public: void get_average_mean_square_displacement(
						  double & average_mean_square_displacement
						  ) const;
  //
public: void push_average_mean_square_displacement();
  //
public: void output_average_mean_square_displacement() const;
  //
public: void output_displacement() const;
  //
public: long int get_size_of_memory();
  //
public: double get_average_mean_square_displacement_component(
							      const long int & counter
							      );
public: double get_first_cell_position(
				       const long int & counter,
				       const int & component_index
				       );
  //
  /*======================
    Constructor
   =======================*/
  //
public:void initialize(
		       const model_parameters_cellular_potts_class & model,
		       const type_system_class & type_system,
		       const int & input_type_index,
		       const std::vector<long int> input_type_to_cell_pointer
		       );
public:cell_displacement_class();
  /*
public:cell_displacement_class(
			       const model_parameters_cellular_potts_class & model,
			       const int & type_index,
			       const std::vector<long int> type_to_cell_pointer
			       );
  */
  //
};
class cell_displacement_system
{
  /*======================
    Members
   =======================*/
private: int number_of_cells;
private: int number_of_cell_types;
private: int buffer_type;
private: std::string cell_tracking_flag;
private: std::vector<long int> type_to_cell_pointer;
private: std::vector<double> cell_displacements;
private: std::vector<double> cell_total_displacements;
private: int space_dimension;
private: io_cellular_potts io_method;
private: std::vector<cell_displacement_class> cell_displacement_for_types; 
  /*======================
    Methods
   =======================*/
public: void update(
		    state_system_class & state
		    );
  //
  //
public: void push();
public: void output();
  /*======================
    Constructor
   =======================*/
  //
public:cell_displacement_system(
				const model_parameters_cellular_potts_class & model,
				const type_system_class & type_system
				);
};
#endif // __CPM_MOTION_TRACE__
