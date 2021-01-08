#ifndef __INIT_CONFIG__
#define __INIT_CONFIG__
#include "cellular_potts_definition.hpp"
#include "cellular_potts_simulation.hpp"
#include "cellular_potts_cell.hpp"
#include "cellular_potts_adhesion.hpp"
#include "cellular_potts_type.hpp"
#include "cellular_potts_site.hpp"
class init_configuration_class{
private: int type_index;
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
  //
  /*======================
    Methods
   =======================*/
  //
private: init_configuration_class();
};
class init_configuration_system{
  /*======================
    Definitions
   =======================*/
private:long int number_of_cells;
private:int  number_of_cell_types;
private:int number_of_neighbor_sites;
private:long long int number_of_sites;
private:long int buffer_cell;
private:int buffer_type; 
private:int space_dimension;
  /*======================
    Arrays
   =======================*/
public:std::vector<init_configuration_class> configuration_setting;
  /*======================
    Arrays
   =======================*/
private: init_configuration_class work_init_configuration;
  /*======================
    Methods
   =======================*/
private: io_cellular_potts io_method;
  //
public:initialize_configuration(
				const model_parameters_cellular_potts_class & model,
				const site_system_class & site_system,
				const cell_system_class & cell_system,
				const type_system_class & type_system
				);
  //
private: initialize_parameter(
			      const model_parameters_cellular_potts_class & model,
			      const type_system_class & type_system
			      );
  //
private: std::string set_structure_type_item(
					     int type_index,
					     std::string child
					     );
};
#endif // __INIT_CONFIG__
