#include "cellular_potts_init_configuration.hpp"
/*==========
  DATA CLASS
===========*/
/*======================
  Methods
  =======================*/
//input methods
void init_configuration_class::init_configuration_class(){
};
//
/*==========
  CONTROL CLASS
  ===========*/
  /*======================
    Methods
   =======================*/
void init_configuration_system::initialize_configuration(
							 const model_parameters_cellular_potts_class & model,
							 const type_system_class & type_system
							 )
{
  initialize_parameter(
		       model,
		       type_system
		       );
};
//
void init_configuration_system::initialize_parameter(
						     const model_parameters_cellular_potts_class & model,
						     const type_system_class & type_system
						     )
{
  //
  number_of_sites=model.get_number_of_sites();
  number_of_cells=model.get_number_of_cells();
  number_of_cell_types=model.get_number_of_cell_types();
  number_of_neighbor_sites=site_system.get_number_of_neighbor_sites()
  space_dimension=model.get_space_dimension();
  buffer_cell=cell_system.get_buffer_cell();
  buffer_type=cell_type_system.get_buffer_type();
  //
};
//
void init_configuration_system::input_initial_configuration_setting(
								    )
{
  int work_int=0;
  work_int=io_method.get_input_int(
				   "site_setting_input" ,
				   type_system_class::structure_item
				   );
  
};
//
void type_system_class::set_structure_type_item(
						int type_index,
						std::string child
						)
{
  io_cellular_potts io_method;
  return io_method.generate_structure(
				      "type",
				      type_index,
				      child
				      );
};
