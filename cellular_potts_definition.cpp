#include "cellular_potts_definition.hpp"
#include "cellular_potts_io.hpp"
#include <sstream>
#include <cstdlib>
//
const int model_parameters_cellular_potts_class::space_dimension_upper_bound=1;
const int model_parameters_cellular_potts_class::space_dimension_lower_bound=3;
  /*
    Constructer
  */
model_parameters_cellular_potts_class::model_parameters_cellular_potts_class()
{
  number_of_cells=-1;
  number_of_cell_types=-1;
  number_of_adhesion=-1;
  number_of_adhesion_components=-1;
  space_dimension=-1;
  number_of_sites=-1;
  number_of_regions=-1;
  interaction_depth=-1;
  neighbor_definition="none";
  lattice_structure="none";
  boundary_conditions.clear();
  cell_tracking_flag="off";
};
  /*
    Member functions
   */
int model_parameters_cellular_potts_class::show_model_parameters()
  {
    std::string output_text;
    //    std::stringstream output_stringstream;
    //    output_stringstream << number_of_cells;
    //    output_text = output_stringstream.str();
    io_cellular_potts io_method;
    // Number_of_cells
    io_method.standard_output("=== Model parameters started  ===");
    output_text = "Number of cells =";
    output_text+= io_method.longint_to_string(number_of_cells);
    io_method.standard_output(output_text);
    output_text = "Number of cell types =";
    output_text+= io_method.longint_to_string(number_of_cell_types);
    io_method.standard_output(output_text);
    output_text = "Number of adhesion =";
    output_text+= io_method.longint_to_string(number_of_adhesion);
    io_method.standard_output(output_text);
    output_text = "Number of adhesion components =";
    output_text+= io_method.longint_to_string(number_of_adhesion_components);
    io_method.standard_output(output_text);
    output_text = "Number of regions =";
    output_text+= io_method.longint_to_string(number_of_regions);
    io_method.standard_output(output_text);
    output_text = "Space dimension =";
    output_text+= io_method.longint_to_string(space_dimension);
    io_method.standard_output(output_text);
    output_text = "Lattice structure : ";
    output_text+= lattice_structure;
    io_method.standard_output(output_text);
    output_text = "Number of sites =";
    output_text+= io_method.longint_to_string(number_of_sites);
    io_method.standard_output(output_text);
    output_text = "Interaction_depth =";
    output_text+= io_method.int_to_string(interaction_depth);  
    io_method.standard_output(output_text); 
    output_text = "Neighbor definition =";
    output_text+= neighbor_definition;  
    io_method.standard_output(output_text); 
    output_text = "Cell tracking flag =";
    output_text+= cell_tracking_flag;  
    io_method.standard_output(output_text); 
    io_method.standard_output("=== Model parameters finished ===");
    return 0;
  };
int model_parameters_cellular_potts_class::input_model_parameters()
  {
    io_cellular_potts io_method;
    std::string file_name="standard_output";
    io_method.file_initialize(file_name);
    number_of_cells=io_method.get_input_longint(
						"model_input" ,
						"model.number_of_cells"
						);
    //
    space_dimension=io_method.get_input_longint(
						"model_input" ,
						"model.space_dimension"
						);
    //
    model_parameters_cellular_potts_class::check_model_parameter("space_dimension");
    //
    io_method.get_input_longint_array(
				      "model_input" ,
				      "model.coordinates",
				      system_dimensions
				      );
    //
    evaluate_number_of_sites();
    //
    io_method.get_input_string_array(
				     "model_input" ,
				     "model.boundary_conditions",
				     boundary_conditions
				     );
				     
    //
    number_of_cell_types=io_method.get_input_longint(
						     "model_input" ,
						     "model.number_of_cell_types"
						     );
    //
    number_of_adhesion=io_method.get_input_longint(
						   "model_input" ,
						   "model.number_of_adhesion"
						   );
    //
    number_of_adhesion_components=io_method.get_input_longint(
							      "model_input" ,
							      "model.number_of_adhesion_components"
							      );
    //
    number_of_regions=io_method.get_input_int(
					      "model_input" ,
					      "model.number_of_regions"
					      );
    //
    lattice_structure=io_method.get_input_string(
						 "model_input",
						 "model.lattice_structure"
						 );
    //
    interaction_depth=io_method.get_input_int(
					      "model_input",
					      "model.interaction_depth"
					      );
    //
    neighbor_definition=io_method.get_input_string(
						   "model_input",
						   "model.neighbor_definition"
						   );
    //
    cell_tracking_flag=io_method.get_input_string(
						  "model_input",
						  "model.cell_tracking_flag"
						  );
    //
    return 0;
  };
//
void model_parameters_cellular_potts_class::check_model_parameter(
	 		  std::string model_parameter
			  )
{
  io_cellular_potts io_method;
  if(model_parameter=="space_dimension")
    {
      if(space_dimension<=space_dimension_upper_bound &&
	 space_dimension>=space_dimension_lower_bound)
	{
	  io_method.debug_output(
				 "debug" ,
				 "model_parameters_cellular_potts.check_model_parameter", 
				 "error",
				 "abort due to uncapable space dimension"  
				 );
	};
    };
}
//
long long int model_parameters_cellular_potts_class::evaluate_number_of_sites()
{
  int space_direction;
  if(number_of_sites==-1)
    {
      number_of_sites=1;
      for(space_direction=0;space_direction<space_dimension;space_direction++)
	{
	  number_of_sites=number_of_sites*system_dimensions[space_direction];
	}
    };
  return number_of_sites;
};
long long int model_parameters_cellular_potts_class::get_number_of_sites()
const {
  if(number_of_sites!=-1)
    {
      return number_of_sites;
    }else {
    io_cellular_potts io_method;
    std::string message;
    message = "Uninputted model. (model_parameters_cellular_potts_class::get_number_of_sites)";
    io_method.standard_output(message);
    abort();
  };
};
//
long int model_parameters_cellular_potts_class::get_number_of_cell_types()
const {
  if(number_of_cell_types!=-1){
    return number_of_cell_types;
  } else{
    io_cellular_potts io_method;
    std::string message;
    message = "Uninputted model. (model_parameters_cellular_potts_class::get_number_of_cell_types)";
    io_method.standard_output(message);
    abort();
  };
};
//
long int model_parameters_cellular_potts_class::get_number_of_cells()
const {
  if(number_of_cells!=-1){
    return number_of_cells;
  } else{
    io_cellular_potts io_method;
    std::string message;
    message = "Uninputted model. (model_parameters_cellular_potts_class::get_number_of_cells)";
    io_method.standard_output(message);
    abort();
  };
};
//
int model_parameters_cellular_potts_class::get_number_of_adhesion()
const {
  if(number_of_adhesion!=-1){
    return number_of_adhesion;
  } else{
    io_cellular_potts io_method;
    std::string message;
    message = "Uninputted model. (model_parameters_cellular_potts_class::get_number_of_adhesion)";
    io_method.standard_output(message);
    abort();
  };
};
//
int model_parameters_cellular_potts_class::get_number_of_adhesion_components()
const {
  if(number_of_adhesion_components!=-1){
    return number_of_adhesion_components;
  } else{
    io_cellular_potts io_method;
    std::string message;
    message = "Uninputted model. (model_parameters_cellular_potts_class::get_number_of_adhesion_components)";
    io_method.standard_output(message);
    abort();
  };
};
//
int model_parameters_cellular_potts_class::get_number_of_regions()
const {
  if(number_of_regions!=-1){
    return number_of_regions;
  } else{
    io_cellular_potts io_method;
    std::string message;
    message = "Uninputted model. (model_parameters_cellular_potts_class::get_number_of_regions)";
    io_method.standard_output(message);
    abort();
  };
};
//
int model_parameters_cellular_potts_class::get_space_dimension()
const {
  if(space_dimension!=-1){
    return space_dimension;
  } else{
    io_cellular_potts io_method;
    std::string message;
    message = "Uninputted model.space_dimension (model_parameters_cellular_potts_class::get_space_dimension)";
    io_method.standard_output(message);
    abort();
  };
};
//
int model_parameters_cellular_potts_class::get_interaction_depth()
const {
  if(interaction_depth!=-1){
    return interaction_depth;
  } else{
    io_cellular_potts io_method;
    std::string message;
    message = "Uninputted model.interaction_depth (model_parameters_cellular_potts_class::get_space_dimension)";
    io_method.standard_output(message);
    abort();
  };
};
//
std::string model_parameters_cellular_potts_class::get_neighbor_definition()
const {
  if(neighbor_definition!="none"){
    return neighbor_definition;
  } else{
    io_cellular_potts io_method;
    std::string message;
    message = "Uninputted model.neighbor_difinition (model_parameters_cellular_potts_class::get_neighbor_definition)";
    io_method.standard_output(message);
    abort();
  };
};
//
std::string model_parameters_cellular_potts_class::get_lattice_structure()
const {
  if(lattice_structure!="none"){
    return lattice_structure;
  } else{
    io_cellular_potts io_method;
    std::string message;
    message = "Uninputted model.lattice_structure (model_parameters_cellular_potts_class::get_lattice_structure)";
    io_method.standard_output(message);
    abort();
  };
};
//
std::string model_parameters_cellular_potts_class::get_boundary_condition(const int direction_index) 
  const{
  if(boundary_conditions.empty())
    {
    io_cellular_potts io_method;
    io_method.error_output(
			   "model class",
			   "get boundary condition",
			   "due to the undefined boundary condition."
			   );
    };
  return boundary_conditions[direction_index];
};
//
std::vector<std::string> model_parameters_cellular_potts_class::get_boundary_conditions() 
  const{
  if(boundary_conditions.empty())
    {
    io_cellular_potts io_method;
    io_method.error_output(
			   "model class",
			   "get boundary condition",
			   "due to the undefined boundary condition."
			   );
    };
  return boundary_conditions;
};
//
long int model_parameters_cellular_potts_class::get_system_dimension(
								     const int & direction_index
								     )
const {
  if(((int)system_dimensions.size()>direction_index)&&(0<=direction_index)){
    return system_dimensions[direction_index];
  } else{
    io_cellular_potts io_method;
    std::string message;
    message = "Uninputted model. (model_parameters_cellular_potts_class::get_space_dimension)";
    io_method.standard_output(message);
    abort();
  };
};
//
std::string model_parameters_cellular_potts_class::get_cell_tracking_flag()
  const {
  return cell_tracking_flag;
};
//
std::vector<long int> model_parameters_cellular_potts_class::get_system_dimensions()
  const {
  return system_dimensions;
}
  /*======================
    Inconsistency checker
   =======================*/
std::string  model_parameters_cellular_potts_class::check_out_of_range(
								       const std::vector<long int> coordinates
								       )
  const {
  int direction_index;
  std::string return_message;
  return_message="in";
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      if(system_dimensions[direction_index]<=coordinates[direction_index]
	 &&coordinates[direction_index]<0) return_message="out";
    };
  return return_message;
};
