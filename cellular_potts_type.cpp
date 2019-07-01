#include "cellular_potts_type.hpp"
#include <cstdlib>
  /*======================
    Constructor
   =======================*/
  type_cellular_potts_class::type_cellular_potts_class()
  {
    Number_of_cells=-1;
    Natural_volume=-1;
    Natural_perimeter=-1;
    Species="none";
  };
/*======================
    Method
 =======================*/
//
std::string type_cellular_potts_class::set_type_double_value(
							     const std::string & value_name,
							     const double & value
							     )
{
  std::string return_message="Unexpected error";
  if(value_name=="Balk_modulus")
    {
      Balk_modulus=value;
      return_message = "Successfully done";
    };
  if(value_name=="Persistent_time")
    {
      Persistent_time=value;
      return_message = "Successfully done";
    };
  if(value_name=="Adhesion_sensitivity")
    {
      Adhesion_sensitivity=value;
      return_message = "Successfully done";
    };
  if(value_name=="Polarity_sensitivity")
    {
      Polarity_sensitivity=value;
      return_message = "Successfully done";
    };
  if(value_name=="Quadrapolarity_sensitivity")
    {
      Quadrapolarity_sensitivity=value;
      return_message = "Successfully done";
    };
  if(value_name=="Field_sensitivity")
    {
      Field_sensitivity=value;
      return_message = "Successfully done";
    };
  if(value_name=="Dipolar_Adhesion_Coupling_basal")
    {
      Dipolar_Adhesion_Coupling_basal=value;
      return_message = "Successfully done";
    };
  if(value_name=="Quadrapolar_Adhesion_Coupling_basal")
    {
      Quadrapolar_Adhesion_Coupling_basal=value;
      return_message = "Successfully done";
    };
  if(return_message=="Unexpected error")
    {
      io_cellular_potts io_method;
      std::string message;
      message = "In input, ";
      message+= return_message + value_name;
      message+= " is found(type_cellular_potts_class::set_type_double_value).";
      io_method.standard_output(message);
      std::abort();
    }
  return return_message;
};
//
//
std::string type_cellular_potts_class::set_type_longint_value(
							      const std::string & value_name,
							      const int & value
							      )
{
  std::string return_message="Unexpected error";
  if(value_name=="Number_of_cells")
    {
    Number_of_cells=value;
    return_message = "Successfully done";
    };
  if(value_name=="Natural_volume")
    {
      Natural_volume=value;
      return_message = "Successfully done";
     };
  if(value_name=="Natural_perimeter")
    {
      Natural_perimeter=value;
      return_message = "Successfully done";
    };
  if(return_message=="Unexpected error")
    {
      io_cellular_potts io_method;
      std::string message;
      message = "In input, ";
      message+= return_message + value_name;
      message+= " is found(type_cellular_potts_class::set_type_longint_value).";
      io_method.standard_output(message);
      std::abort();
    }
  return return_message;
};
//
std::string type_cellular_potts_class::set_type_string_value(
							     const std::string & value_name,
							     const std::string & value
							     )
{
  std::string return_message="Unexpected error";
  if(value_name=="Species")
    {
    Species=value;
    return_message = "Successfully done";
    };
  if(return_message=="Unexpected error")
    {
      io_cellular_potts io_method;
      std::string message;
      message = "In input, ";
      message+= return_message;
      message+= " is found(type_cellular_potts_class::set_type_string_value).";
      io_method.standard_output(message);
      std::abort();
    }
  return return_message;
};
//
std::string type_cellular_potts_class::add_type_double_vector_value(
								    const std::string & value_name,
								    const std::vector<double> & value
								    )
{
  std::string return_message="Unexpected error";
  std::vector<double>::const_iterator index=value.begin();
  if(value_name=="Natural_polarity")
    {
      Natural_polarity.clear();
	while(index != value.end())
	  {
	    Natural_polarity.push_back(*index);
	    ++index;
	  };
	return_message = "Successfully inputed";
    };
  if(value_name=="Polar_coupling_table")
    {
      Polar_coupling_table.clear();
	while(index != value.end())
	  {
	    Polar_coupling_table.push_back(*index);
	    ++index;
	  };
	return_message = "Successfully inputed";
    };
  if(value_name=="Isotropic_Adhesion_Coupling_constant")
    {
      Isotropic_Adhesion_Coupling_constant.clear();
	while(index != value.end())
	  {
	    Isotropic_Adhesion_Coupling_constant.push_back(*index);
	    ++index;
	  };
	return_message = "Successfully inputed";
    };
  if(value_name=="Dipolar_Adhesion_Coupling_constant")
    {
      Dipolar_Adhesion_Coupling_constant.clear();
	while(index != value.end())
	  {
	    Dipolar_Adhesion_Coupling_constant.push_back(*index);
	    ++index;
	  };
	return_message = "Successfully inputed";
    };
  if(value_name=="Quadrapolar_Adhesion_Coupling_constant")
    {
      Quadrapolar_Adhesion_Coupling_constant.clear();
	while(index != value.end())
	  {
	    Quadrapolar_Adhesion_Coupling_constant.push_back(*index);
	    ++index;
	  };
	return_message = "Successfully inputed";
    };
  if(return_message=="Unexpected error")
    {
      io_cellular_potts io_method;
      std::string message;
      message = "In input, ";
      message+= return_message;
      message+= " is found (type_cellular_potts_class::add_type_double_vector_value).";
      io_method.standard_output(message);
      std::abort();
    }
      return return_message;
};
//
/*======================
    Functions
 =======================*/
void type_system_class::initialize_type(
					const model_parameters_cellular_potts_class & model
					)
{
  // Temporary variables
  int work_int;
  long int work_longint;
  std::string work_string;
  double work_double;
  std::vector<double> work_vector_double;
  // Temporary classes
  type_cellular_potts_class work_type;
  toolbox tool;
  io_cellular_potts io_method;
  std::string error;
  // 
  if(!cell_type.empty()) cell_type.clear();
  int type_index = 0;
  int number_of_cell_types=model.get_number_of_cell_types();
  while(type_index<number_of_cell_types)
    {
      // Number_of_cells
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "number_of_cells"
						 );
      work_int=io_method.get_input_int(
				       "cell_type_input" ,
				       type_system_class::structure_item
				       );
      work_type.set_type_longint_value(
				       "Number_of_cells",
				       work_int
				       );
      // Species
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "species"
						 );
      work_string=io_method.get_input_string(
					     "cell_type_input" ,
					     type_system_class::structure_item
					     );
      if(work_string=="buffer") buffer_type=type_index;
      work_type.set_type_string_value(
				      "Species",
				      work_string
				      );
      // volume_turm.Natural_volume
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "volumes.natural_volume"
						 );
      work_longint=io_method.get_input_longint(
					       "cell_type_input" ,
					       type_system_class::structure_item
					       );
      work_type.set_type_longint_value(
				       "Natural_volume",
				       work_longint
				       );
      // volume_turm.balk_modulus
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "volumes.balk_modulus"
						 );
      work_double=io_method.get_input_double(
					     "cell_type_input" ,
					     type_system_class::structure_item
					     );
      work_type.set_type_double_value(
				      "Balk_modulus",
				      work_double
				      );
      // Persistent_time
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "persistence.persistent_time"
						 );
      work_double=io_method.get_input_double(
					     "cell_type_input" ,
					     type_system_class::structure_item
					     );
      work_type.set_type_double_value(
				      "Persistent_time",
				      work_double
				      );
      //
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "persistence.adhesion_sensitivity"
						 );
      work_double=io_method.get_input_double(
					     "cell_type_input" ,
					     type_system_class::structure_item
					     );
      work_type.set_type_double_value(
				      "Adhesion_sensitivity",
				      work_double
				      );
      //
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "persistence.polarity_sensitivity"
						 );
      work_double=io_method.get_input_double(
					     "cell_type_input" ,
					     type_system_class::structure_item
					     );
      work_type.set_type_double_value(
				      "Polarity_sensitivity",
				      work_double
				      );
      //
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "persistence.quadrapolarity_sensitivity"
						 );
      work_double=io_method.get_input_double(
					     "cell_type_input" ,
					     type_system_class::structure_item
					     );
      work_type.set_type_double_value(
				      "Quadrapolarity_sensitivity",
				      work_double
				      );
      //
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "external_field.field_sensitivity"
						 );
      work_double=io_method.get_input_double(
					     "cell_type_input" ,
					     type_system_class::structure_item
					     );
      work_type.set_type_double_value(
				      "Field_sensitivity",
				      work_double
				      );
      // Natural_perimeter
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "natural_perimeter"
						 );
      work_int=io_method.get_input_int(
				       "cell_type_input" ,
				       type_system_class::structure_item
				       );
      work_type.set_type_longint_value(
				       "Natural_perimeter",
				       work_int
				       );
      //Natural_polarity.absolute_value
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "natural_polarity.absolute_value"
						 );
      work_double=io_method.get_input_double(
					     "cell_type_input" ,
					     type_system_class::structure_item
					     );
      //Natural_polarity.directions
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "natural_polarity.directions"
						 );
      work_vector_double.clear();
      io_method.get_input_double_array(
				       "cell_type_input",
				       type_system_class::structure_item,
				       work_vector_double
				       );
      // 
      tool.normalize(work_vector_double,work_double);
      // 
      work_type.add_type_double_vector_value(
					     "Natural_polarity",
					     work_vector_double
					     );
	   //
      //Coupling_constant
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "isotropic_adhesion_coupling"
						 );
      work_vector_double.clear();
      io_method.get_input_double_array(
				       "cell_type_input",
				       type_system_class::structure_item,
				       work_vector_double
				       );
      work_type.add_type_double_vector_value(
					     "Isotropic_Adhesion_Coupling_constant",
					     work_vector_double
					     );
      // dipolar
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "dipolar_adhesion_coupling"
						 );
      work_vector_double.clear();
      io_method.get_input_double_array(
				       "cell_type_input",
				       type_system_class::structure_item,
				       work_vector_double
				       );
      work_type.add_type_double_vector_value(
					     "Dipolar_Adhesion_Coupling_constant",
					     work_vector_double
					     );
      //
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "dipolar_adhesion_basal"
						 );
      work_double=io_method.get_input_double(
					     "cell_type_input" ,
					     type_system_class::structure_item
					     );
      work_type.set_type_double_value(
				      "Dipolar_Adhesion_Coupling_basal",
				      work_double
				      );
      // quadrapolar
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "quadrapolar_adhesion_coupling"
						 );
      work_vector_double.clear();
      io_method.get_input_double_array(
				       "cell_type_input",
				       type_system_class::structure_item,
				       work_vector_double
				       );
      work_type.add_type_double_vector_value(
					     "Quadrapolar_Adhesion_Coupling_constant",
					     work_vector_double
					     );
      //
      type_system_class::set_structure_type_item(
						 cell_type.size(),
						 "quadrapolar_adhesion_basal"
						 );
      work_double=io_method.get_input_double(
					     "cell_type_input" ,
					     type_system_class::structure_item
					     );
      work_type.set_type_double_value(
				      "Quadrapolar_Adhesion_Coupling_basal",
				      work_double
				      );
      //
      work_vector_double.clear();
      cell_type.push_back(work_type);
      //
      type_index++;
      //
    };
  //
  std::vector< std::vector<double> > coupling;
  work_vector_double.clear();
  for(type_index=0;type_index<number_of_cell_types;type_index++)
    {
      work_vector_double.push_back(0.0);
    };
  for(type_index=0;type_index<number_of_cell_types;type_index++)
    {
      cell_type[type_index].get_adhesion_couplings(
						   "isotropic",
						   work_vector_double
						   );
      coupling.push_back(work_vector_double);
    };
  if(
     tool.symmetry_check_double(coupling)!="symmetric"
     )
    {
      io_method.error_output(
			     "type_system_class",
			     "initialize_type",
			     "due to irregular istropic adhesion symmetry."
			     );
    };
  coupling.clear();
  for(type_index=0;type_index<number_of_cell_types;type_index++)
    {
      cell_type[type_index].get_adhesion_couplings(
						   "dipolar",
						   work_vector_double
						   );
      coupling.push_back(work_vector_double);
    };
  if(
     tool.symmetry_check_double(coupling)!="symmetric"
     )
    {
      io_method.error_output(
			     "type_system_class",
			     "initialize_type",
			     "due to irregular dipolar adhesion symmetry."
			     );
    };
  coupling.clear();
  for(type_index=0;type_index<number_of_cell_types;type_index++)
    {
      cell_type[type_index].get_adhesion_couplings(
						   "quadrapolar",
						   work_vector_double
						   );
      coupling.push_back(work_vector_double);
    };
  if(
     tool.symmetry_check_double(coupling)!="symmetric"
     )
    {
      io_method.error_output(
			     "type_system_class",
			     "initialize_type",
			     "due to irregular quadrapolar adhesion symmetry."
			     );
    };
  //
  long int total=0;
  for(type_index=0;type_index<number_of_cell_types;type_index++)
    {
      total+=cell_type[type_index].get_number_of_cells();
    }
  if(total>model.get_number_of_cells())
    {
      io_method.error_output(
			     "type_system_class",
			     "initialize_type",
			     "due to inconsistency between total number and types for cells."
			     );
    };
  //
};
//
void type_cellular_potts_class::show()
const {
  io_cellular_potts io_method;
  std::string message;
  int component_index;
  message = "Number_of_cells:";
  message+= io_method.longint_to_string(Number_of_cells);
  io_method.standard_output(message);
  message = "Species:";
  message+= Species;
  io_method.standard_output(message);
  message = "Natural_volume:";
  message+= io_method.longint_to_string(Natural_volume);
  io_method.standard_output(message);
  message = "Natural_perimeter:";
  message+= io_method.longint_to_string(Natural_perimeter);
  io_method.standard_output(message);
  message = "Persistent time:";
  message+= io_method.double_to_string(Persistent_time);
  io_method.standard_output(message);
  message = "Natural_polarity: (";
  for(component_index=0;component_index<(int)Natural_polarity.size();component_index++)
    {
      if(component_index!=0) message += ",";
      message += io_method.double_to_string(Natural_polarity.at(component_index));
    };
  message+=")";
  io_method.standard_output(message);
  message = "Adhesion_sensitivity:";
  message+= io_method.double_to_string(Adhesion_sensitivity);
  io_method.standard_output(message);
  message = "Polarity_sensitivity:";
  message+= io_method.double_to_string(Polarity_sensitivity);
  io_method.standard_output(message);
  message = "Quadrapolarity_sensitivity:";
  message+= io_method.double_to_string(Quadrapolarity_sensitivity);
  io_method.standard_output(message);
  message = "Field_sensitivity:";
  message+= io_method.double_to_string(Field_sensitivity);
  io_method.standard_output(message);
  message = "Isotoropic_Coupling_constant: (";
  for(component_index=0;component_index<(int)Isotropic_Adhesion_Coupling_constant.size();component_index++)
    {
      if(component_index!=0) message += ",";
      message += io_method.double_to_string(Isotropic_Adhesion_Coupling_constant.at(component_index));
    };
  message+=")";
  io_method.standard_output(message);
  message = "Dipolar_Coupling_constant: (";
  for(component_index=0;component_index<(int)Dipolar_Adhesion_Coupling_constant.size();component_index++)
    {
      if(component_index!=0) message += ",";
      message += io_method.double_to_string(Dipolar_Adhesion_Coupling_constant.at(component_index));
    };
  message+=")";
  io_method.standard_output(message);
  message = "Dipolar_Coupling_basal: (";
  message += io_method.double_to_string(Dipolar_Adhesion_Coupling_basal);
  message+=")";
  io_method.standard_output(message);
  message = "Quadrapolar_Coupling_constant: (";
  for(component_index=0;component_index<(int)Quadrapolar_Adhesion_Coupling_constant.size();component_index++)
    {
      if(component_index!=0) message += ",";
      message += io_method.double_to_string(Quadrapolar_Adhesion_Coupling_constant.at(component_index));
    };
  message+=")";
  io_method.standard_output(message);
  message = "Quadrapolar_Coupling_basal: (";
  message += io_method.double_to_string(Quadrapolar_Adhesion_Coupling_basal);
  message+=")";
  io_method.standard_output(message);
};
//
/*
const std::vector<double> type_cellular_potts_class::get_adhesion_couplings(const std::string adhesion_type)
  const {
  std::vector<double> work_vector;
  if(adhesion_type=="isotropic")
    {
      work_vector=Isotropic_Adhesion_Coupling_constant;
    }
  else if(adhesion_type=="dipolar")
    {
      work_vector=Dipolar_Adhesion_Coupling_constant;
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "type_cellular_potts_class",
			     "get_adhesion_couplings",
			     "due to undefined adhesion type"
			     );
    };
  return work_vector;
};
*/
void type_cellular_potts_class::get_adhesion_couplings(
						       const std::string & adhesion_type,
						       std::vector<double> & coupling_constants
						       )						
  const {
  if(adhesion_type=="isotropic")
    {
      if(coupling_constants.size()==Isotropic_Adhesion_Coupling_constant.size())
	{
	  copy(
	       Isotropic_Adhesion_Coupling_constant.begin(),
	       Isotropic_Adhesion_Coupling_constant.end(),
	       coupling_constants.begin()
	       );
	}
      else
	{
	io_cellular_potts io_method;
	io_method.error_output(
			       "type_cellular_potts_class",
			       "get_adhesion_couplings",
			       "due to difference in size for Isotropic_Adhesion_Coupling_constant"
			       );
	};
    }
  else if(adhesion_type=="dipolar")
    {
      if(coupling_constants.size()==Dipolar_Adhesion_Coupling_constant.size())
	{
	  copy(
	       Dipolar_Adhesion_Coupling_constant.begin(),
	       Dipolar_Adhesion_Coupling_constant.end(),
	       coupling_constants.begin()
	       );
	}
      else
	{
	io_cellular_potts io_method;
	io_method.error_output(
			       "type_cellular_potts_class",
			       "get_adhesion_couplings",
			       "due to difference in size for Dipolar_Adhesion_Coupling_constant"
			       );
	};
    }
  else if(adhesion_type=="quadrapolar")
    {
      if(coupling_constants.size()==Quadrapolar_Adhesion_Coupling_constant.size())
	{
	  copy(
	       Quadrapolar_Adhesion_Coupling_constant.begin(),
	       Quadrapolar_Adhesion_Coupling_constant.end(),
	       coupling_constants.begin()
	       );
	}
      else
	{
	io_cellular_potts io_method;
	io_method.error_output(
			       "type_cellular_potts_class",
			       "get_adhesion_couplings",
			       "due to difference in size for Quadrapolar_Adhesion_Coupling_constant"
			       );
	};
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "type_cellular_potts_class",
			     "get_adhesion_couplings",
			     "due to undefined adhesion type"
			     );
    };	
  //
};
//
void type_cellular_potts_class::get_adhesion_basal(
						   const std::string & adhesion_type,
						   double & basal_value
						   ) 
  const {
  if(adhesion_type=="dipolar")
    {
      basal_value=Dipolar_Adhesion_Coupling_basal;
    }
  else if(adhesion_type=="quadrapolar")
    {
      basal_value=Quadrapolar_Adhesion_Coupling_basal;
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "type_cellular_potts_class",
			     "get_adhesion_basal",
			     "due to undefined adhesion type"
			     );
    };	
  //
};
//
void type_system_class::show_typelist()
{
  long int type_index;
  std::string message;
  io_cellular_potts io_method;
  //  
  message = "=== Type List Output ===";
  io_method.standard_output(message);
  for(type_index=0;type_index<(int)cell_type.size();type_index++)
    {
      message = "Cell type:[";
	message+= io_method.longint_to_string(type_index);
	message+= "]";
      io_method.standard_output(message);
      cell_type[type_index].show();
    };
};
//
//
void type_system_class::show_type(const int & type_index)
{
  std::string message;
  io_cellular_potts io_method;
  //  
  message = "=== Type List Output ===";
  io_method.standard_output(message);
  message = "Cell type:[";
  message+= io_method.longint_to_string(type_index);
  message+= "]";
  io_method.standard_output(message);
  cell_type[type_index].show();
};
//
void type_system_class::set_structure_type_item(
						int type_index,
						std::string child
						)
{
  io_cellular_potts io_method;
  structure_item=io_method.generate_structure(
					      "type",
					      type_index,
					      child
					      );
};
//
long int type_system_class::get_number_of_cells(
						const int & type_index
						)
  const {
  if((int)cell_type.size()>type_index&&type_index>=0)
    {
      return cell_type[type_index].get_number_of_cells();
    }
  else
    {
    return long_int_error_number;
  };
};
//
long int type_cellular_potts_class::get_number_of_cells()
  const {
  return Number_of_cells;
};
//
/*
const std::vector<double> type_system_class::get_adhesion_couplings(
								    const int type_index, 
								    const std::string adhesion_type
								    )
  const {
  // debug
  int neighbor_index;
  std::vector<double> data;
  for(neighbor_index=0;neighbor_index<(int)data.size();neighbor_index++)
    {
      data=cell_type[type_index].get_adhesion_couplings(adhesion_type);
      fprintf(stderr,"data:%d,%s,%f\n",type_index,adhesion_type.c_str(),data[neighbor_index]);
    };
  //debug
  return cell_type[type_index].get_adhesion_couplings(adhesion_type);
};
*/
void type_system_class::get_adhesion_couplings(
					       const int & type_index, 
					       const std::string & adhesion_type,
					       std::vector<double> & adhesion_couplings
					       ) 
const {
  cell_type[type_index].get_adhesion_couplings(
					       adhesion_type,
					       adhesion_couplings
					       );
};
//
void type_system_class::get_adhesion_basal(
					   const int & type_index, 
					   const std::string & adhesion_type,
					   double & basal_value
					   ) 
const {
  cell_type[type_index].get_adhesion_basal(
					   adhesion_type,
					   basal_value
					   );
};
//
long int type_cellular_potts_class::get_natural_volume()
  const {
  return Natural_volume;
};
//
double type_cellular_potts_class::get_balk_modulus()
  const {
  return Balk_modulus;
};
//
double type_cellular_potts_class::get_persistent_time()
  const {
  return Persistent_time;
};
//
double type_cellular_potts_class::get_adhesion_sensitivity()
  const {
  return Adhesion_sensitivity;
};
//
double type_cellular_potts_class::get_polarity_sensitivity()
  const {
  return Polarity_sensitivity;
};
//
double type_cellular_potts_class::get_quadrapolarity_sensitivity()
  const {
  return Quadrapolarity_sensitivity;
};
//
double type_cellular_potts_class::get_field_sensitivity()
  const {
  return Field_sensitivity;
};
//
double type_system_class::calculate_cell_scale()
  const {
  int number_of_types=(int)cell_type.size();
  int type_index;
  int non_buffer_number=0;
  std::vector<double> volumes(number_of_types,0.0);
  for(type_index=0;type_index<number_of_types;type_index++)
    { 
      if(type_index!=buffer_type)
	{
	  non_buffer_number=non_buffer_number+1;
	  volumes[type_index]=(double)cell_type[type_index].get_natural_volume();
	};
    };
  return std::accumulate(volumes.begin(),volumes.end(),0.0)/(double)non_buffer_number;
};
//
long int type_system_class::get_natural_volume(const int & type_index)
  const {
  return cell_type[type_index].get_natural_volume();
};
//
std::vector<long int> type_system_class::get_natural_volumes()
  const {
  std::vector<long int> work_vector;
  int type_index;
  for(type_index=0;type_index<(int)cell_type.size();type_index++)
    {
      work_vector.push_back(cell_type[type_index].get_natural_volume());
    };
  return work_vector;
};
//
double type_system_class::get_balk_modulus(const int & type_index)
  const {
  return cell_type[type_index].get_balk_modulus();
};
//
std::vector<double> type_system_class::get_balk_moduluses()
  const {
  std::vector<double> work_vector;
  int type_index;
  for(type_index=0;type_index<(int)cell_type.size();type_index++)
    {
    work_vector.push_back(cell_type[type_index].get_balk_modulus());
    };
  return work_vector;
};
//
double type_system_class::get_persistent_time(const int & type_index)
  const {
  return cell_type[type_index].get_persistent_time();
};
//
void type_system_class::get_persistent_times(
						   std::vector<double> & persistence_times
						   )
  const {
  int type_index;
  for(type_index=0;type_index<(int)cell_type.size();type_index++)
    {
      persistence_times[type_index]=cell_type[type_index].get_persistent_time();
    };
};
//
void type_system_class::get_adhesion_sensitivities(
						   std::vector<double> & adhesion_sensitivities
						   )
  const {
  int type_index;
  for(type_index=0;type_index<(int)cell_type.size();type_index++)
    {
      adhesion_sensitivities[type_index]=cell_type[type_index].get_adhesion_sensitivity();
    };
};
//
double type_system_class::get_adhesion_sensitivity(const int & type_index)
  const {
  return cell_type[type_index].get_adhesion_sensitivity();
};
//
double type_system_class::get_polarity_sensitivity(const int & type_index)
  const {
  return cell_type[type_index].get_polarity_sensitivity();
};
//
double type_system_class::get_quadrapolarity_sensitivity(const int & type_index)
  const {
  return cell_type[type_index].get_quadrapolarity_sensitivity();
};
//
double type_system_class::get_field_sensitivity(const int & type_index)
  const {
  return cell_type[type_index].get_field_sensitivity();
};
//
int type_system_class::get_buffer_type() 
const {
  return buffer_type;
};
//
void type_system_class::set_parameter_double(
					     const int & type_index,
					     const std::string & parameter_identifier,
					     const double & value
					     )
{
  cell_type[type_index].set_type_double_value(
					      parameter_identifier,
					      value
					      );
};
//
void type_system_class::set_parameter_double_vector(
						    const int & type_index,
						    const std::string & parameter_identifier,
						    const std::vector<double> & values
						    )
{
  cell_type[type_index].add_type_double_vector_value(
						     parameter_identifier,
						     values
						     );
};
//
type_system_class::type_system_class()
{
  buffer_type=-1; // default buffer type
  long_int_error_number=-1;
  int_error_number=-1;
};
