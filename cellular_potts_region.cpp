# include "cellular_potts_region.hpp"
/*======================
  Methods
  =======================*/
void region_cellular_potts_class::set_region_identifier(
							const std::string & input_region_identifier
							)
{
  region_identifier=input_region_identifier;
};
//
std::string region_cellular_potts_class::get_region_identifier()
  const {
  return region_identifier;
};
//
void region_cellular_potts_class::set_boundary(
					       const std::string & boundary_identifier,
					       const std::vector<long int> & input_boundary
					       )
{
  std::vector<long int>::const_iterator direction_index=input_boundary.begin();
  if(boundary_identifier=="upper")
    {
      upper_boundary.clear();
      while(direction_index!=input_boundary.end())
	{
	  upper_boundary.push_back(*direction_index);
	  direction_index++;
	};
    }
  else if(boundary_identifier=="lower")
    {
      lower_boundary.clear();
      while(direction_index!=input_boundary.end())
	{
	  lower_boundary.push_back(*direction_index);
	  direction_index++;
	};
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "region_cellular_potts_class",
			     "set_boundary",
			     "irrigal boundary identidier"
			     );
    };
};
//
long int region_cellular_potts_class::get_boundary(
						   const std::string & boundary_identifier,
						   const int & direction_index
						   )
  const {
  long int return_value=-1;
  if(boundary_identifier=="upper")
    {
      return_value=upper_boundary[direction_index];
    }
  else if(boundary_identifier=="lower")
    {
      return_value=lower_boundary[direction_index];
    };
  return return_value;
};
//
void region_cellular_potts_class::set_external_field(
						     const std::vector<double> & input_field
						     )
{
  external_field.clear();
  std::vector<double>::const_iterator direction_index=input_field.begin();
  while(direction_index!=input_field.end())
    {
      external_field.push_back(*direction_index);
      direction_index++;
    }
};
//
double region_cellular_potts_class::get_field(
					      int & direction_index
					      )
  const {
  return external_field[direction_index];
};
//
void region_cellular_potts_class::show()
  const {
  io_cellular_potts io_method;
  std::string message;
  message = "region_identifier = " + region_identifier;
  io_method.standard_output(message);
  message = "Lower Boundary = (";
  for(int direction_index=0;direction_index<(int)lower_boundary.size();direction_index++)
    {
      message += io_method.longint_to_string(lower_boundary[direction_index]);
      if(direction_index<(int)lower_boundary.size()-1) message += ",";
    };
  message+= ")" ;
  io_method.standard_output(message);
  message = "Upper Boundary = ("; 
  for(int direction_index=0;direction_index<(int)upper_boundary.size();direction_index++)
    {
      message += io_method.longint_to_string(upper_boundary[direction_index]);
      if(direction_index<(int)upper_boundary.size()-1) message += ",";
    };
  message+= ")" ;
  io_method.standard_output(message);
  message = "Field = (";
  for(int direction_index=0;direction_index<(int)external_field.size();direction_index++)
    {
      message += io_method.double_to_string(external_field[direction_index]);
      if(direction_index<(int)external_field.size()-1) message += ",";
    };
  message+= ")" ;
  io_method.standard_output(message);
};
  /*======================
    Constructor
   =======================*/
  region_cellular_potts_class::region_cellular_potts_class()
  {
  };  
//
void region_system_class::initialize_region()
{
  int region_index;
  std::string work_string;
  std::vector<long int> work_vector_longint(space_dimension,-1);
  std::vector<double> work_vector_double(space_dimension,0.0);
  region_cellular_potts_class work_region;
  for(region_index=0;region_index<number_of_regions;region_index++)
    {
      // boundary identififer
      set_structure_region_item(
				region_index,
				"region_identifier"
				);
      //
      work_string=io_method.get_input_string(
					     "region_input" ,
					     structure_item
					     );
      //
      work_region.set_region_identifier(work_string);
      // upper boundary
      set_structure_region_item(
				region_index,
				"boundary.upper"
				);
      //
      work_vector_longint.clear();
      //
      io_method.get_input_longint_array(
					"region_input",
					structure_item,
					work_vector_longint
					);
      //
      work_region.set_boundary(
			       "upper",
			       work_vector_longint
			       );
      // lower boundary
      set_structure_region_item(
				region_index,
				"boundary.lower"
				);
      //
      work_vector_longint.clear();
      //
      io_method.get_input_longint_array(
					"region_input",
					structure_item,
					work_vector_longint
					);
      //
      work_region.set_boundary(
			       "lower",
			       work_vector_longint
			       );
      // field
      set_structure_region_item(
				region_index,
				"external_field"
				);
      //
      work_vector_double.clear();
      //
      io_method.get_input_double_array(
				       "region_input",
				       structure_item,
				       work_vector_double
				       );
      //
      work_region.set_external_field(
				     work_vector_double
				     );
      //
      regions.push_back(work_region);
    };
};
//
void region_system_class::set_structure_region_item(
						    const int & type_index,
						    const std::string & child
						    )
{
  io_cellular_potts io_method;
  structure_item=io_method.generate_structure(
					      "region.region_type",
					      type_index,
					      child
					      );
};
//
void region_system_class::get_field(
				    const site_system_class & site_system,
				    const long long int & site_index,
				    std::vector<double> & field
				    )
{
  for(int direction_index=0;direction_index<space_dimension;direction_index++)
    {
      field[direction_index]=0.0;
    };
  for(int region_index=0;region_index<number_of_regions;region_index++)
    {
      if(region_identify(region_index,site_system,site_index)==in_region_flag)
	{
	  for(int direction_index=0;direction_index<space_dimension;direction_index++)
	    {
	      field[direction_index]+=regions[region_index].get_field(direction_index);
	    };
	}
    };
};
//
int region_system_class::region_identify(
					 const int region_index,
					 const site_system_class & site_system,
					 const long long int & site_index
					 )
{
  site_system.get_site_coordinates(
				   site_index,
				   work_vector_for_coordinate
				   );
  int direction_index;
  int return_value=error_region_flag;
  int check=0;
  if(regions[region_index].get_region_identifier()=="in")
    {
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  if(
	     (
	      regions[region_index].get_boundary("upper",direction_index)
	      >=
	      work_vector_for_coordinate[direction_index]
	      )
	     &&
	     (
	      regions[region_index].get_boundary("lower",direction_index)
	      <=
	      work_vector_for_coordinate[direction_index]
	      )
	     ) check++; 
	};
    }
  else if(regions[region_index].get_region_identifier()=="out")
    {
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  if(
	     (
	      regions[region_index].get_boundary("upper",direction_index)
	      <=
	      work_vector_for_coordinate[direction_index]
	      )
	     ||
	     (
	      regions[region_index].get_boundary("lower",direction_index)
	      >=
	      work_vector_for_coordinate[direction_index]
	      )
	     ) check++; 
	};
    };
  if((check==space_dimension&&regions[region_index].get_region_identifier()=="in")
     ||
     (check>=1&&regions[region_index].get_region_identifier()=="out")) 
    {
      return_value=in_region_flag;
    }
  else
    {
      return_value=out_region_flag;
    };
  //
  return return_value;
};
//
void region_system_class::show_regions()
{
  std::string message;
  for(int region_index=0;region_index<number_of_regions;region_index++)
    {
      message="=== Region[" + io_method.int_to_string(region_index) + "]===";
      io_method.standard_output(message);
      regions[region_index].show();
    }
};
//
void region_system_class::show_region(
				      const int & region_index
				      )
{
  std::string message;
  message="=== Region[" + io_method.int_to_string(region_index) + "]===";
  io_method.standard_output(message);
  regions[region_index].show();
};
//
void region_system_class::get_region_field(
					   const int & region_index,
					   std::vector<double> & field
					   )
const {
  for(int direction_index=0;direction_index<space_dimension;direction_index++)
    {
      field[direction_index]=regions[region_index].get_field(direction_index);
    };
};
//
void region_system_class::set_region_field(
					   const int & region_index,
					   const std::vector<double> & field
					   )
{
  regions[region_index].set_external_field(field);
};
//
void region_system_class::test_region(
				      const site_system_class & site_system,
				      const long long int & site_index
				      )
{
  std::vector<double> work_field(space_dimension,0.0);
  site_system.get_site_coordinates(
				   site_index,
				   work_vector_for_coordinate
				   );
  get_field(
	    site_system,
	    site_index,
	    work_field
	    );
  std::string message;
  message = "debug: site[";
  message+= io_method.longlongint_to_string(site_index);
  message+= "],Cod(";
  for(int index=0;index<space_dimension;index++)
    {
      message+=io_method.longint_to_string(work_vector_for_coordinate[index]);
      if(index<space_dimension-1) message+=",";
    };
  message+= "),Field(";
  for(int index=0;index<space_dimension;index++)
    {
      message+=io_method.double_to_string(work_field[index]);
      if(index<space_dimension-1) message+=",";
    };
  message+= ")";
  if((int)(work_field[0])!=0||(int)(work_field[1])!=0) fprintf(stderr,"%s\n",message.c_str());
};
//
region_system_class::region_system_class(
					 const model_parameters_cellular_potts_class & model
					 )
{
  number_of_regions=model.get_number_of_regions();
  space_dimension=model.get_space_dimension();
  int direction_index;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      work_vector_for_coordinate.push_back(0);
    };
  in_region_flag=1;
  out_region_flag=0;
  error_region_flag=-1;
  initialize_region();
  show_regions();
};
