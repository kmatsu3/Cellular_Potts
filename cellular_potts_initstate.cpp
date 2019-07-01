#include "cellular_potts_initstate.hpp"
void initial_region_class::set_longint_vector_value(
						    const std::string & value_identifier,
						    const std::vector<long int> & input_vector
						    )
{
  if(value_identifier=="cell_indicators")
    {
      for(int pointer_index=0;pointer_index<2;pointer_index++)
	{
	  cell_indicators[pointer_index]=input_vector[pointer_index];
	};
    }
  else if(value_identifier=="generation_origin")
    {
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  generation_origin[direction_index]=input_vector[direction_index];
	};
    }
  else if(value_identifier=="generation_array")
    {
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  generation_array[direction_index]=input_vector[direction_index];
	};
    }
  else if(value_identifier=="initial_cell_dimensions")
    {
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	 initial_cell_dimensions[direction_index]=input_vector[direction_index];
	};
    }
  else if(value_identifier=="initial_cell_separations")
    {
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	 initial_cell_separations[direction_index]=input_vector[direction_index];
	};
    };
};
//
void initial_region_class::get_longint_vector_value(
						    const std::string & value_identifier,
						    std::vector<long int> & output_vector
						    )
{
  if(value_identifier=="cell_indicators")
    {
      for(int pointer_index=0;pointer_index<2;pointer_index++)
	{
	  output_vector[pointer_index]=cell_indicators[pointer_index];
	};
    }
  else if(value_identifier=="generation_origin")
    {
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  output_vector[direction_index]=generation_origin[direction_index];
	};
    }
  else if(value_identifier=="generation_array")
    {
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  output_vector[direction_index]=generation_array[direction_index];
	};
    }
  else if(value_identifier=="initial_cell_dimensions")
    {
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  output_vector[direction_index]=initial_cell_dimensions[direction_index];
	};
    }
  else if(value_identifier=="initial_cell_separations")
    {
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  output_vector[direction_index]=initial_cell_separations[direction_index];
	};
    };
};
//
void initial_region_class::set_string_value(
					    const std::string & value_identifier,
					    const std::string & input_value
					    )
{
  if(value_identifier=="generation_method_for_configuration")
    {
      generation_method_for_configuration=input_value;
    }
  else if(value_identifier=="generation_type_for_configuration")
    {
      generation_type_for_configuration=input_value;
    }
  else if(value_identifier=="generation_method_for_polarity")
    {
      generation_method_for_configuration=input_value;
    }
  else if(value_identifier=="generation_type_for_polarity")
    {
      generation_type_for_configuration=input_value;
    }
};
//
void initial_region_class::set_double_vector_value(
						   const std::string & value_identifier,
						   const std::vector<double> & input_vector
						   )
{
  if(value_identifier=="initial_polarity")
    {
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	 initial_polarity[direction_index]=input_vector[direction_index];
	};
    };
};
//
void initial_region_class::get_string_value(
					    const std::string & value_identifier,
					    std::string & output_value
					    )
{
  if(value_identifier=="generation_method_for_configuration")
    {
      output_value=generation_method_for_configuration;
    }
  else if(value_identifier=="generation_type_for_configuration")
    {
      output_value=generation_type_for_configuration;
    }
  else if(value_identifier=="generation_method_for_polarity")
    {
      output_value=generation_method_for_configuration;
    }
  else if(value_identifier=="generation_type_for_polarity")
    {
      output_value=generation_type_for_configuration;
    }
};
//
void initial_region_class::get_double_vector_value(
						   const std::string & value_identifier,
						   std::vector<double> & output_vector
						   )
{
  if(value_identifier=="initial_polarity")
    {
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  output_vector[direction_index]=initial_polarity[direction_index];
	};
    };
};
//
initial_region_class::initial_region_class(
					   const model_parameters_cellular_potts_class & model
					   )
{
  space_dimension=model.get_space_dimension();
  cell_indicators.push_back(-1);
  cell_indicators.push_back(-1);
  for(int direction_index=0;direction_index<space_dimension;direction_index++)
    {
      generation_origin.push_back(-1);
      generation_array.push_back(-1);
      initial_cell_dimensions.push_back(-1);
      initial_cell_separations.push_back(-1);
      initial_polarity.push_back(1.0/sqrt(float(space_dimension)));
    };
  default_type_for_configuration="random";
  default_type_for_polarity="random";
};
  /*======================
    Methods for control system
   =======================*/
void initial_state_system_class::initialize_initial_region(
							   const model_parameters_cellular_potts_class & model
							   )
{
  initial_region.clear();
  initial_region_class work_region(model);
  int work_int;
  std::string work_string;
  std::vector<long int> work_longint_vector(2,-1);
  std::vector<long int> work_longint_space_vector(space_dimension,-1);
  structure_item="configuration_setting.number_of_regions";
  work_int=io_method.get_input_int(
				   "site_setting_input" ,
				   structure_item
				   );
  //
  for(int region_index=0;region_index<number_of_regions;region_index++)
    {
      // cell indicator
      set_structure_region_item(
				region_index,
				"cell_indicator.initial"
				);
      //
      work_longint_vector[0]=io_method.get_input_longint(
							 "site_setting_input",
							 structure_item
							 );
      //
      set_structure_region_item(
				region_index,
				"cell_indicator.final"
				);
      //
      work_longint_vector[1]=io_method.get_input_longint(
							 "site_setting_input",
							 structure_item
							 );
      //
      work_region.set_longint_vector_value(
					   "cell_indicators",
					   work_longint_vector
					   );
      // generation_origin
      set_structure_region_item(
				region_index,
				"origin"
				);
      //
      io_method.get_input_longint_array(
					"site_setting_input" ,
					structure_item,
					 work_longint_space_vector
					);
      //
      work_region.set_longint_vector_value(
					   "generation_origin",
					   work_longint_space_vector
					   );
      // generation_array
      set_structure_region_item(
				region_index,
				"array"
				);
      //
      io_method.get_input_longint_array(
					"site_setting_input" ,
					structure_item,
					 work_longint_space_vector
					);
      //
      work_region.set_longint_vector_value(
					   "generation_origin",
					   work_longint_space_vector
					   );
      // generation_array
      work_longint_space_vector.clear();
      //
      set_structure_region_item(
				region_index,
				"array"
				);
      //
      io_method.get_input_longint_array(
					"site_setting_input" ,
					structure_item,
					 work_longint_space_vector
					);
      //
      work_region.set_longint_vector_value(
					   "generation_origin",
					   work_longint_space_vector
					   );
      // cell dimensions
      work_longint_space_vector.clear();
      //
      set_structure_region_item(
				region_index,
				"cell_dimensions"
				);
      //
      io_method.get_input_longint_array(
					"site_setting_input" ,
					structure_item,
					 work_longint_space_vector
					);
      //
      work_region.set_longint_vector_value(
					   "cell_dimensions",
					   work_longint_space_vector
					   );
      // cell separations
      work_longint_space_vector.clear();
      //
      set_structure_region_item(
				region_index,
				"cell_separations"
				);
      //
      io_method.get_input_longint_array(
					"site_setting_input" ,
					structure_item,
					 work_longint_space_vector
					);
      //
      work_region.set_longint_vector_value(
					   "cell_separations",
					   work_longint_space_vector
					   );
      // generation method for configuration
      set_structure_region_item(
				region_index,
				"configuration.generation_method"
				);
      //      
      work_string=io_method.get_input_string(
					     "site_setting_input" ,
					     structure_item
					     );
      //
      work_region.set_string_value(
				   "generation_method_for_configuration",
				   work_string
				   );
      // generation type for configuration
      set_structure_region_item(
				region_index,
				"configuration.generation_type"
				);
      //      
      work_string=io_method.get_input_string(
					     "site_setting_input" ,
					     structure_item
					     );
      //
      work_region.set_string_value(
				   "generation_type_for_configuration",
				   work_string
				   );
      // generation method for configuration
      set_structure_region_item(
				region_index,
				"polarity.generation_method"
				);
      //      
      work_string=io_method.get_input_string(
					     "site_setting_input" ,
					     structure_item
					     );
      //
      work_region.set_string_value(
				   "generation_method_for_polarity",
				   work_string
				   );
      // generation type for configuration
      set_structure_region_item(
				region_index,
				"polarity.generation_type"
				);
      //      
      work_string=io_method.get_input_string(
					     "site_setting_input" ,
					     structure_item
					     );
      //
      work_region.set_string_value(
				   "generation_type_for_polarity",
				   work_string
				   );
      //
    };
};
//
void initial_state_system_class::random_checkerboard(
						     std::vector<long int> & configuration,
						     const model_parameters_cellular_potts_class & model,
						     const cell_system_class & cell_system,
						     const type_system_class & cell_type_system,
						     const site_system_class & site_system,
						     simulation_system_class & simulation,
						     const int & region_index
						     )
{
  //
  cell_region_end.clear();
  //
  initial_region[region_index].get_longint_vector_value(
							"cell_indicators",
							cell_indicators
							);
  initial_region[region_index].get_longint_vector_value(
							"generation_origin",
							generation_origin
							);
  initial_region[region_index].get_longint_vector_value(
							"generation_array",
							generation_array
							);
  initial_region[region_index].get_longint_vector_value(
							"initial_cell_dimensions",
							initial_cell_dimensions
							);
  initial_region[region_index].get_longint_vector_value(
							"initial_cell_separations",
							initial_cell_separations
							);
  //
  for(int direction_index=0;direction_index<space_dimension;direction_index++)
    {
      cell_region_end.push_back(
				generation_origin[direction_index]
				+generation_array[direction_index]
				*(initial_cell_dimensions[direction_index]
				  +initial_cell_separations[direction_index])
				);
    };
  //
  if(model.check_out_of_range(cell_region_end)=="out")
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "state class",
			     "random_checkerboard",
			     "due to out of range of initial cell position"
			     );
    };
  //
  for(long long int site_index=0;site_index<number_of_sites;site_index++)
    {
      configuration[site_index]=site_to_cell(
					     model,
					     cell_system,
					     cell_type_system,
					     site_system,
					     site_index,
					     "checkerboard"
					     );
    };
  // now construction for randomization of cell position
  random_shuffle_for_cells(
			   configuration,
			   simulation
			   );
  //
};
//
long int initial_state_system_class::site_to_cell(
						  const model_parameters_cellular_potts_class & model,
						  const cell_system_class & cell_system,
						  const type_system_class & cell_type_system,
						  const site_system_class & site_system,
						  const long long int & site_index,
						  const std::string & pattern
						  )
{
  std::vector<long int> coordinates;
  long int return_value;
  int space_dimension=model.get_space_dimension();
  unsigned long work_unsignedlong=(unsigned long)space_dimension;
  int component_index;
  boost::dynamic_bitset<> out_of_range_flag; 
  out_of_range_flag.resize(work_unsignedlong);
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      coordinates.push_back(0);
    };
  //
  coordinates=site_system.get_coordinates(site_index);
  //
  if(pattern=="checkerboard")
    {
      std::vector<long int> cell_coordinates;
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  cell_coordinates.push_back(0);
	};
      out_of_range_flag.reset();
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  coordinates[component_index]
	    =coordinates[component_index]
	    -generation_origin[component_index];
	  if(coordinates[component_index]<0)
	    {
	      out_of_range_flag.set(component_index);
	    } 
	  else if(
		  coordinates[component_index]
		  >=
		  (initial_cell_dimensions[component_index]
		   +initial_cell_separations[component_index])
		  *generation_array[component_index]
		  )
	    {
	      out_of_range_flag.set(component_index);
	    } 
	  else 
	    {
	      cell_coordinates[component_index]
		=(long int)floor(
				 (double)coordinates[component_index]
				 /
				 (double)(
					  initial_cell_dimensions[component_index]
					  +initial_cell_separations[component_index]
					  )
				 );
	      if(
		 coordinates[component_index]
		 -cell_coordinates[component_index]
		 *(
		   initial_cell_dimensions[component_index]
		   +initial_cell_separations[component_index]
		   )
		 > initial_cell_dimensions[component_index]
		 )
		{
		  out_of_range_flag.set(component_index);
		}
	    };
	  if(out_of_range_flag.any())
	    {
	      return_value=cell_system.get_buffer_cell();
	    } 
	  else 
	    {
	      return_value = cell_coordinate_to_index(model,cell_system,cell_coordinates);
	    };
	};
    }else{
      io_cellular_potts io_method;
      io_method.error_output(
			     "state_system_class",
			     "site_to_cell",
			     "Undifined initial cell pattern input."
			     );
      return_value=-1;
  };
  return return_value;
};
//
long long int initial_state_system_class::cell_coordinate_to_index(
								   const model_parameters_cellular_potts_class & model,
								   const cell_system_class & cell_system,
								   const std::vector<long int> & cell_coordinates
								   ) const 
{
  int component_index;
  long long int dimension_size=1;
  int work_int;
  long long int return_value=0;
  std::string out_of_system="off";
  int space_dimension=model.get_space_dimension();
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      work_int=space_dimension-component_index-1;
      return_value += (long long int)cell_coordinates.at(work_int)*dimension_size;
      dimension_size = dimension_size * (long long int)(generation_array[work_int]);
    };
  if (return_value<(long long int)cell_indicators[1]-(long long int)cell_indicators[0]+1)
    {
      return return_value+cell_indicators[0];
    }
  else
    {
      return cell_system.get_buffer_cell();
    };
};
//
void initial_state_system_class::random_polar(
					      std::vector<double> & cell_polarities,
					      const model_parameters_cellular_potts_class & model,
					      const cell_system_class       & cell_system,
					      const type_system_class       & cell_type_system,
					      const site_system_class       & site_system,
					      simulation_system_class & simulation
					      )
{
  int space_dimension=model.get_space_dimension();
  long int cell_index;
  int component_index;
  std::vector<double> work_vector(space_dimension,-1);
  long int number_of_cells=model.get_number_of_cells();
  toolbox tool;
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  work_vector[component_index]=2.0*simulation.random_generator_for_initializer()-1.0;
	}
	//      work_vector=simulation.get_initializer_random_number((long long int)space_dimension);
      tool.normalize(work_vector,1.0);
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  cell_polarities[cell_index*space_dimension+component_index]=work_vector[component_index];
	};
    };
};
//
void initial_state_system_class::random_shuffle_type_of_cell(
							     const model_parameters_cellular_potts_class & model,
							     cell_system_class & cell_system,
							     const type_system_class       & cell_type_system,
							     simulation_system_class & simulation						   
							     )
  const {
  long int cell_index;
  std::vector<int> work_types(number_of_cells,0);
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      work_types[cell_index]=cell_system.get_type(cell_index);
    };
  shuffle_table_generator_for_type(model,cell_type_system,simulation,work_types);
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      cell_system.set_type(cell_index,work_types[cell_index]);
    };
};
//
void initial_state_system_class::random_shuffle_for_cells(
							  std::vector<long int> & configuration,
							  simulation_system_class & simulation
							  )
  {
    std::vector<long int>::iterator site_index=configuration.begin();
    long int max_value_of_cells=0;
    while(site_index!=configuration.end())
      {
	if(max_value_of_cells<(*site_index))
	  {
	    max_value_of_cells=(*site_index);
	  }
	site_index++;
      };
    //
    std::vector<long int> shuffle_table;
    long int pivot;
    for(pivot=0;pivot<max_value_of_cells;pivot++)
      {
	shuffle_table.push_back(0);
      };
    //
    shuffle_table_generator_long_integer(
					 simulation,
					 shuffle_table
					 );
    //
    shuffle_table.push_back(max_value_of_cells);
    //
    site_index=configuration.begin();
    long int work_integer;
    //
    while(site_index!=configuration.end())
      {
	work_integer=shuffle_table.at(*site_index);
	(*site_index)=work_integer;
	site_index++;
      };
  };
//
void initial_state_system_class::shuffle_table_generator_long_integer(
								      simulation_system_class & simulation,
								      std::vector<long int> & shuffle_table
								      )
{
  std::list<long int> work_numbers;
  long int pivot=0;;
  long int target;
  //
  std::vector<long int>::iterator work_iterator=shuffle_table.begin();
  while(work_iterator!=shuffle_table.end())
    {
      work_numbers.push_back(pivot);
      pivot++;
      work_iterator++;
    };
  //
  work_iterator=shuffle_table.begin();
  std::list<long int>::iterator sub_work_iterator;
  while(work_iterator!=shuffle_table.end())
    {
      sub_work_iterator=work_numbers.begin();
      target=(int)(simulation.random_generator_for_initializer()*(double)work_numbers.size());
      for(pivot=0;pivot<target;pivot++)
	{
	  sub_work_iterator++;
	};
      (*work_iterator)=(*sub_work_iterator);
      work_numbers.erase(sub_work_iterator);
      work_iterator++;
    }
};
//
void initial_state_system_class::shuffle_table_generator_for_type(
								  const model_parameters_cellular_potts_class & model,
								  const type_system_class       & cell_type_system,
								  simulation_system_class & simulation,
								  std::vector<int> & shuffle_table
								  )
  const {
  std::list<int> work_cell_types(shuffle_table.size(),0);
  std::list<int>::iterator work_iterator;
  long long int cell_index;
  long long int pivot;
  long long int target;
  for(cell_index=0;cell_index<(int)shuffle_table.size();cell_index++)
    {
      work_cell_types.push_back(shuffle_table[cell_index]);
    };
  for(cell_index=0;cell_index<(int)shuffle_table.size();cell_index++)
    {
      target=(int)(simulation.random_generator_for_initializer()*(double)work_cell_types.size());
      work_iterator=work_cell_types.begin();
      for(pivot=0;pivot<target;pivot++)
	{
	  work_iterator++;
	};
      shuffle_table[cell_index]=(*work_iterator);
      work_cell_types.erase(work_iterator);
    };
  
}
//
void initial_state_system_class::generate_configuration(
							std::vector<long int> & configuration,
							const model_parameters_cellular_potts_class & model,
							const cell_system_class & cell_system,
							const type_system_class & cell_type_system,
							const site_system_class & site_system,
							simulation_system_class & simulation
							)
{
  fill_buffer(configuration);
  //
  for(int region_index=0;region_index<number_of_regions;region_index++)
    {
      generate_for_each_region(
			       configuration,
			       model,
			       cell_system,
			       cell_type_system,
			       site_system,
			       simulation,
			       region_index
			       );
    };
};
//
void initial_state_system_class::fill_buffer(
					     std::vector<long int> & configuration
					     )
{
  for(long long int site_index=0; site_index<(long long int)configuration.size();site_index++)
    {
      configuration[site_index]=buffer_cell;
    };
};
//
void initial_state_system_class::generate_for_each_region(
							  std::vector<long int> & configuration,
							  const model_parameters_cellular_potts_class & model,
							  const cell_system_class & cell_system,
							  const type_system_class & cell_type_system,
							  const site_system_class & site_system,
							  simulation_system_class & simulation,
							  const int & region_index
							  )
{
  random_checkerboard(
		      configuration,
		      model,
		      cell_system,
		      cell_type_system,
		      site_system,
		      simulation,
		      region_index
		      );
};
//
void initial_state_system_class::set_structure_region_item(
							   const int & region_index,
							   std::string child
							   )
{
  structure_item=io_method.generate_structure(
					      "configuration_setting.region",
					      region_index,
					      child
					      );
};
  /*======================
    Constructor for control system
   =======================*/
initial_state_system_class::initial_state_system_class(
						       const model_parameters_cellular_potts_class & model,
						       const cell_system_class & cell_system
						       )
{
  space_dimension=model.get_space_dimension();
  number_of_cells=model.get_number_of_cells();
  number_of_sites=model.get_number_of_sites();
  buffer_cell=cell_system.get_buffer_cell();
  system_dimensions.clear();
  for(int direction_index=0;direction_index<space_dimension;direction_index++)
    {
      system_dimensions.push_back(model.get_system_dimension(direction_index));
    };
};
