#include "cellular_potts_state.hpp"
  /*======================
    Methods
   =======================*/
//input methods
void state_system_class::initialize_state(
					  const model_parameters_cellular_potts_class & model,
					  cell_system_class & cell_system,
					  const type_system_class & cell_type_system,
					  site_system_class & site_system,
					  simulation_system_class & simulation
					  )
{
  std::string error;
  io_cellular_potts io_method;
  io_method.standard_output("*=== State is initialized !");
  input_initial_configuration_setting();
  generate_configuration(
			 model,
			 cell_system,
			 cell_type_system,
			 site_system,
			 simulation
			 );
  // checker
  //
  error=state_consistency_checker(
				  model,
				  cell_system,
				  cell_type_system,
				  site_system,
				  "zero_volume"
				  );
 //fprintf(stderr,"ok?");
 //
  error=state_consistency_checker(
				  model,
				  cell_system,
				  cell_type_system,
				  site_system,
				  "nearst_connection"
				  );
  // macroscopic quantities
  cell_system.calculate_cell_volumes(
				     model,
				     configuration,
				     cell_volumes
				     );
  //
  cell_scale=pow(cell_type_system.calculate_cell_scale(),1.0/(double)(model.get_space_dimension()));
  //;
  generate_polarity(
		    model,
		    cell_system,
		    cell_type_system,
		    site_system,
		    simulation
		    );
  //
  cell_system.assign_cell_origins(
				  model,
				  configuration,
				  cell_origins
				  );
  //
  //  fprintf(stderr,"ok in?\n");
  //
  cell_system.calculate_cell_total_positions(
					     model,
					     site_system,
					     configuration,
					     cell_origins,
					     cell_total_positions
					     );
  //
  cell_system.show_cell_positions(
				  model,
				  site_system,
				  configuration,
				  cell_origins,
				  cell_total_positions,
				  cell_volumes
				  );
  //  fprintf(stderr,"ok out?\n");
  // set random number
  allocate_random_number_memory(
				model,
				simulation
				);
};
//
void state_system_class::advance_time(
				      const model_parameters_cellular_potts_class & model,
				      cell_system_class & cell_system,
				      const type_system_class & cell_type_system,
				      const adhesion_system_class & adhesion_system,
				      site_system_class & site_system,
				      region_system_class & region_system,
				      simulation_system_class & simulation
				      )
{
  generate_random_number(model,simulation);
  // debug
  //clock_t start_time, end_time, total_time;
  // clock_t sub_start_time, sub_end_time, sub_total_time;
  //  long long int counter;
  //
  long long int number_of_flips=simulation.get_number_of_flips();
  long long int flip_index=0;
  long long int candidate_site=0;
  int candidate_neighbor=0;
  long int candidate_cell;
  //  fprintf(stderr,"check_0;%lld",flip_index);
  cell_system.assign_cell_origins(
				  model,
				  configuration,
				  cell_origins
				  );
  //
  //  fprintf(stderr,"ok in?\n");
  //
  cell_system.calculate_cell_total_positions(
					     model,
					     site_system,
					     configuration,
					     cell_origins,
					     cell_total_positions
					     );
  //
  local_state_class present_state(
				  model,
				  site_system,
				  cell_system,
				  cell_type_system
				  );
  //
  local_state_class candidate_state(
				    model,
				    site_system,
				    cell_system,
				    cell_type_system
				    );
  //
  polarity_motion_class polarity_motion(
					model,
					cell_system,
					cell_type_system
					);
  // memolize cell center positions in class polarity motion
  polarity_motion.initialize_polarities(
					model,
					cell_system,
					site_system,
					cell_polarities,
					cell_total_positions,
					cell_volumes,
					cell_origins
					);
  //
  copy_volume(
  			 cell_volumes,
			 original_cell_volumes
			 );
  /*
  // debug begin
  cell_system.check_cell_weight_polarity(
     	 model,
         site_system,
         configuration,
         cell_origins,
         cell_total_positions,
         original_cell_volumes
         );
  // debug end
  */
  //io_method.standard_output("debug3"+io_method.longlongint_to_string(cell_volumes[0]));
  //io_method.standard_output("debug3"+io_method.longlongint_to_string(original_cell_volumes[0]));
  //
  //
  init_work();
  //
  //counter=0;
  //total_time=clock()-clock();
  /*
  candidate_state.get_cpu_time(
			       cumulative_time,
			       "init"
			       );
  */
  //
  for(flip_index=0;flip_index<number_of_flips;flip_index++)
    {
      //
      candidate_site=random_for_site_choice[flip_index];
      candidate_neighbor=random_for_neighbor_choice[flip_index];
      //
      candidate_cell=configuration[
				   site_system.get_neighbor_site(
								 candidate_site,
								 candidate_neighbor
								 )
				   ];
      //
      if(
	 configuration[candidate_site]!=candidate_cell
	 &&
	 cell_system.mobile_table[configuration[candidate_site]]
	 &&
	 cell_system.mobile_table[candidate_cell]
	 )
	{
	  // debug
	  //start_time=clock();
	  //counter=counter+1;
	  //
	  present_state.initialize_site(
					candidate_site,
					configuration,
					cell_polarities,
					external_field,
					finite_field_flag,
					model,
					cell_system,
					cell_type_system,
					site_system
					);
	  candidate_state.initialize_site(
					  candidate_site,
					  configuration,
					  cell_polarities,
					  external_field,
					  finite_field_flag,
					  model,
					  cell_system,
					  cell_type_system,
					  site_system
					  );
	  //
	  present_state.initialize_local_state(
					       configuration[candidate_site],
					       configuration,
					       cell_volumes,
					       original_cell_volumes,
					       cell_polarities,
					       cell_total_positions,
					       cell_origins,
					       model,
					       cell_system,
					       cell_type_system,
					       site_system
					       );
	  //
	  candidate_state.initialize_local_state(
						 candidate_cell,
						 configuration,
						 cell_volumes,
						 original_cell_volumes,
						 cell_polarities,
						 cell_total_positions,
						 cell_origins,
						 model,
						 cell_system,
						 cell_type_system,
						 site_system
						 );
	  //
	  present_adhesion_energy=present_state.get_local_isotropic_adhesion_energy();
	  candidate_adhesion_energy=candidate_state.get_local_isotropic_adhesion_energy();
	  //
	  present_state.set_product_polarity();
	  candidate_state.set_product_polarity();
	  // polar adhesion
	  present_adhesion_energy=present_adhesion_energy+present_state.get_local_dipolar_adhesion_energy();
	  //
	  /*
	  candidate_state.get_cpu_time(
				       cumulative_time,
				       "start"
				       );
	  */
	  //
	  candidate_adhesion_energy=candidate_adhesion_energy+candidate_state.get_local_dipolar_adhesion_energy();
	  //
	  /*
	  candidate_state.get_cpu_time(
				       cumulative_time,
				       "add"
				       );
	  */
	  // debug
	  present_adhesion_energy+=present_state.get_local_adhesion_energy(adhesion_system);
	  candidate_adhesion_energy+=candidate_state.get_local_adhesion_energy(adhesion_system);
	  //
	  polarity_driving_difference=present_state.get_local_polarity_driving_energy()
	    -candidate_state.get_local_polarity_driving_energy();
	  //
	  // external_field
	  present_state.set_product_external_field();
	  candidate_state.set_product_external_field();
	  field_driving_difference=present_state.get_local_field_driving_energy()
	    -candidate_state.get_local_field_driving_energy();
	  //
	  //
	  present_volume_energy_difference=present_state.get_volume_energy_difference(configuration);
	  candidate_volume_energy_difference=candidate_state.get_volume_energy_difference(configuration);
	  //
	  //	  present_state.show_local_site_list(
	  //					       model,
	  //					       site_system,
	  //					       configuration
	  //					     );
	  //	  candidate_state.show_local_site_list(
	  //					       model,
	  //					       site_system,
	  //					       configuration
	  //					       );
	  //
	  update_flag=Metropholis_check_update(
					       get_energy_difference(),
					       random_for_trial[flip_index]
					       ) ;
	  //
	  if(update_flag==success_return) 
	    {
	      update(model,site_system,candidate_site,candidate_cell);
	      // update work
	    };
	  //
	  //end_time=clock();
	  //total_time=total_time+(end_time-start_time);
	};
      //
    };
  //
  //  fprintf(stderr,"tot: mcs/time[sec]: %lld/%f\n",counter,(double)(total_time)/CLOCKS_PER_SEC);
  //
  /*
  candidate_state.get_cpu_time(
			       cumulative_time,
			       "get"
			       );
  fprintf(stderr,"sub: mcs/time[sec]: %lld/%f\n",counter,(double)(cumulative_time)/CLOCKS_PER_SEC);
  */
  //
  cell_system.calculate_cell_total_positions(
					     model,
					     site_system,
					     configuration,
					     cell_origins,
					     cell_total_positions
					     );
  //
  get_adhesion_field(
		     model,
		     cell_system,
		     site_system,
		     adhesion_system
		     );
  //
  polarity_motion.update_polarities(
				    model,
				    cell_system,
				    site_system,
				    cell_total_positions,
				    cell_volumes,
				    cell_origins,
				    adhesion_field,
				    cell_polarities
				    );
  //
      //
  polarity_motion.get_displacements(
				    model,
				    cell_system,
				    site_system,
				    cell_total_positions,
				    cell_volumes,
				    cell_origins,
				    cell_displacements
				    );
};
//
void state_system_class::input_initial_configuration_setting()
  {
    io_cellular_potts io_method;
    /* for configuration */
    generation_method_for_configuration=io_method.get_input_string(
						 "site_setting_input" ,
						 "configuration_setting.method"
						 );
    //
    generation_type_for_configuration=io_method.get_input_string(
					       "site_setting_input" ,
					       "configuration_setting.type"
					       );
    //
    io_method.get_input_longint_array(
				      "site_setting_input" ,
				      "configuration_setting.origin",
				      generation_origin
				      );
    //
    io_method.get_input_longint_array(
				      "site_setting_input" ,
				      "configuration_setting.array",
				      generation_array
				      );
    //
    io_method.get_input_longint_array(
				      "site_setting_input" ,
				      "configuration_setting.cell_dimensions",
				      initial_cell_dimensions
				      );
    //
    io_method.get_input_longint_array(
				      "site_setting_input" ,
				      "configuration_setting.cell_separations",
				      initial_cell_separations
				      );
    /*for polarity*/
    generation_method_for_polarity=io_method.get_input_string(
							      "site_setting_input" ,
							      "polarity_setting.method"
							      );
    //
    generation_type_for_polarity=io_method.get_input_string(
					       "site_setting_input" ,
					       "polarity_setting.type"
					       );
    //
  };
void state_system_class::generate_configuration(
						const model_parameters_cellular_potts_class & model,
						const cell_system_class & cell_system,
						const type_system_class & cell_type_system,
						site_system_class & site_system,
						simulation_system_class & simulation
						)
{
  if(generation_method_for_configuration=="generation")
    {
      if(generation_type_for_configuration=="default")
	generation_type_for_configuration="random_checkerboard";
      if(generation_type_for_configuration=="random_checkerboard") 
	{
	  random_checkerboard(
			      model,
			      cell_system,
			      cell_type_system,
			      site_system,
			      simulation
			      );
	}
      else if(generation_type_for_configuration=="normal_checkerboard")
	{
	  normal_checkerboard(
			      model,
			      cell_system,
			      cell_type_system,
			      site_system,
			      simulation
			      );
	}
      else
	{
	    io_method.error_output(
				   "state_system_class",
				   "generate_configuration",
				   "undifined generation_type_for_configuration"
				   );
	};
    }
  else if(generation_method_for_configuration=="input")
    {
      configuration_light_read(model, site_system);
    }
  else
    {
      io_method.error_output(
			     "state_system_class",
			     "generate_configuration",
			     "undifined generation_method_for_configuration"
			     );
    };
};
//
void state_system_class::load_configuration(
					    const model_parameters_cellular_potts_class & model,
					    const site_system_class & site_system
					    )
{
  long long int site_index;
  long int cell_index;
  int direction_index;
  long long int pivot_site_index;
  std::vector<long int> load_coordinate(space_dimension,0);
  std::string structure_item;
  for(site_index=0;site_index<number_of_sites;site_index++)
    {
      set_structure_site_item(
			      site_index,
			      "cell",
			      structure_item
			      );
      cell_index=io_method.get_input_longint(
					     "configuration_input",
					     structure_item
					     );
      //
      set_structure_site_item(
			      site_index,
			      "coordinates",
			      structure_item
			      );
      load_coordinate.clear();
      io_method.get_input_longint_array(
					"configuration_input",
					structure_item,
					load_coordinate
					);
      //
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  load_coordinate[direction_index]=load_coordinate[direction_index]
	    +load_configuration_origin[direction_index];
	};
      //
      pivot_site_index=site_system.coordinate_to_site(
						      model,
						      load_coordinate
						      );
      //
      configuration[site_index]=cell_index;
    };
};
//
void state_system_class::set_structure_site_item(
						 const long long int & site_index,
						 const std::string & child,
						 std::string & structure_item
						 )
{
  structure_item=io_method.generate_structure(
					      "configuration.site",
					      site_index,
					      child
					      );
};
//
void state_system_class::normal_checkerboard(
					     const model_parameters_cellular_potts_class & model,
					     const cell_system_class & cell_system,
					     const type_system_class & cell_type_system,
					     const site_system_class & site_system,
					     simulation_system_class & simulation
					     )
{
  std::vector<long int> system_dimensions;
  int space_dimension=model.get_space_dimension();
  int direction_index;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      system_dimensions.push_back(model.get_system_dimension(direction_index));
    };
  //
  std::vector<long int> cell_region_end;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
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
			     "normal_checkerboard",
			     "due to out of range of initial cell position"
			     );
    };
  //
  long long int site_index;
  long long int number_of_sites=model.get_number_of_sites();
  for(site_index=0;site_index<number_of_sites;site_index++)
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
};
//
void state_system_class::random_checkerboard(
					     const model_parameters_cellular_potts_class & model,
					     const cell_system_class & cell_system,
					     const type_system_class & cell_type_system,
					     const site_system_class & site_system,
					     simulation_system_class & simulation
					     )
{
  std::vector<long int> system_dimensions;
  int space_dimension=model.get_space_dimension();
  int direction_index;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      system_dimensions.push_back(model.get_system_dimension(direction_index));
    };
  //
  std::vector<long int> cell_region_end;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
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
  long long int site_index;
  long long int number_of_sites=model.get_number_of_sites();
  for(site_index=0;site_index<number_of_sites;site_index++)
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
			   simulation
			   );
  //
};
//
long int state_system_class::site_to_cell(
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
long long int state_system_class::cell_coordinate_to_index(
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
  int number_of_cells=model.get_number_of_cells();
  if (return_value<number_of_cells)
    {
    return return_value;
    }else{
    return cell_system.get_buffer_cell();
  };
};
// output methods
void state_system_class::plot_mid_state(
					const model_parameters_cellular_potts_class & model,
					const cell_system_class & cell_system,
					const type_system_class & cell_type_system,
					const site_system_class & site_system,
					const long int & time_index,
					const int & sweep_step
					) 
{
  output_configuration(
		       model,
		       cell_system,
		       site_system,
		       time_index,
		       sweep_step
		       );
  output_polarity(
		  model,
		  site_system,
		  cell_system,
		  time_index,
		  sweep_step
		  );
  output_displacement(
		      model,
		      site_system,
		      cell_system,
		      time_index,
		      sweep_step
		      );
};
//
void state_system_class::finalize_state(
					const model_parameters_cellular_potts_class & model,
					const cell_system_class & cell_system,
					const type_system_class & cell_type_system,
					const site_system_class & site_system
					) 
{
  output_configuration(
		       model,
		       cell_system,
		       site_system,
		       state_finilization_step,
		       state_finilization_step
		       );
  std::string file_header="cell_position_";
  output_cell_position(
		       file_header,
		       state_finilization_step,
		       state_finilization_step
		       );
  output_polarity(
		  model,
		  site_system,
		  cell_system,
		  state_finilization_step,
		  state_finilization_step
		  );
  output_displacement(
		      model,
		      site_system,
		      cell_system,
		      state_finilization_step,
		      state_finilization_step
		      );
};
//
void state_system_class::input_configuration_output_setting()
{
    io_cellular_potts io_method;
    io_method.get_input_string_array(
				     "model_input" ,
				     "io_mode.configuration_output_type",
				     configuration_output_type
				     );
    io_method.get_input_string_array(
				     "model_input" ,
				     "io_mode.polarity_output_type",
				     polarity_output_type
				     );
};
//
void state_system_class::generate_polarity(
					   const model_parameters_cellular_potts_class & model,
					   const cell_system_class       & cell_system,
					   const type_system_class       & cell_type_system,
					   site_system_class       & site_system,
					   simulation_system_class & simulation
					   )
{
  if(generation_method_for_polarity=="generation")
    {
      if(generation_type_for_polarity=="default") 
	{
	  //    generation_type_for_polarity="random";
	  generation_type_for_polarity=default_type_for_polarity.c_str();
	}
      //fprintf(stderr,"type:%s,%s\n",generation_method_for_polarity.c_str(),default_type_for_polarity.c_str());
      if(generation_type_for_polarity=="random") random_polar(
							      model,
							      cell_system,
							      cell_type_system,
							      site_system,
							      simulation
							      );
      //      fprintf(stderr,"type:%s,%s\n",generation_method_for_polarity.c_str(),default_type_for_polarity.c_str());
      //abort();
    };
  if(generation_method_for_polarity=="input")
    {
      //      load_polarity(model,site_system);
      polarity_light_read(model,site_system);
    };
};
//
void state_system_class::load_polarity(
				       const model_parameters_cellular_potts_class & model,
				       const site_system_class & site_system
				       )
{
  long int cell_index;
  long int input_number_of_cells;
  int direction_index;
  std::string structure_item;
  std::vector<double> work_vector_polarity;
  //
  input_number_of_cells=io_method.get_input_longint(
						    "polarity_input",
						    "polarity.number_of_cells"
						    );
  //
  for(cell_index=0;cell_index<input_number_of_cells;cell_index++)
    {
      set_structure_polarity_item(
				  cell_index,
				  "components",
				  structure_item
				  );
      work_vector_polarity.clear();
      io_method.get_input_double_array(
				       "polarity_input",
				       structure_item,
				       work_vector_polarity
				       );
      //
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  cell_polarities[cell_index*space_dimension+direction_index]
	    =work_vector_polarity[direction_index];
	};
      //
    };
};
//
void state_system_class::set_structure_polarity_item(
						     const long int & cell_index,
						     const std::string & child,
						     std::string & structure_item
						     )
{
  structure_item=io_method.generate_structure(
					      "polarity.cell",
					      cell_index,
					      child
					      );
};
//
void state_system_class::random_polar(
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
void state_system_class::random_shuffle_type_of_cell(
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
}
//
void state_system_class::random_shuffle_for_cells(
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
void state_system_class::shuffle_table_generator_long_integer(
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
void state_system_class::shuffle_table_generator_for_type(
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
void state_system_class::output_configuration(
					      const model_parameters_cellular_potts_class & model,
					      const cell_system_class & cell_system,
					      const site_system_class & site_system,
					      const long int & time_index,
					      const int & sweep_step
					      ) const
{
  long long int site_index;
  long long int number_of_sites=model.get_number_of_sites();
  int component_index;
  int space_dimension=model.get_space_dimension();
  std::string data_identifier;
  io_cellular_potts io_method;
  boost::property_tree::ptree instance_property_tree;
  int iotype_index;
  for(iotype_index=0;iotype_index<(int)configuration_output_type.size();iotype_index++)
    {
      if(configuration_output_type[iotype_index]=="xml"&&time_index==state_finilization_step)
	{
  //
	  //
	  for(component_index=0;component_index<space_dimension;component_index++)
	    {
	      data_identifier="configuration.system_dimensions.component";
	      instance_property_tree.add(data_identifier,model.get_system_dimension(component_index));
	    }
	  std::vector<long int> coordinates;
	  for(component_index=0;component_index<space_dimension;component_index++)
	    {
	      coordinates.push_back(0);
	    }
	  for(site_index=0;site_index<number_of_sites;site_index++)
	    {
	      data_identifier=io_method.generate_structure_longlongint(
								       "configuration.site",
								       site_index,
								       "cell"
								       );
	      instance_property_tree.put(data_identifier,configuration.at(site_index));
	      coordinates=site_system.get_coordinates(site_index);
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  data_identifier
		    =io_method.generate_structure_longlongint(
							      "configuration.site",
							      site_index,
							      "coordinates.component"
							      );
		  instance_property_tree.add(data_identifier,coordinates[component_index]);
		}
	    };
	  //
	  const int indent=2;
	  // for boost 1.66
	   boost::property_tree::xml_parser::xml_writer_settings<char> 
	     settings = 
	   boost::property_tree::xml_parser::xml_writer_make_settings<char>(' ', indent);
	   boost::property_tree::xml_parser::write_xml(
	  					      io_method.get_filename("configuration_output"),
						      instance_property_tree,
						      std::locale(),
						      settings
						      );
	   /*
	  // for boost 1.55
	  boost::property_tree::xml_parser::write_xml(
						      io_method.get_filename("configuration_output"),
						      instance_property_tree,
						      std::locale(),
						      boost::property_tree::xml_writer_make_settings<std::string>( 
														 ' ', 
														 indent
														 )
	 					      );*/
	  //
	};
      //
      std::string message;
      //
      if(configuration_output_type[iotype_index]=="gnuplot")
	{
	  int plane_dimension=site_system.get_plane_dimension();
	  if(space_dimension>=plane_dimension)
	    {
	      std::vector<long int> coordinates;
	      std::vector<long int> plane_coordinates;
	      std::vector<long int> inplane_coordinates;
	      std::vector<std::string> plane_identifier;
	      std::vector<std::string> outplane_identifier;
	      long int plane_index;
	      long int outplane_index;
	      long int inplane_index;
	      long int plane_size=1;
	      long long int site_index;
	      long long int pivot_site_index;
	      long int number_of_planes;
	      std::string message;
	      std::vector<std::string> array_data;
	      long int inplane_x_dimension=model.get_system_dimension(0);
	      long int inplane_y_dimension=model.get_system_dimension(1);
	      std::vector<std::vector<long int> > plotdata(inplane_x_dimension, std::vector<long int>(inplane_y_dimension));
	      std::string configuration_plot="confgplot";
	      gnuplot_driver gplot(configuration_plot);
	      long int number_of_cells=model.get_number_of_cells();
	      //
	      for(plane_index=0;plane_index<space_dimension;plane_index++)
		{
		  plane_identifier.push_back("off");
		};
	      plane_identifier[0]="on";
	      plane_identifier[1]="on";
	      //
	      for(plane_index=0;plane_index<space_dimension;plane_index++)
		{
		  outplane_identifier.push_back("on");
		};
	      outplane_identifier[0]="off";
	      outplane_identifier[1]="off";
	      //
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  if(plane_identifier[component_index]=="on") plane_size=plane_size*model.get_system_dimension(component_index);
		};
	      //
	      number_of_planes=1;
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  if(plane_identifier[component_index]=="off") number_of_planes=number_of_planes*model.get_system_dimension(component_index);
		};
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  coordinates.push_back(0);
		  array_data.push_back("");
		}
	      for(component_index=0;component_index<space_dimension-plane_dimension;component_index++)
		{
		  plane_coordinates.push_back(0);
		}
	      for(component_index=0;component_index<plane_dimension;component_index++)
		{
		  inplane_coordinates.push_back(0);
		}
	      for(plane_index=0;plane_index<number_of_planes;plane_index++)
		{
		  site_system.plane_to_sub_coordinates(
						       plane_index, 
						       model,
						       plane_identifier,
						       plane_coordinates
						       );
		  outplane_index=0;
		  for(component_index=0;component_index<space_dimension;component_index++)
		    {
		      if(plane_identifier[component_index]=="off") 
			{
			  coordinates[component_index]=plane_coordinates[outplane_index];
			  outplane_index++;
			};
		    };
		  for(site_index=0;site_index<plane_size;site_index++)
		    {
		      site_system.plane_to_sub_coordinates(
							   site_index, 
							   model,
							   outplane_identifier,
							   inplane_coordinates
							   );
		      inplane_index=0;
		      for(component_index=0;component_index<space_dimension;component_index++)
			{
			if(plane_identifier[component_index]=="on") 
			  {
			    coordinates[component_index]=inplane_coordinates[inplane_index];
			    inplane_index++;
			  };
			};
		      pivot_site_index=site_system.coordinate_to_site(
								      model,
								      coordinates
								      );
		      if(buffer_cell!=configuration[pivot_site_index])
			{
			  plotdata[coordinates[0]][coordinates[1]]
			    =(
			      configuration[pivot_site_index]
			      +2*number_of_cells
			      *(number_of_cell_types-cell_system.get_type(
									  configuration[pivot_site_index]
									  )-1
				)
			      );
			  //			  if(plotdata[coordinates[0]][coordinates[1]]=0)
			  //			    {
			  //  message=io_method.longint_to_string(configuration[pivot_site_index]);//
			  //			      io_method.error_output("","",message);
			  //}
			}
		      else 
			{
			  plotdata[coordinates[0]][coordinates[1]]
			    =2*number_of_cells*number_of_cell_types;
			};
		    };
		  gplot.output_pm3d_data(
					 plotdata,
					 plane_index,
					 time_index,
					 sweep_step,
					 model.get_system_dimension(0),
					 model.get_system_dimension(1)
					 );
		};
	      gplot.set_axis(
			     0,
			     0.0,
			     (double)model.get_system_dimension(0),
			     (double)model.get_system_dimension(0)/4.0,
			     "x"
			     );
	      gplot.set_axis(
			     1,
			     0.0,
			     (double)model.get_system_dimension(1),
			     (double)model.get_system_dimension(1)/4.0,
			     "y"
			     );
	      gplot.set_axis(
			     2,
			     0.0,
			     (double)(2*number_of_cells*number_of_cell_types),
			     (double)number_of_cells/4.0,
			     "cell"
			     );
	      gplot.make_plot_file(
				   number_of_planes-1,
				   "png",
				   "splot",
				   time_index,
				   sweep_step
				   );
	    };
	};
    };
};
//
void state_system_class::output_cell_position(
					      const std::string & file_header,
					      const long int & time_index,
					      const int & sweep_step
					      )
{
  std::string message;
  std::string file_name = file_header
    +io_method.longint_to_format_string(time_index,"%04d")
    +io_method.longint_to_format_string(sweep_step,"%04d")
    + ".txt";
  io_method.file_initialize(file_name);
  for(long int cell_index=0; cell_index<number_of_cells; cell_index++)
    {
      message =io_method.double_to_string(
					  (double)cell_total_positions[cell_index*space_dimension]
					  /(double)cell_volumes[cell_index]
					  );
      for(int component_index=1;component_index<space_dimension; component_index++)
	{
	  message+=" ";
	  message+=io_method.double_to_string(
					      (double)cell_total_positions[cell_index*space_dimension+component_index]
					      /(double)cell_volumes[cell_index]
					      );
	};
      message+="\n";
      io_method.output_message(message, file_name);
    };
}
//
//
/*
const void state_system_class::output_polarity(
					       const model_parameters_cellular_potts_class & model,
					       const site_system_class & site_system,
					       const cell_system_class & cell_system,
					       const long int & time_index
					       ) const
{
  long long int site_index;
  long int cell_index;
  long int number_of_cells=model.get_number_of_cells();
  int component_index;
  int space_dimension=model.get_space_dimension();
  std::string data_identifier;
  io_cellular_potts io_method;
  boost::property_tree::ptree instance_property_tree;
  int iotype_index;
  for(iotype_index=0;iotype_index<(int)polarity_output_type.size();iotype_index++)
    {
      if(polarity_output_type[iotype_index]=="xml"&&time_index==state_finilization_step)
	{
	  //
	  for(component_index=0;component_index<space_dimension;component_index++)
	    {
	      data_identifier="configuration.system_dimensions.component";
	      instance_property_tree.add(data_identifier,model.get_system_dimension(component_index));
	    }
	  std::vector<long int> coordinates;
	  for(component_index=0;component_index<space_dimension;component_index++)
	    {
	      coordinates.push_back(0);
	    }
	  for(site_index=0;site_index<number_of_sites;site_index++)
	    {
	      data_identifier=io_method.generate_structure_longlongint("configuration.site",site_index,"cell");
	      instance_property_tree.put(data_identifier,configuration.at(site_index));
	      coordinates=site_system.get_coordinates(site_index);
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  data_identifier=io_method.generate_structure_longlongint("configuration.site",site_index,"coordinates.component");
		  instance_property_tree.add(data_identifier,coordinates[component_index]);
		}
	    };
	  //
	  const int indent=2;
	  boost::property_tree::xml_parser::write_xml(
						      io_method.get_filename("configuration_output"),
						      instance_property_tree,
						      std::locale(),
						      boost::property_tree::xml_parser::xml_writer_make_settings( 
														 ' ', 
														 indent, 
														 boost::property_tree::xml_parser::widen<char>("utf-8")
														 )
						      );
	  //
	};
      if(configuration_output_type[iotype_index]=="gnuplot")
	{
	  int plane_dimension=site_system.get_plane_dimension();
	  if(space_dimension>=plane_dimension)
	    {
	      std::vector<long int> coordinates;
	      std::vector<long int> plane_coordinates;
	      std::vector<long int> inplane_coordinates;
	      std::vector<std::string> plane_identifier;
	      std::vector<std::string> outplane_identifier;
	      long int plane_index;
	      long int outplane_index;
	      long int inplane_index;
	      long int plane_size=1;
	      long long int site_index;
	      long long int pivot_site_index;
	      long int number_of_planes;
	      std::string message;
	      std::vector<std::string> array_data;
	      long int inplane_x_dimension=model.get_system_dimension(0);
	      long int inplane_y_dimension=model.get_system_dimension(1);
	      std::vector<std::vector<long int> > plotdata(inplane_x_dimension, std::vector<long int>(inplane_y_dimension));
	      std::string configuration_plot="confgplot";
	      gnuplot_driver gplot(configuration_plot);
	      long int number_of_cells=model.get_number_of_cells();
	      //
	      for(plane_index=0;plane_index<space_dimension;plane_index++)
		{
		  plane_identifier.push_back("off");
		};
	      plane_identifier[0]="on";
	      plane_identifier[1]="on";
	      //
	      for(plane_index=0;plane_index<space_dimension;plane_index++)
		{
		  outplane_identifier.push_back("on");
		};
	      outplane_identifier[0]="off";
	      outplane_identifier[1]="off";
	      //
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  if(plane_identifier[component_index]=="on") plane_size=plane_size*model.get_system_dimension(component_index);
		};
	      //
	      number_of_planes=1;
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  if(plane_identifier[component_index]=="off") number_of_planes=number_of_planes*model.get_system_dimension(component_index);
		};
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  coordinates.push_back(0);
		  array_data.push_back("");
		}
	      for(component_index=0;component_index<space_dimension-plane_dimension;component_index++)
		{
		  plane_coordinates.push_back(0);
		}
	      for(component_index=0;component_index<plane_dimension;component_index++)
		{
		  inplane_coordinates.push_back(0);
		}
	      for(plane_index=0;plane_index<number_of_planes;plane_index++)
		{
		  site_system.plane_to_sub_coordinates(
						       plane_index, 
						       model,
						       plane_identifier,
						       plane_coordinates
						       );
		  outplane_index=0;
		  for(component_index=0;component_index<space_dimension;component_index++)
		    {
		      if(plane_identifier[component_index]=="off") 
			{
			  coordinates[component_index]=plane_coordinates[outplane_index];
			  outplane_index++;
			};
		    };
		  for(site_index=0;site_index<plane_size;site_index++)
		    {
		      site_system.plane_to_sub_coordinates(
							   site_index, 
							   model,
							   outplane_identifier,
							   inplane_coordinates
							   );
		      inplane_index=0;
		      for(component_index=0;component_index<space_dimension;component_index++)
			{
			if(plane_identifier[component_index]=="on") 
			  {
			    coordinates[component_index]=inplane_coordinates[inplane_index];
			    inplane_index++;
			  };
			};
		      pivot_site_index=site_system.coordinate_to_site(
								      model,
								      coordinates
								      );
		      plotdata[coordinates[0]][coordinates[1]]=configuration[pivot_site_index];
		    };
		  gplot.output_pm3d_data(
					 plotdata,
					 plane_index,
					 time_index,
					 model.get_system_dimension(0),
					 model.get_system_dimension(1)
					 );
		};
	      gplot.set_axis(
			     0,
			     0.0,
			     (double)model.get_system_dimension(0),
			     (double)model.get_system_dimension(0)/4.0,
			     "x"
			     );
	      gplot.set_axis(
			     1,
			     0.0,
			     (double)model.get_system_dimension(1),
			     (double)model.get_system_dimension(1)/4.0,
			     "y"
			     );
	      gplot.set_axis(
			     2,
			     0.0,
			     (double)number_of_cells,
			     (double)number_of_cells/4.0,
			     "cell"
			     );
	      gplot.make_plot_file(
				   number_of_planes-1,
				   "png",
				   "splot",
				   time_index
				   );
	    };
	};
    };
};
//
//
*/
 //
void state_system_class::output_polarity(
					 const model_parameters_cellular_potts_class & model,
					 const site_system_class & site_system,
					 const cell_system_class & cell_system,
					 const long int & time_index,
					 const int & sweep_step
					 ) const
{
  long int cell_index;
  long int number_of_cells=model.get_number_of_cells();
  int component_index;
  int space_dimension=model.get_space_dimension();
  std::string data_identifier;
  io_cellular_potts io_method;
  boost::property_tree::ptree instance_property_tree;
  int iotype_index;
  std::vector<double> components(space_dimension,0.0);
  for(iotype_index=0;iotype_index<(int)polarity_output_type.size();iotype_index++)
    { 
      if(polarity_output_type[iotype_index]=="xml"&&time_index==state_finilization_step)
	{
  //
	  data_identifier="polarity.number_of_cells";
	  instance_property_tree.add(data_identifier,model.get_number_of_cells());
	  data_identifier="polarity.number_of_components";
	  instance_property_tree.add(data_identifier,model.get_space_dimension());
	  for(cell_index=0;cell_index<number_of_cells;cell_index++)
	    {
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  data_identifier
		    =io_method.generate_structure_longint(
							  "polarity.cell",
							  cell_index,
							  "components.component"
							  );
		  instance_property_tree.add(
					     data_identifier,
					     cell_polarities[cell_index*space_dimension+component_index]
					     );
		}
	    };
	  //
	  const int indent=2;
	  // for boost 1.66
	   boost::property_tree::xml_parser::xml_writer_settings<char> 
	    settings = 
	   boost::property_tree::xml_parser::xml_writer_make_settings<char> (' ', indent);
	   boost::property_tree::xml_parser::write_xml(
	  					      io_method.get_filename("polarity_output"),
	  				      instance_property_tree,
	  				      std::locale(),
	  				      settings
	  				      );
	  // for boost 1.55 
	  /*
	  boost::property_tree::xml_parser::write_xml(
						      io_method.get_filename("polarity_output"),
						      instance_property_tree,
						      std::locale(),
						      boost::property_tree::xml_parser::xml_writer_make_settings<std::string>( 
															      ' ', 2
															      //															      indent, 
															      //boost::property_tree::xml_parser::widen<char>("utf-8")
															       )
						      );*/
	};
      if(polarity_output_type[iotype_index]=="gnuplot")
	{
	  int plane_dimension=site_system.get_plane_dimension();
	  if(space_dimension>=plane_dimension)
	    {
	      std::vector<long int> coordinates;
	      std::vector<long int> plane_coordinates;
	      std::vector<long int> inplane_coordinates;
	      std::vector<std::string> plane_identifier;
	      std::vector<std::string> outplane_identifier;
	      long int plane_index;
	      long int outplane_index;
	      long int inplane_index;
	      long int plane_size=1;
	      long long int site_index;
	      long long int pivot_site_index;
	      long int number_of_planes;
	      std::string message;
	      std::vector<std::string> array_data;
	      std::string configuration_plot="confgplot";
	      gnuplot_driver gplot(configuration_plot);
	      long int number_of_cells=model.get_number_of_cells();
	      std::vector<std::vector<double> > arrows;
	      std::vector<std::vector<double> > cell_coordinates;
	      std::vector<long int> translate_cell_types;
	      std::vector<double> work_vector_double(space_dimension,0.0);
	      std::vector<int> cell_stacks(number_of_cells,0); 
	      //
	      for(plane_index=0;plane_index<space_dimension;plane_index++)
		{
		  plane_identifier.push_back("off");
		};
	      plane_identifier[0]="on";
	      plane_identifier[1]="on";
	      //
	      for(plane_index=0;plane_index<space_dimension;plane_index++)
		{
		  outplane_identifier.push_back("on");
		};
	      outplane_identifier[0]="off";
	      outplane_identifier[1]="off";
	      //
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  if(plane_identifier[component_index]=="on") plane_size=plane_size*model.get_system_dimension(component_index);
		};
	      //
	      number_of_planes=1;
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  if(plane_identifier[component_index]=="off") number_of_planes=number_of_planes*model.get_system_dimension(component_index);
		};
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  coordinates.push_back(0);
		  array_data.push_back("");
		}
	      for(component_index=0;component_index<space_dimension-plane_dimension;component_index++)
		{
		  plane_coordinates.push_back(0);
		}
	      for(component_index=0;component_index<plane_dimension;component_index++)
		{
		  inplane_coordinates.push_back(0);
		}
	      for(plane_index=0;plane_index<number_of_planes;plane_index++)
		{
      		  site_system.plane_to_sub_coordinates(
		  				       plane_index, 
		  				       model,
		  				       plane_identifier,
		  				       plane_coordinates
		  				       );
		  outplane_index=0;
		  //
		  for(component_index=0;component_index<space_dimension;component_index++)
		    {
		      if(plane_identifier[component_index]=="off") 
			{
			  coordinates[component_index]=plane_coordinates[outplane_index];
			  outplane_index++;
			};
		    };
		  // initialization
		  cell_coordinates.clear();
		  translate_cell_types.clear();
		  arrows.clear();
		  //
		  for(site_index=0;site_index<plane_size;site_index++)
		    {
		      site_system.plane_to_sub_coordinates(
							   site_index, 
							   model,
							   outplane_identifier,
							   inplane_coordinates
							   );
		      inplane_index=0;
		      for(component_index=0;component_index<space_dimension;component_index++)
			{
			if(plane_identifier[component_index]=="on") 
			  {
			    coordinates[component_index]=inplane_coordinates[inplane_index];
			    inplane_index++;
			  };
			};
		      pivot_site_index=site_system.coordinate_to_site(
								      model,
								      coordinates
								      );
		      if(0<=configuration[pivot_site_index]&&configuration[pivot_site_index]<number_of_cells)
			{
			  if(cell_stacks[configuration[pivot_site_index]]==0)
			    {
			      cell_stacks[configuration[pivot_site_index]]=1;
			      cell_system.calculate_cell_position_in_system(
									    model,
									    site_system,
									    cell_origins,
									    cell_total_positions,
									    cell_volumes,
									    configuration[pivot_site_index],
									    work_vector_double
									    );
			      cell_coordinates.push_back(work_vector_double);
			      translate_cell_types.push_back(cell_types[configuration[pivot_site_index]]);
			      //			      for(component_index=0;component_index<space_dimension;component_index++)
			      //	{
			      //	  cell_coordinates[configuration[pivot_site_index]][component_index]
			      //	    =work_vector_double[component_index];
			      //	};
			      for(component_index=0;component_index<space_dimension;component_index++)
				{
				  work_vector_double[component_index]
				    =cell_polarities[configuration[pivot_site_index]*space_dimension+component_index];
				};
			      arrows.push_back(work_vector_double);
			    };
			};
		    };
		  gplot.output_vectors_with_type_on_loading_file(
								 cell_coordinates,
								 arrows,
								 translate_cell_types,
								 plane_index,
								 time_index,
								 sweep_step,
								 space_dimension,
								 cell_scale
								 );
		};
	    };
	};
    };
};
 //
void state_system_class::output_displacement(
					     const model_parameters_cellular_potts_class & model,
					     const site_system_class & site_system,
					     const cell_system_class & cell_system,
					     const long int & time_index,
					     const int & sweep_step
					     ) 
{
  long int cell_index;
  long int number_of_cells=model.get_number_of_cells();
  int component_index;
  int space_dimension=model.get_space_dimension();
  std::string data_identifier;
  io_cellular_potts io_method;
  boost::property_tree::ptree instance_property_tree;
  int iotype_index;
  std::vector<double> components(space_dimension,0.0);
  for(long int cell_index=0;
      cell_index<number_of_cells;
      cell_index++)
    {
      for(int direction_index=0;
	  direction_index<space_dimension;
	  direction_index++)
	{
	  plot_displacements[cell_index*space_dimension+direction_index]=
	    (double)cell_total_positions[cell_index*space_dimension+direction_index]
	    /(double)cell_volumes[cell_index]
	    -plot_memory_position[cell_index*space_dimension+direction_index];
	};
    };
  for(long int cell_index=0;
      cell_index<number_of_cells;
      cell_index++)
    {
      for(int direction_index=0;
	  direction_index<space_dimension;
	  direction_index++)
	{
	  plot_memory_position[cell_index*space_dimension+direction_index]=
	    (double)cell_total_positions[cell_index*space_dimension+direction_index]
	    /(double)cell_volumes[cell_index];
	};
    };
  for(iotype_index=0;iotype_index<(int)polarity_output_type.size();iotype_index++)
    { 
      if(polarity_output_type[iotype_index]=="xml"&&time_index==state_finilization_step)
	{
  //
	  data_identifier="displacement.number_of_cells";
	  instance_property_tree.add(data_identifier,model.get_number_of_cells());
	  data_identifier="displacement.number_of_components";
	  instance_property_tree.add(data_identifier,model.get_space_dimension());
	  for(cell_index=0;cell_index<number_of_cells;cell_index++)
	    {
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  data_identifier
		    =io_method.generate_structure_longint(
							  "displacement.cell",
							  cell_index,
							  "components.component"
							  );
		  instance_property_tree.add(
					     data_identifier,
					     plot_displacements[cell_index*space_dimension+component_index]
					     );
		}
	    };
	  //
	  const int indent=2;
	  // for boost 1.66
	   boost::property_tree::xml_parser::xml_writer_settings<char> 
	    settings = 
	   boost::property_tree::xml_parser::xml_writer_make_settings<char> (' ', indent);
	   boost::property_tree::xml_parser::write_xml(
	  					      io_method.get_filename("displacement_output"),
	  				      instance_property_tree,
	  				      std::locale(),
	  				      settings
	  				      );
	  // for boost 1.55 
	  /*
	  boost::property_tree::xml_parser::write_xml(
						      io_method.get_filename("polarity_output"),
						      instance_property_tree,
						      std::locale(),
						      boost::property_tree::xml_parser::xml_writer_make_settings<std::string>( 
															      ' ', 2
															      //															      indent, 
															      //boost::property_tree::xml_parser::widen<char>("utf-8")
															       )
						      );*/
	};
      if(polarity_output_type[iotype_index]=="gnuplot")
	{
	  int plane_dimension=site_system.get_plane_dimension();
	  if(space_dimension>=plane_dimension)
	    {
	      std::vector<long int> coordinates;
	      std::vector<long int> plane_coordinates;
	      std::vector<long int> inplane_coordinates;
	      std::vector<std::string> plane_identifier;
	      std::vector<std::string> outplane_identifier;
	      long int plane_index;
	      long int outplane_index;
	      long int inplane_index;
	      long int plane_size=1;
	      long long int site_index;
	      long long int pivot_site_index;
	      long int number_of_planes;
	      std::string message;
	      std::vector<std::string> array_data;
	      std::string configuration_plot="displplot";
	      gnuplot_driver gplot(configuration_plot);
	      long int number_of_cells=model.get_number_of_cells();
	      std::vector<std::vector<double> > arrows;
	      std::vector<std::vector<double> > cell_coordinates;
	      std::vector<long int> translate_cell_types;
	      std::vector<double> work_vector_double(space_dimension,0.0);
	      std::vector<int> cell_stacks(number_of_cells,0); 
	      //
	      for(plane_index=0;plane_index<space_dimension;plane_index++)
		{
		  plane_identifier.push_back("off");
		};
	      plane_identifier[0]="on";
	      plane_identifier[1]="on";
	      //
	      for(plane_index=0;plane_index<space_dimension;plane_index++)
		{
		  outplane_identifier.push_back("on");
		};
	      outplane_identifier[0]="off";
	      outplane_identifier[1]="off";
	      //
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  if(plane_identifier[component_index]=="on") 
		    plane_size=plane_size*model.get_system_dimension(component_index);
		};
	      //
	      number_of_planes=1;
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  if(plane_identifier[component_index]=="off") 
		    number_of_planes
		      =number_of_planes*model.get_system_dimension(component_index);
		};
	      for(component_index=0;component_index<space_dimension;component_index++)
		{
		  coordinates.push_back(0);
		  array_data.push_back("");
		}
	      for(component_index=0;component_index<space_dimension-plane_dimension;component_index++)
		{
		  plane_coordinates.push_back(0);
		}
	      for(component_index=0;component_index<plane_dimension;component_index++)
		{
		  inplane_coordinates.push_back(0);
		}
	      for(plane_index=0;plane_index<number_of_planes;plane_index++)
		{
      		  site_system.plane_to_sub_coordinates(
		  				       plane_index, 
		  				       model,
		  				       plane_identifier,
		  				       plane_coordinates
		  				       );
		  outplane_index=0;
		  //
		  for(component_index=0;component_index<space_dimension;component_index++)
		    {
		      if(plane_identifier[component_index]=="off") 
			{
			  coordinates[component_index]
			    =plane_coordinates[outplane_index];
			  outplane_index++;
			};
		    };
		  // initialization
		  cell_coordinates.clear();
		  translate_cell_types.clear();
		  arrows.clear();
		  //
		  for(site_index=0;site_index<plane_size;site_index++)
		    {
		      site_system.plane_to_sub_coordinates(
							   site_index, 
							   model,
							   outplane_identifier,
							   inplane_coordinates
							   );
		      inplane_index=0;
		      for(component_index=0;component_index<space_dimension;component_index++)
			{
			if(plane_identifier[component_index]=="on") 
			  {
			    coordinates[component_index]
			      =inplane_coordinates[inplane_index];
			    inplane_index++;
			  };
			};
		      pivot_site_index=site_system.coordinate_to_site(
								      model,
								      coordinates
								      );
		      if(0<=configuration[pivot_site_index]&&configuration[pivot_site_index]<number_of_cells)
			{
			  if(cell_stacks[configuration[pivot_site_index]]==0)
			    {
			      cell_stacks[configuration[pivot_site_index]]=1;
			      cell_system.calculate_cell_position_in_system(
									    model,
									    site_system,
									    cell_origins,
									    cell_total_positions,
									    cell_volumes,
									    configuration[pivot_site_index],
									    work_vector_double
									    );
			      cell_coordinates.push_back(work_vector_double);
			      translate_cell_types.push_back(cell_types[configuration[pivot_site_index]]);
			      for(component_index=0;component_index<space_dimension;component_index++)
				{
				  work_vector_double[component_index]
				    =plot_displacements[configuration[pivot_site_index]*space_dimension
						     +component_index];
				};
			      arrows.push_back(work_vector_double);
			    };
			};
		    };
		  gplot.output_vectors_with_type_on_loading_file(
								 cell_coordinates,
								 arrows,
								 translate_cell_types,
								 plane_index,
								 time_index,
								 sweep_step,
								 space_dimension,
								 cell_scale
								 );
		};
	    };
	};
    };
};
//
std::string state_system_class::state_consistency_checker(
							  const model_parameters_cellular_potts_class & model,
							  const cell_system_class & cell_system,
							  const type_system_class & cell_type_system,
							  const site_system_class & site_system,
							  const std::string & job_flag
							  ) 
const {
  std::string return_value="ok";
  if(job_flag=="zero_volume")
    {
      long int cell_index;
      long long int site_index;
      long int number_of_cells=model.get_number_of_cells();
      long long int number_of_sites=model.get_number_of_sites();
      std::vector<int> counter;
      for(cell_index=0;cell_index<number_of_cells;cell_index++)
	{
	  counter.push_back(0);
	};
      for(site_index=0;site_index<number_of_sites;site_index++)
	{
	  if(
	     configuration[site_index] < number_of_cells
	     &&
	     configuration[site_index] >= 0
	     )
	  counter[configuration[site_index]]++;
	};
      for(cell_index=0;cell_index<number_of_cells;cell_index++)
	{
	  if((counter[cell_index]==0)&&(cell_system.get_type(cell_index)!=cell_type_system.get_buffer_type()))
	    {
	      io_cellular_potts io_method;
	      std::string message;
	      message  ="No volume of the cell ";
	      message +=io_method.longlongint_to_string(cell_index);
	      io_method.error_output(
				     "state_system_class",
				     "state_consistency_checker",
				     message
				     );
	      return_value="error";
	    }
	};
    };
  if(job_flag=="nearst_connection")
    {
      long long int site_index;
      int neighbor_index;
      long long int number_of_sites=model.get_number_of_sites();
      int coordinate_number=site_system.get_coordinate_number(model);
      long long int neighbor_site;
      int false_flag;
      for(site_index=0;site_index<number_of_sites;site_index++)
	{
	  false_flag=0;
	  for(neighbor_index=0;neighbor_index<coordinate_number;neighbor_index++)
	    {
	      neighbor_site=site_system.get_nearest_neighbor_site(model,site_index,neighbor_index);
	      if(configuration[neighbor_site]==configuration[site_index]) false_flag=1;
	    };
	  if(false_flag==0)
	    {
	      io_cellular_potts io_method;
	      std::string message;
	      message  ="Absence of same state neighbor at site ";
	      message +=io_method.longlongint_to_string(site_index);
	      message +="(state_system_class::state_consistency_checker)";
	      io_method.standard_output(message);
	      return_value="error";
	    };
	};
    };
  return return_value;
};				     
//
void state_system_class::allocate_random_number_memory(
						       const model_parameters_cellular_potts_class & model,
						       simulation_system_class & simulation
						       )
{
  long long int number_of_flips=simulation.get_number_of_flips();
  long long int random_index;
  for(random_index=0;random_index<number_of_flips;random_index++)
    {
      random_for_trial.push_back(0.0);
      random_for_site_choice.push_back(0);
      random_for_neighbor_choice.push_back(0);
    };
};
//
void state_system_class::generate_random_number(
						const model_parameters_cellular_potts_class & model,
						simulation_system_class & simulation
						)
{
  long long int number_of_random_numbers=simulation.get_number_of_flips();
  random_for_trial=simulation.get_trial_random_number(number_of_random_numbers);
  random_for_site_choice=simulation.get_site_choice_random_number(number_of_random_numbers);
  random_for_neighbor_choice=simulation.get_neighbor_choice_random_number(number_of_random_numbers);
};
//
//const std::vector<long long int> state_system_class::get_cell_volumes() 
//  const {
//  return cell_volumes;
//};
void state_system_class::get_cell_volumes(
					  std::vector<long long int> & volumes
					  ) 
  const {
  long int cell_index;
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      volumes[cell_index]=cell_volumes[cell_index];
    }
}; 
//
void  state_system_class::get_cell_polarities(
					      std::vector<double> & polarities 
					      ) 
  const {
  long int cell_index;
  int direction_index;
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  polarities[cell_index*space_dimension+direction_index]
	    =cell_polarities[cell_index*space_dimension+direction_index];
	};
    };
}; 
//
void state_system_class::get_sorted_cell_polarities(
						    std::vector<double> & polarities 
						    ) 
  const {
  long int cell_index;
  int direction_index;
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  polarities[direction_index*number_of_cells+cell_index]
	    =cell_polarities[cell_index*space_dimension+direction_index];
	};
    };
}; 
//
void state_system_class::get_cell_displacements(
						std::vector<double> & output_cell_displacements
						) 
  const {
  int index;
  for(index=0;index<number_of_cells*space_dimension;index++)
    {
       output_cell_displacements[index]=cell_displacements[index];
    };
}; 
//
void state_system_class::get_sorted_cell_displacements(
						       std::vector<double> & output_sorted_cell_displacements
						       ) 
  const {
  int cell_index;
  int direction_index;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      for(cell_index=0;cell_index<number_of_cells;cell_index++)
	{
	  output_sorted_cell_displacements[number_of_cells*direction_index+cell_index]
	    =cell_displacements[cell_index*space_dimension+direction_index];
	};
    };
}; 
//
void state_system_class::set_model_parameter(
					     const double & parameter,
					     const std::string & parameter_identifier
					     )
{
  if(parameter_identifier=="beta")
    {
      beta=parameter;
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "state_system_class",
			     "set_model_parameter",
			     "due to undefined input"
			     );
    };
};
  //
void state_system_class::set_model_parameter_array(
						   const std::vector<double> & parameters,
						   const std::string & parameter_identifier
						   )
{
      io_cellular_potts io_method;
      io_method.error_output(
			     "state_system_class",
			     "set_model_parameter",
			     "due to undefined input"
			     );
};
//
/*
inline double state_system_class::normalized_product(
						     const std::vector<double> & vector_1,
						     const std::vector<double> & vector_2
						     ) 
{
  //
  double norm_1     =std::inner_product(vector_1.begin(),vector_1.end(),vector_1.begin(),0.0);
  double norm_2     =std::inner_product(vector_2.begin(),vector_2.end(),vector_2.begin(),0.0);
  double product_12 =std::inner_product(vector_1.begin(),vector_1.end(),vector_2.begin(),0.0);
  return product_12/pow(norm_1*norm_2,0.5);
  //
};
*/
//
inline double state_system_class::external_product_sign(
							const std::vector<double> & vector_1,
							const std::vector<double> & vector_2
							) 
{
  //
  double external_product_12(0.0);
  if(space_dimension==2)
    {
      external_product_12 =vector_1[0]*vector_2[1]-vector_1[1]*vector_2[0];
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "state_system_class",
			     "product_sign",
			     "due to under construction!"
			     );
    }
  return copysign(1.0,external_product_12);
  //
};
//
void state_system_class::set_adhesion_coupling(
					       const type_system_class & cell_type_system
					       )
{
  int type_index;
  int neighbor_type_index;
  double work_double;
  for(
      type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      cell_type_system.get_adhesion_couplings(
					      type_index,
					      "isotropic",
					      work_vector_coupling
					      );
      // isotropic
      for(neighbor_type_index=0;
	  neighbor_type_index<number_of_cell_types;
	  neighbor_type_index++)
	{
	  isotropic_adhesion_couplings[type_index][neighbor_type_index]
	    =work_vector_coupling[neighbor_type_index];
	};
      // dipolar
      cell_type_system.get_adhesion_couplings(
					      type_index,
					      "dipolar",
					      work_vector_coupling
					      );
      for(neighbor_type_index=0;
	  neighbor_type_index<number_of_cell_types;
	  neighbor_type_index++)
	{
	  dipolar_adhesion_couplings[type_index][neighbor_type_index]
	    =work_vector_coupling[neighbor_type_index];
	};
      //
      cell_type_system.get_adhesion_basal(
					  type_index,
					  "dipolar",
					  work_double
					  );
      //
      dipolar_adhesion_basals[type_index]=work_double;
      // quadrapolar
      cell_type_system.get_adhesion_couplings(
					      type_index,
					      "quadrapolar",
					      work_vector_coupling
					      );
      for(neighbor_type_index=0;
	  neighbor_type_index<number_of_cell_types;
	  neighbor_type_index++)
	{
	  quadrapolar_adhesion_couplings[type_index][neighbor_type_index]
	    =work_vector_coupling[neighbor_type_index];
	};
      //
      cell_type_system.get_adhesion_basal(
					  type_index,
					  "quadrapolar",
					  work_double
					  );
      //
      quadrapolar_adhesion_basals[type_index]=work_double;
    };
};
//
void state_system_class::set_external_field(
					    const site_system_class & site_system,
					    region_system_class & region_system
					    )
{
  for(long long int site_index=0;site_index<number_of_sites;site_index++)
    {
      region_system.get_field(site_system,site_index,work_external_field);
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  external_field[site_index*space_dimension+direction_index]=work_external_field[direction_index];
	};
      if(tool.check_finite_vector(external_field)==tool.get_value("true")) finite_field_flag[site_index]=1;
    };
};
//
void state_system_class::copy_volume(
					const std::vector<long long int> & origin_volumes,
					std::vector<long long int> & copy_volumes
					)
{
	int counter=0;
	std::vector<long long int>::const_iterator cell_index=origin_volumes.begin();  
	while(cell_index!=origin_volumes.end())
	{
		copy_volumes[counter]=(*cell_index);
		//debug
		//io_method.standard_output("debug4:"+io_method.int_to_string(*cell_index));
		//
		counter++;
		cell_index++;
	};

};
//
void state_system_class::init_work()
{
  polarity_work=0.0;
  field_work=0.0;
}
//
void state_system_class::add_work()
{
  polarity_work-=polarity_driving_difference;
  field_work-=field_work+field_driving_difference;
}
  /*======================
    Constructor
   =======================*/
state_system_class::state_system_class(
				       const model_parameters_cellular_potts_class & model,
				       const site_system_class & site_system,
				       const cell_system_class & cell_system,
				       const type_system_class & cell_type_system
				       )
{
  // allocate configuration memory
  true_value=1;
  false_value=0;
  number_of_sites=model.get_number_of_sites();
  number_of_cells=model.get_number_of_cells();
  number_of_cell_types=model.get_number_of_cell_types();
  number_of_neighbor_sites=site_system.get_number_of_neighbor_sites();
  long long int init_site_number=(long long int)configuration.size();
  long long int site_index;
  space_dimension=model.get_space_dimension();
  buffer_cell=cell_system.get_buffer_cell();
  buffer_type=cell_type_system.get_buffer_type();
  for(site_index=0;site_index<number_of_sites;site_index++)
    {
      if(site_index>=init_site_number)
	{
	  configuration.push_back(0);
	}else{
	configuration[site_index]=0;
      };
    };
  int direction_index;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      load_configuration_origin.push_back(0);
    };
  // allocate cell order parameters
  long int number_of_cells=model.get_number_of_cells();
  long int init_cell_number=cell_volumes.size();
  long int cell_index;
  int type_index;
  double work_double;
  original_cell_volumes.clear();
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      if(cell_index>=init_cell_number)
	{
	  cell_volumes.push_back(-1);
	  original_cell_volumes.push_back(-1);
	}else{
	cell_volumes[cell_index]=-1;
	original_cell_volumes.push_back(-1);
      };
    };
  init_cell_number=cell_perimaters.size();
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      if(cell_index>=init_site_number)
	{
	  cell_perimaters.push_back(-1);
	}else{
	cell_perimaters[cell_index]=-1;
      };
    };
  init_cell_number=cell_displacements.size();
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      if(cell_index*space_dimension>=init_site_number)
	{
	  for(direction_index=0;direction_index<space_dimension;direction_index++)
	    {
	      cell_displacements.push_back(-1);
	    };
	}
      else
	{
	  for(direction_index=0;direction_index<space_dimension;direction_index++)
	    {
	      cell_displacements[cell_index*space_dimension+direction_index]=-1;
	    };
	};
    };
  int space_dimension=model.get_space_dimension();
  init_cell_number=cell_origins.size();
  int component_index;
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      for(component_index=0;component_index<space_dimension;component_index++)
	if(cell_index*((long int)space_dimension+component_index)>=init_site_number)
	{
	  cell_origins.push_back(-1);
	}else{
	cell_origins[cell_index*space_dimension+component_index]=-1;
	};
    };
  //
  init_cell_number=cell_polarities.size();
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      for(component_index=0;component_index<space_dimension;component_index++)
	if(cell_index*((long int)space_dimension+component_index)>=init_site_number)
	{
	  cell_polarities.push_back(-1);
	}else{
	cell_polarities[cell_index*space_dimension+component_index]=-1;
	};
    };
  init_cell_number=cell_total_positions.size();
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      for(component_index=0;component_index<space_dimension;component_index++)
	if(cell_index*((long int)space_dimension+component_index)>=init_site_number)
	  {
	    cell_total_positions.push_back(-10000000);
	  }
	else
	  {
	    cell_total_positions[cell_index*space_dimension+component_index]=-1;
	  };
    };
  //
  for(site_index=0;site_index<number_of_sites;site_index++)
    {
      polar_product.push_back(-100000000.0);
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  external_field.push_back(0.0);
	}
      finite_field_flag.push_back(0);
    };
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      cell_types.push_back(cell_system.get_type(cell_index));
    };
  cell_types.push_back(cell_type_system.get_buffer_type());
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      work_relative_coordinate.push_back(-10000000.0);
      work_cell_polarity.push_back(-10000000.0);
      work_long_int.push_back(0);
      work_external_field.push_back(0.0);
    };
  //
  for(type_index=0;type_index<number_of_cell_types;type_index++)
    {
      work_vector_coupling.push_back(0.0);
    };
  //
  for(
      type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      cell_type_system.get_adhesion_couplings(
					      type_index,
					      "isotropic",
					      work_vector_coupling
					      );
      isotropic_adhesion_couplings.push_back(
					     work_vector_coupling
					     );
      // dipolar
      cell_type_system.get_adhesion_couplings(
					      type_index,
					      "dipolar",
					      work_vector_coupling
					      );
      dipolar_adhesion_couplings.push_back(
					   work_vector_coupling
					   );
      //
      cell_type_system.get_adhesion_basal(
					  type_index,
					  "dipolar",
					  work_double
					  );
      dipolar_adhesion_basals.push_back(
					work_double
					);
      // quadrapolar
      cell_type_system.get_adhesion_couplings(
					      type_index,
					      "quadrapolar",
					      work_vector_coupling
					      );
      quadrapolar_adhesion_couplings.push_back(
					       work_vector_coupling
					       );
      //
      cell_type_system.get_adhesion_basal(
					  type_index,
					  "quadrapolar",
					  work_double
					  );
      quadrapolar_adhesion_basals.push_back(
					    work_double
					    );
    };
  //
  configuration_output_type.clear();
  polarity_output_type.clear();
  input_configuration_output_setting();
  // initialization default parameters
  default_type_for_configuration="random_checkerboard";
  default_type_for_polarity="random";
  //state_finilization_step=std::numeric_limits<long>::max();
  state_finilization_step=9999;
  // work memory for plot
  plot_memory_position.clear();
  plot_displacements.clear();
  //
  for(
      cell_index=0;
      cell_index<number_of_cells;
      cell_index++
      )
    {
      for(
	  direction_index=0;
	  direction_index<space_dimension;
	  direction_index++
	  )
	{
	  adhesion_field.push_back(0.0);
	  plot_memory_position.push_back(0.0);
	  plot_displacements.push_back(0.0);
	};
    };
  for(
      direction_index=0;
      direction_index<space_dimension;
      direction_index++
      )
    {
      work_adhesion_field.push_back(0.0);
    };
  //
};

  /*======================
    Methods for local_state_class
    =======================*/
void local_state_class::initialize_site(
					const long long int & site_index,
					const std::vector<long int> & configuration,
					const std::vector<double> & cell_polarity,
					const std::vector<double> & field,
					const std::vector<int> & finite_field_flag,
					const model_parameters_cellular_potts_class & model,
					const cell_system_class & cell_system,
					const type_system_class & cell_type_system,
					const site_system_class & site_system
					)
{ 
  int neighbor_index;
  int direction_index;
  local_site=site_index;
  site_system.get_neighbor_sites(
				 site_index,
				 neighbor_sites
				 );
  // get local field
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      external_field[direction_index]=field[space_dimension*local_site+direction_index];
    };
  field_flag=finite_field_flag[local_site];
  // get neighbor site data
  for(neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
    {
      //
      neighbor_cells[neighbor_index]=configuration[neighbor_sites[neighbor_index]];
      neighbor_types[neighbor_index]=cell_system.get_type(neighbor_cells[neighbor_index]);
      //
      if(neighbor_cells[neighbor_index]!=cell_system.get_buffer_cell())
	{
	  for(direction_index=0;direction_index<space_dimension;direction_index++)
	    {
	      work_vector_double[direction_index]=cell_polarity[neighbor_cells[neighbor_index]*space_dimension+direction_index];
	    };
	  neighbor_polarities[neighbor_index]=work_vector_double;
	}
      else
	{
	  for(direction_index=0;direction_index<space_dimension;direction_index++)
	    {
	      work_vector_double[direction_index]=0.0;
	    };
	  neighbor_polarities[neighbor_index]=work_vector_double;
	};
    };

};
//
void local_state_class::initialize_local_state(
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
					       )
{
  /*//////////////////////////
  // for point-site local data
  *///////////////////////////
  int direction_index;
  cell=cell_index;
  //  buffer_cell=cell_system.get_buffer_cell();
  cell_type=cell_system.get_type(cell);
  //  fprintf(stderr,"Check: %ld, %lld in \n",cell,local_site);
  if(cell==buffer_cell)
    {
      // if the cell is a buffer, we neglect the procedure here. 
    }
  else if(configuration[local_site]==cell)
    {
    	//debug
  //fprintf(stderr, "debug2:%lld¥n",original_cell_volumes[cell]);
      /*//::::::::::::::::::::::::::::::::::
      // the case of present configuration
      *///::::::::::::::::::::::::::::::::::
      // set coordinate of the pointed site
      site_system.origin_shift(
			       model,
			       cell_origins[cell],
			       local_site,
			       point_vector
			       );
      // set volume of cell
      volume=cell_volumes[cell];
      // set total coordinate for calculation of center of volume 
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  total_coordinates[direction_index]
	    =cell_total_positions[cell*space_dimension+direction_index];
	};
  //
    } 
  else
    {
      /*//::::::::::::::::::::::::::::::::::
      // the case of candidate configuration
      *///::::::::::::::::::::::::::::::::::
      // set coordinate of the pointed site
      site_system.origin_shift(
			       model,
			       cell_origins[cell],
			       local_site,
			       point_vector
			       );
      // set volume of cell
      volume=cell_volumes[cell]+1;
      // set total coordinate for calculation of center of volume 
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  total_coordinates[direction_index]
	    =cell_total_positions[cell*space_dimension+direction_index];
// deleted due to adiabatic approx.
//	    +point_vector[direction_index];
	};
    };
  //
  /*//////////////////////////
  // for polarity & external field
  *///////////////////////////
  if(cell==buffer_cell)
    {
      // if the cell is a buffer, we neglect the procedure here. 
      // because buffer has no polarity and therefore nessecity for calcluation of relative position from center of mass.
    }
  else
    {
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  relative_coordinates[direction_index]
	    =(double)point_vector[direction_index]
	    -(double)total_coordinates[direction_index]
	    /(double)original_cell_volumes[cell];
	};
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  polarities[direction_index]=cell_polarities[space_dimension*cell+direction_index];
	};
    };
  // for neighbor
  int neighbor_index;
  //  long long int neighbor_volume;
  //
  /*//////////////////////////
  // for non-local data needed for calculation cell-cell intraction
  *///////////////////////////
  for(neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
    {
      //
      if(neighbor_cells[neighbor_index]==cell_system.get_buffer_cell())
	{
	  // if the neighbor cell is a buffer, we neglect the procedure here. 
	  // because buffer has no polarity and therefore nessecity for calcluation of relative position from center of mass.
	}
      else if(neighbor_cells[neighbor_index]!=cell)
	{
	  //
	  site_system.origin_shift(
				   model,
				   cell_origins[neighbor_cells[neighbor_index]],
				   neighbor_sites[neighbor_index],
				   point_vector
				   );
	  //
	  if(
	     (configuration[local_site]!=cell)
	     &&(neighbor_cells[neighbor_index]==configuration[local_site])
	     &&(configuration[local_site]!=cell_system.get_buffer_cell())
	     )
	    {
	      // if neighbor cell has a same type of present state (cell),
	      // we should consider displacement of cell center of mass.
	      // here I impliment the exceptional precedure.
	      site_system.origin_shift(
				       model,
				       cell_origins[configuration[local_site]],
				       local_site,
				       original_vector
				       );
// deleted due to adiabatic approx
//	      neighbor_volume=cell_volumes[neighbor_cells[neighbor_index]]-1;
//	      neighbor_volume=cell_volumes[neighbor_cells[neighbor_index]];
	      //
	      for(direction_index=0;direction_index<space_dimension;direction_index++)
		{
		  work_vector[direction_index]
		    =cell_total_positions[neighbor_cells[neighbor_index]*space_dimension+direction_index];
// deleted due to adiabatic approx
//		    -original_vector[direction_index];
		};
	      neighbor_total_coordinates[neighbor_index]=work_vector;
	      //
	    }
	  else 
	    {
	      //	      neighbor_volume=cell_volumes[neighbor_cells[neighbor_index]];
	      //
	      for(direction_index=0;direction_index<space_dimension;direction_index++)
		{
		  work_vector[direction_index]
		    =cell_total_positions[neighbor_cells[neighbor_index]*space_dimension+direction_index];
		};
	      neighbor_total_coordinates[neighbor_index]=work_vector;
	      //
	    };
	  //
	  for(direction_index=0;direction_index<space_dimension;direction_index++)
	    {
	      work_vector_double[direction_index]
		=(double)point_vector[direction_index]
		-(double)neighbor_total_coordinates[neighbor_index][direction_index]
// deleted due to adiabatic approx
//		/(double)neighbor_cell_volume;
		/(double)original_cell_volumes[neighbor_cells[neighbor_index]];
	    };
	  neighbor_relative_coordinates[neighbor_index]=work_vector_double;
	}
      else
	{
	  for(direction_index=0;direction_index<space_dimension;direction_index++)
	    {
	      work_vector_double[direction_index]=0.0;
	    };
	  neighbor_relative_coordinates[neighbor_index]=work_vector_double;
	};
    };
};
//
void local_state_class::debug_output_local_state(
						 const std::string & data_identifier
						 ) 
const {
	int index;
	if(data_identifier=="cell")
	{
		fprintf(stderr,"Cell: %ld, Neighbor cells:", cell);
			for(index=0;index<(int)number_of_neighbor_sites;index++)
			{
				fprintf(stderr,"%ld,",neighbor_cells[index]);
			};
		fprintf(stderr,"\n");
    };
	//	toolbox tool;
	if(data_identifier=="projected_polarities")
	{
				fprintf
				(stderr,"%f,",
		             tool.normalized_product(
					 polarities,
					 relative_coordinates
					 	)
				);
			for(index=0;index<(int)space_dimension;index++)
			{
				fprintf(stderr,"%f,%f/",
		             polarities[index],
					 relative_coordinates[index]
					);
			};
		fprintf(stderr,"::");
			for(index=0;index<(int)number_of_neighbor_sites;index++)
			{
				if((neighbor_cells[index]!=buffer_cell)&&(neighbor_cells[index]!=cell))
				{
				fprintf(stderr,"%f,",
		             tool.normalized_product(
					 neighbor_polarities[index],
					 neighbor_relative_coordinates[index]
					 	)
					);
				}
				else if(neighbor_cells[index]==cell)
				{
					fprintf(stderr,"hom,");
				}
				else
				{
					fprintf(stderr,"buf,");
				};
			};
		fprintf(stderr,"\n");
	};
	if(data_identifier=="coordinates")
	{
		int component_index;
		fprintf(stderr,"(");
		for(component_index=0;component_index<space_dimension;component_index++)
		{
			fprintf(stderr,"%f,",relative_coordinates[component_index]);
		};
		fprintf(stderr,")\n");
			for(index=0;index<(int)number_of_neighbor_sites;index++)
			{
			fprintf(stderr,"(%d:%ld/",index,neighbor_cells[index]);
			for(component_index=0;component_index<space_dimension;component_index++)
				{
				fprintf(stderr,"%f,",neighbor_relative_coordinates[index][component_index]);
				};
			fprintf(stderr,")\n");
			};
	};
};
//
void local_state_class::show_local_state(
					 const std::vector<long int> & configuration
					 )
  const{
  io_cellular_potts io_method;
  std::string message;
  message ="==== local state/ ";
  if(configuration[local_site]==cell)
    {
      message+="status: present";
    }
  else
    {
      message+="status: candidate";
    };
  message+=" ====";
  io_method.standard_output(message);
  //
  message ="# on site information";
  io_method.standard_output(message);
  //
  message ="Site/Cel/Vol/Pol/TCod/RCod:";
  message+=io_method.longlongint_to_string(local_site);
  message+="/";
  message+=io_method.longint_to_string(cell);
  message+="/";
  message+=io_method.longint_to_string(volume);
  message+="/";
  message+=io_method.double_array_to_string(polarities,",");
  message+="/";
  message+=io_method.longlongint_array_to_string(total_coordinates,",");
  message+="/";
  message+=io_method.double_array_to_string(relative_coordinates,",");
  io_method.standard_output(message);
  //
  message ="# on neighbor information";
  io_method.standard_output(message);
  //
  int neighbor_index;
  message ="Sits:";
  message+=io_method.longlongint_array_to_string(neighbor_sites,",");
  io_method.standard_output(message);
  message ="Cels:";
  message+=io_method.longint_array_to_string(neighbor_cells,",");
  io_method.standard_output(message);
  message ="Typs:";
  message+=io_method.longint_array_to_string(neighbor_types,",");
  io_method.standard_output(message);
  message ="TCods:";
  for(neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
    {
      message+=io_method.longlongint_array_to_string(neighbor_total_coordinates[neighbor_index],",");
    }
  io_method.standard_output(message);
  message ="Pols:";
  for(neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
    {
      message+=io_method.double_array_to_string(neighbor_polarities[neighbor_index],",");
    }
  io_method.standard_output(message);
  message ="RCods:";
  for(neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
    {
      message+=io_method.double_array_to_string(neighbor_relative_coordinates[neighbor_index],",");
    }
  io_method.standard_output(message);  
};
//
void local_state_class::show_local_site_list(
					     const model_parameters_cellular_potts_class &model,
					     const site_system_class & site_system,
					     const std::vector<long int> & configuration
					     )
  const{
  std::vector<long long int> neighbor_sites_on_site_system(neighbor_sites.size());
  //  neighbor_sites_on_site_system=site_system.get_neighbor_sites(local_site);
  site_system.get_neighbor_sites(
				 local_site,
				 neighbor_sites_on_site_system
				 );
  std::vector<long int> local_coordinates(space_dimension);
  site_system.get_site_coordinates(
				   local_site,
				   local_coordinates
				   );
  std::vector<double> check_vector(space_dimension,0.0);
  io_cellular_potts io_method;
  //  toolbox tool;
  int component_index;
  std::string output_message;
  output_message ="***Local coordinates ";
  output_message+=io_method.longint_array_to_string(
						    local_coordinates,
						    ","
						    );
  output_message+=" : cell=";
  output_message+=io_method.longint_to_string(cell);
  std::vector<double> cell_coordinates(space_dimension);
  recalculate_cell_center(
			  model,
			  site_system,
			  cell,
			  configuration,
			  cell_coordinates
			  );
  output_message+=io_method.double_array_to_string(
						   cell_coordinates,
						   ","
						   );
  output_message+= " : rel/";
  output_message+=io_method.double_array_to_string(
						   relative_coordinates,
						   ","
						   );
  //
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      check_vector[component_index]=local_coordinates[component_index]-cell_coordinates[component_index];
    };
  if(tool.comparison_double_vector(check_vector,relative_coordinates)==tool.get_value("false")&&(cell!=buffer_cell))
    {
      io_method.standard_output(output_message);
    };
  int neighbor_index;
  for(neighbor_index=0;neighbor_index<(int)neighbor_sites_on_site_system.size();neighbor_index++)
    {
      output_message = "Neighbor ";
      output_message+= io_method.int_to_string(neighbor_index);
      output_message+= "(";
      output_message+= io_method.int_to_string(neighbor_cells[neighbor_index]);
      output_message+= ":";
      output_message+= io_method.int_to_string(cell);
      output_message+= ")";
      output_message+= " : cod/";
      site_system.get_site_coordinates(
				       neighbor_sites[neighbor_index],
				       local_coordinates
				       );
      output_message+= io_method.longint_array_to_string(
							 local_coordinates,
							 ","
							 );
      output_message+= " : cell/";
      output_message+= io_method.longint_to_string(neighbor_cells[neighbor_index]);
      recalculate_cell_center(
			      model,
			      site_system,
			      neighbor_cells[neighbor_index],
			      configuration,
			      cell_coordinates
			      );
      output_message+=io_method.double_array_to_string(
						       cell_coordinates,
						       ","
						       );
      output_message+= " : rel/";
      output_message+=io_method.double_array_to_string(
						       neighbor_relative_coordinates[neighbor_index],
						       ","
						       );
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  check_vector[component_index]=local_coordinates[component_index]-cell_coordinates[component_index];
	};
      if(
	 (tool.comparison_double_vector(check_vector,neighbor_relative_coordinates[neighbor_index])==tool.get_value("false"))
	 &&(cell!=buffer_cell)
	 &&(neighbor_cells[neighbor_index]!=buffer_cell)
	 &&(cell!=neighbor_cells[neighbor_index])
	 )
	{
	  io_method.standard_output(output_message);
	};
    };
};
//
void local_state_class::recalculate_cell_center(
						const model_parameters_cellular_potts_class &model,
						const site_system_class & site_system,
						const long int & cell_index,
						const std::vector<long int> & configuration,
						std::vector<double> & center_coordinates
						)
const {
  long long int site_index;
  int component_index;
  long long int number_of_sites=model.get_number_of_sites();
  long long int cell_volume=0;
  std::vector<long int> local_coordinates(space_dimension,0);
  std::vector<long int> summarized_coordinates(space_dimension,0);
  for(site_index=0;site_index<number_of_sites;site_index++)
    {
      if((configuration[site_index]==cell_index)&&(site_index!=local_site))
	{
	  site_system.get_site_coordinates(
					   site_index,
					   local_coordinates
					   );
	  for(component_index=0;component_index<space_dimension;component_index++)
	    {
	      summarized_coordinates[component_index]=summarized_coordinates[component_index]+local_coordinates[component_index];
	    };
	  cell_volume++;
	}else if((cell==cell_index)&&(site_index==local_site))
	{
	  site_system.get_site_coordinates(
					   site_index,
					   local_coordinates
					   );
	  for(component_index=0;component_index<space_dimension;component_index++)
	    {
	      summarized_coordinates[component_index]=summarized_coordinates[component_index]+local_coordinates[component_index];
	    };
	  cell_volume++;
	};
    };
  for(component_index=0;component_index<space_dimension;component_index++)
    {
     center_coordinates[component_index]=(double)summarized_coordinates[component_index]/(double)cell_volume;
    };
};
  /*======================
    Constructor for local_state_class
   =======================*/
local_state_class::local_state_class(
				     const model_parameters_cellular_potts_class & model,
				     const site_system_class & site_system,
				     const cell_system_class & cell_system,
				     const type_system_class & cell_type_system
				     )
{
  false_value=0;
  true_value=1;
  // array dimensions are set
  space_dimension=model.get_space_dimension();
  number_of_neighbor_sites=site_system.get_number_of_neighbor_sites();
  number_of_cells=model.get_number_of_cells();
  number_of_regions=model.get_number_of_regions();
  /* prepare working memory */
  std::vector<double> initialization_vector_double(space_dimension,0.0);
  std::vector<int> initialization_vector_int(space_dimension,0);
  std::vector<long int> initialization_vector_longint(space_dimension,0);
  std::vector<long long int> initialization_vector_longlongint(space_dimension,0);
  //
  /* for on site */
  // all initialization
  local_site=-1;
  cell=-1;
  volume=-1;
  // initialization of vectors;
  int direction_index;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      polarities.push_back(0.0);
      total_coordinates.push_back(0);
      relative_coordinates.push_back(0.0);
      point_vector.push_back(0);
      work_vector.push_back(0);
      work_vector_double.push_back(0.0);
      original_vector.push_back(0);
      external_field.push_back(0.0);
    };
  // for regions
  int region_index;
  for(region_index=0;region_index<number_of_regions;region_index++)
    {
    };
  // for neighbors
  int neighbor_index;
  for(neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
    {
      neighbor_sites.push_back(-1);
      neighbor_cells.push_back(-1);
      neighbor_types.push_back(-1);
      neighbor_total_coordinates.push_back(initialization_vector_longlongint);
      neighbor_polarities.push_back(initialization_vector_double);
      //      isotropic_adhesion_coupling_zeros.push_back(initialization_vector_int);
      //      dipolar_adhesion_coupling_zeros.push_back(initialization_vector_int);
      //      quadrapolar_adhesion_coupling_zeros.push_back(initialization_vector_int);
      neighbor_relative_coordinates.push_back(initialization_vector_double);
      neighbor_shift.push_back(site_system.get_neighbor_coordinates(neighbor_index));
      work_vector_double_for_neighbors.push_back(0.0);
    };
  int type_index;
  //  int sub_index;
  int number_of_cell_types=model.get_number_of_cell_types();
  double work_double;
  natural_volumes.clear();
  balk_moduluses.clear();
  polarity_sensitivity.clear();
  quadrapolarity_sensitivity.clear();
  buffer_type=cell_type_system.get_buffer_type();
  buffer_cell=cell_system.get_buffer_cell();
  //
  for(type_index=0;type_index<number_of_cell_types;type_index++)
    {
      work_vector_coupling.push_back(0.0);
    };
  isotropic_adhesion_couplings.clear();
  dipolar_adhesion_couplings.clear();
  dipolar_adhesion_basals.clear();
  quadrapolar_adhesion_couplings.clear();
  quadrapolar_adhesion_basals.clear();
  for(
      type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      // isotropic
      cell_type_system.get_adhesion_couplings(
					      type_index,
					      "isotropic",
					      work_vector_coupling
					      );
      isotropic_adhesion_couplings.push_back(
					     work_vector_coupling
					     );
      // dipolar
      cell_type_system.get_adhesion_couplings(
					      type_index,
					      "dipolar",
					      work_vector_coupling
					      );
      dipolar_adhesion_couplings.push_back(
					   work_vector_coupling
					   );
      //
      cell_type_system.get_adhesion_basal(
					  type_index,
					  "dipolar",
					  work_double
					  );
      dipolar_adhesion_basals.push_back(
					work_double
					);
      // quadrapolar
      cell_type_system.get_adhesion_couplings(
					      type_index,
					      "quadrapolar",
					      work_vector_coupling
					      );
      quadrapolar_adhesion_couplings.push_back(
					       work_vector_coupling
					       );
      //
      cell_type_system.get_adhesion_basal(
					  type_index,
					  "quadrapolar",
					  work_double
					  );
      quadrapolar_adhesion_basals.push_back(
					    work_double
					    );
      // 
      /*      for(sub_index=0;sub_index<number_of_cell_types;sub_index++)
	{
	  if(
	     tool.finite_abs_check(
				   isotropic_adhesion_couplings[type_index][sub_index]
				   )
	     ==
	     tool.get_value("true")
	     )
	    {
	      isotropic_adhesion_coupling_zeros[type_index][sub_index]
		=tool.get_value("false");
	    }
	  else
	    {
	      isotropic_adhesion_coupling_zeros[type_index][sub_index]
		=tool.get_value("true");
	    };
	  if(
	     tool.finite_abs_check(
				   dipolar_adhesion_couplings[type_index][sub_index]
				   )
	     ==
	     tool.get_value("true")
	     )
	    {
	      dipolar_adhesion_coupling_zeros[type_index][sub_index]
		=tool.get_value("false");
	    }
	  else
	    {
	      dipolar_adhesion_coupling_zeros[type_index][sub_index]
		=tool.get_value("true");
	    };
	  /*
	  if(
	     tool.finite_abs_check(
				   quadrapolar_adhesion_couplings[type_index][sub_index]
				   )
	     ==
	     tool.get_value("true")
	     )
	    {
	      quadrapolar_adhesion_coupling_zeros[type_index][sub_index]
		=tool.get_value("false");
	    }
	  else
	    {
	      quadrapolar_adhesion_coupling_zeros[type_index][sub_index]
		=tool.get_value("true");
	    };
	};
      */
      //
      natural_volumes.push_back(cell_type_system.get_natural_volume(type_index));
      balk_moduluses.push_back(cell_type_system.get_balk_modulus(type_index));
      //
      polarity_sensitivity.push_back(cell_type_system.get_polarity_sensitivity(type_index));
      quadrapolarity_sensitivity.push_back(cell_type_system.get_quadrapolarity_sensitivity(type_index));
      //
      field_sensitivity.push_back(cell_type_system.get_field_sensitivity(type_index));
    };
};
//
int state_system_class::Metropholis_check_update(
						 const double energy_difference,
						 const double random_number
						 )
  const {
  double transition_probability;
  int return_value=fail_return;
  transition_probability=exp(-beta*energy_difference);
  //  	fprintf(stderr,"transprob,random_number, %f,%f,%f,%f\n",transition_probability,energy_difference,beta,random_number);
  if(transition_probability>random_number)
    {
      return_value=success_return;
    };
  return return_value;
};
//
void state_system_class::update(
				const model_parameters_cellular_potts_class &model,
				const site_system_class & site_system,
				const long long int & candidate_site,
				const long int & candidate_cell
				)
{
  //  int direction_index;
  int space_dimension=model.get_space_dimension();
  std::vector<long int> work_vector(space_dimension);
  //
  if(configuration[candidate_site]>=0&&configuration[candidate_site]<number_of_cells)
  cell_volumes[configuration[candidate_site]]--;
  if(candidate_cell>=0&&candidate_cell<number_of_cells)
  cell_volumes[candidate_cell]++;
  //
  configuration[candidate_site]=candidate_cell;
  //
  add_work();
  //
};
//
void state_system_class::configuration_light_read(
						  const model_parameters_cellular_potts_class & model,
						  site_system_class & site_system
						  )
{
  boost::property_tree::ptree instance_property_tree;
  const std::string filename="configuration.txt";
  boost::property_tree::read_xml(filename.c_str(),instance_property_tree);
  std::string data_identifier;
  // std::string message;
  std::vector<long int> input_coordinates(space_dimension,-1);
  int component_index;
  long long int input_site_index;
  for(long long int site_index=0;site_index<number_of_sites;site_index++)
    {
      data_identifier = "configuration.site[";
      data_identifier+= io_method.longlongint_to_string(site_index);
      data_identifier+= "].coordinates";
      io_method.standard_output(data_identifier);
      //
      component_index=0;
      //
      /*
      site_system.coordinate_generater(
				       model,
				       site_index,
				       input_coordinates
				       );
      */
      //
      BOOST_FOREACH(const boost::property_tree::ptree::value_type& child, instance_property_tree.get_child(data_identifier.c_str()))
	{
	  input_coordinates[component_index]
	    =boost::lexical_cast<long int>(
					   boost::lexical_cast<std::string>(child.second.data())
					   );
	  component_index++;
	};
      //     if(site_index!=site_system.coordinate_to_site(model,input_coordinates))
      //	{
      input_site_index=site_system.coordinate_to_site(model,input_coordinates);
      /*
	io_method.error_output(
	"state system class",
	"configuration_read",
	"due to inconsistencey between input coordinates and site index."
	);
      */
      //	}
      //      else
      //	{
      //	  input_site_index=site_system.coordinate_to_site(model,input_coordinates);
      // input_site_index=site_index;
      //	};
      //
      data_identifier = "configuration.site[";
      data_identifier+= io_method.longlongint_to_string(site_index);
      data_identifier+= "].cell";
      //
      /*
      message = "Debug site(";
      message+= io_method.longlongint_to_string(site_index);
      message+= ":::";
      message+= io_method.longlongint_to_string(input_site_index);
      message+= ") at (";
      message+= io_method.longint_to_string(input_coordinates[0]);
      message+= ",";
      message+= io_method.longint_to_string(input_coordinates[1]);
      message+= ")/(";
      site_system.get_site_coordinates(
				       input_site_index,
				       input_coordinates
				       );
      message+= io_method.longint_to_string(input_coordinates[0]);
      message+= ",";
      message+= io_method.longint_to_string(input_coordinates[1]);
      message+= ") cell=";
      */
      if(boost::optional<std::string> instance_data 
	 = instance_property_tree.get_optional<std::string>(data_identifier.c_str())) 
	{
	  configuration[input_site_index]=boost::lexical_cast<long int>(instance_data.get());
	  //  message+= io_method.longint_to_string(configuration[input_site_index]);
	};
      //
      //io_method.standard_output(message);
    };
};
//
void state_system_class::polarity_light_read(
					     const model_parameters_cellular_potts_class & model,
					     site_system_class & site_system
					     )
{
  boost::property_tree::ptree instance_property_tree;
  const std::string filename="polarity.txt";
  boost::property_tree::read_xml(filename.c_str(),instance_property_tree);
  std::string data_identifier;
  int component_index;
  std::string message;
  for(long int cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      data_identifier = "polarity.cell[";
      data_identifier+= io_method.longlongint_to_string(cell_index);
      data_identifier+= "].components";
      //
      component_index=0;
      //
      BOOST_FOREACH(const boost::property_tree::ptree::value_type& child, instance_property_tree.get_child(data_identifier.c_str()))
	{
	  cell_polarities[cell_index*space_dimension+component_index]
	    =boost::lexical_cast<double>(
					 boost::lexical_cast<std::string>(child.second.data())
					 );
	  component_index++;
	};
      /*
      message = "Debug: Cell No.";
      message+= io_method.longint_to_string(cell_index);
      message+= " : (";
      message+= io_method.double_to_string(cell_polarities[cell_index*space_dimension+0]);
      message+= ",";
      message+= io_method.double_to_string(cell_polarities[cell_index*space_dimension+1]);
      message+= ")";
      io_method.standard_output(message);
      */
      //
    };
};
