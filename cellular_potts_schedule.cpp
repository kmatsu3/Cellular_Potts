#include "cellular_potts_schedule.hpp"
model_parameter_schedule_class::model_parameter_schedule_class()
{
};
//
void model_parameter_schedule_class::set_class_identifier(
							  const std::string & identifier
							  )
{
  class_identifier=identifier;
};
//
std::string model_parameter_schedule_class::get_class_identifier() 
  const {
  return class_identifier;
};
//
//
void model_parameter_schedule_class::set_parameter_identifier(
							      const std::string & identifier
							      )
{
  parameter_identifier=identifier;
};
//
std::string model_parameter_schedule_class::get_parameter_identifier() 
  const {
  return parameter_identifier;
};
//
void model_parameter_schedule_class::set_type_identifier(
						const long int & identifier
						)
{
  type_identifier=identifier;
};
//
long int model_parameter_schedule_class::get_type_identifier() 
  const {
  return type_identifier;
};
//
void model_parameter_schedule_class::set_component_identifier(
							      const long int & identifier
							      )
{
  component_identifier=identifier;
};
//
long int model_parameter_schedule_class::get_component_identifier() 
  const {
  return component_identifier;
};
//
void model_parameter_schedule_class::set_schedule_value(
					       const std::string & identifier,
					       const double & value
					       )
{
  if(identifier=="start_value")
    {
      start_value=value;
    } 
  else if (identifier=="incriment")
    {
      incriment=value;
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "model_parameter_schedule_class",
			     "set_schedule_value",
			     "due to undifined indentifier"
			     );
    };
};
//
double model_parameter_schedule_class::get_schedule_value(
							  const std::string & identifier
							  ) const
{
  double return_value=0.0;
  if(identifier=="start_value")
    {
      return_value=start_value;
    } 
  else if (identifier=="incriment")
    {
      return_value=incriment;
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "model_parameter_schedule_class",
			     "set_schedule_value",
			     "due to undifined indentifier"
			     );
    };
  //
  return return_value;
}
//
void model_parameter_schedule_class::output_schedule()
const {
  io_cellular_potts io_method;
  std::string message;
  message = "Class:" + class_identifier;
  io_method.standard_output(message);
  message = "Type:" + io_method.longint_to_string(type_identifier);
  io_method.standard_output(message);
  message = "Parameter:" + parameter_identifier;
  io_method.standard_output(message);
  message = "Parameter Component:" + io_method.int_to_string(component_identifier);
  io_method.standard_output(message);
  message = "Start value:" + io_method.double_to_string(start_value);
  io_method.standard_output(message);
  message = "Incriment:" + io_method.double_to_string(incriment);
  io_method.standard_output(message);
  message = "------------------------";
  io_method.standard_output(message);
} ;
//
void schedule_system_class::monte_carlo(
					const model_parameters_cellular_potts_class & model,
					cell_system_class & cell_system,
					type_system_class & cell_type_system,
					site_system_class & site_system,
					region_system_class & region_system,
					simulation_system_class & simulation,
					state_system_class & state,
					cell_displacement_system & cell_track,
					observation_system_class & macro_data,
					observables_type_system_class & observables,
					shape_system_class & shape_system,
					adhesion_system_class & adhesion_system
					) 
  {
  long long int monte_carlo_steps;
  int sweep_step;
  clock_t start_time, end_time;
  std::string message;
  //  clock_t sub_start_time, sub_end_time;
  // Initialization
  state.initialize_state(
			 model,
			 cell_system,
			 cell_type_system,
			 site_system,
			 simulation
			 );
  //
  for(
      sweep_step=0;
      sweep_step<number_of_sweep_steps;
      sweep_step++
      )
    {
      //
      observables.initialize(
			     number_of_observations,
			     sweep_step,
			     model,
			     cell_system,	
			     cell_type_system
			     );
      shape_system.initialize(
			      model,
			      cell_system,
			      cell_type_system,
			      site_system,
			      sweep_step
			      );
      //
      input_parameters(
		       cell_type_system,
		       region_system,
		       adhesion_system,
		       sweep_step
		       );
      //
      state.set_model_parameter(
				beta,
				"beta"
				);
      state.set_adhesion_coupling(cell_type_system);
      state.set_external_field(site_system,region_system);
      // Main Loop
      start_time=clock();
      message = "sweep no." + io_method.longlongint_to_string(sweep_step) + "start ////////////\n";
      io_method.standard_output(message);
      fprintf(stderr,"sweep no.%d start ////////////\n",sweep_step);
      for(
	  monte_carlo_steps=0;
	  monte_carlo_steps<monte_carlo_steps_for_relaxzation;
	  monte_carlo_steps++
	  )
	{
	  state.advance_time(
			     model,
			     cell_system,
			     cell_type_system,
			     adhesion_system,
			     site_system,
			     region_system,
			     simulation
			     );
	};
      //
      for(
	  monte_carlo_steps=0;
	  monte_carlo_steps<monte_carlo_steps_for_observation;
	  monte_carlo_steps++
	  )
	{
	  state.advance_time(
			     model,
			     cell_system,
			     cell_type_system,
			     adhesion_system,
			     site_system,
			     region_system,
			     simulation
			     );
	  //
	  cell_track.update(state);
	  observables.track(state);
	  //
	  if(period_of_midplots<=1)
	    {
	      state.plot_mid_state(
				   model,
				   cell_system,
				   cell_type_system,
				   site_system,
				   (long int)monte_carlo_steps,
				   sweep_step
				   );
	    }
	  else if(monte_carlo_steps%period_of_midplots==period_of_midplots-1)
	    {
	      state.plot_mid_state(
				   model,
				   cell_system,
				   cell_type_system,
				   site_system,
				   (long int)monte_carlo_steps/period_of_midplots,
				   sweep_step
				   );
	    };
	  if(period_of_observation<=1)
	    {
	      macro_data.observation(
				     model,
				     cell_system,
				     cell_type_system,
				     site_system,
				     state,
				     monte_carlo_steps
				     );
	      //
	      //cell_track.push_average_mean_square_displacement();
	      cell_track.push();
	      //    observables.track_push();
	    }
	  else if(monte_carlo_steps%period_of_observation==period_of_observation-1)
	    {
	      macro_data.observation(
				     model,
				     cell_system,
				     cell_type_system,
				     site_system,
				     state,
				     monte_carlo_steps/period_of_observation
				     );
	      //
	      observables.calculation(state);
	      shape_system.calculate_polar_correlation(
						       state,
						       cell_system
						       );
	      shape_system.calculate_shape_tensor(
						  state,
						  model,
						  site_system,
						  cell_system
						  );
	      //
	      //	      cell_track.push_average_mean_square_displacement();
	      cell_track.push();
	    };
	};
      //data output
      macro_data.output_observation(parameter_titles,model_parameters);
      macro_data.output_average_observation(parameter_titles,model_parameters);
      //
      observables.output(parameter_titles,model_parameters);
      shape_system.output_polar_correlation(parameter_titles,model_parameters);
      shape_system.output_shape_characters(parameter_titles,model_parameters);
      //cell_track.output_average_mean_square_displacement();
      cell_track.output();
      //
      end_time=clock();
      message = "sweep no." + io_method.longlongint_to_string(sweep_step)
	+ "finished: mcs/time[sec]: Relax:"
	+ io_method.longlongint_to_string(monte_carlo_steps_for_relaxzation)
	+ "/Observe:"
	+ io_method.longlongint_to_string(monte_carlo_steps) + "/"
	+ io_method.longlongint_to_string(period_of_observation)
	+ "mcs/"
	+ io_method.double_to_string((double)(end_time-start_time)/CLOCKS_PER_SEC)
	+ "\n";
      io_method.standard_output(message);
      fprintf(stderr,"sweep no.%d finished: mcs/time[sec]: Relax:%lld-Observe:%lld mcs/%f\n",
	      sweep_step,monte_carlo_steps_for_relaxzation,
	      monte_carlo_steps,(double)(end_time-start_time)/CLOCKS_PER_SEC);
    };
};
//

//
void schedule_system_class::input_parameters(
					     type_system_class & cell_type_system,
					     region_system_class & region_system,
					     adhesion_system_class & adhesion_system,
					     const int & sweep_step
					     )
  {
  //
  int parameter_index;
  int type_index;
  int counter=0;
  int component_index;
  double work_double;
  std::string message;
  for (
       parameter_index=0;
       parameter_index<number_of_control_parameters;
       parameter_index++
       )
    {
      if(parameter_schedule[parameter_index].get_class_identifier()=="simulation" )
	{
	  message = "Temperature (Beta) change from" + io_method.double_to_string(temperature)
	    +"("+ io_method.double_to_string(beta) + ")";
	  if(parameter_schedule[parameter_index].get_parameter_identifier()=="temperature")
	    {
	      temperature = incrimentation(parameter_index,sweep_step);
	      beta=1.0/temperature;
	    }
	  if(parameter_schedule[parameter_index].get_parameter_identifier()=="beta")
	    {
	      beta = incrimentation(parameter_index,sweep_step);
	    }
	  message += " to "+ io_method.double_to_string(temperature)
	    +"("+ io_method.double_to_string(beta) + ")";
	  io_method.standard_output(message);
	  //
	  model_parameters[counter]=beta;
	  counter++;
	}
      else if(parameter_schedule[parameter_index].get_class_identifier()=="adhesion" )
	{
	  type_index=parameter_schedule[parameter_index].get_type_identifier();
	  component_index=parameter_schedule[parameter_index].get_component_identifier();
	  message = "adhesion [" + io_method.int_to_string(type_index) + "] is changed as follows";
	  io_method.standard_output(message);
	  message = "Before:" ;
	  io_method.standard_output(message);
	  adhesion_system.show_adhesion();
	  if(parameter_schedule[parameter_index].get_parameter_identifier()=="coupling_constant")
	    {
	      adhesion_system.get_coupling_constant(type_index);
	      work_double
		= incrimentation(parameter_index,sweep_step);
	      //	      io_method.standard_output("Debug:"+io_method.double_to_string(work_double)+","+io_method.int_to_string(type_index)+","+io_method.int_to_string(sweep_step));
	      adhesion_system.set_coupling_constant(type_index,
						    work_double);
	    }
	  adhesion_system.make_typepair_to_adhesion_maps();
	  //
	  message = "After:" ;
	  io_method.standard_output(message);
	   adhesion_system.show_adhesion();
	  //
	  model_parameters[counter]=work_double;
	  counter++;
	}
      else if(parameter_schedule[parameter_index].get_class_identifier()=="region type" )
	{
	  type_index=parameter_schedule[parameter_index].get_type_identifier();
	  component_index=parameter_schedule[parameter_index].get_component_identifier();
	  message = "region [" + io_method.int_to_string(type_index) + "] is changed as follows";
	  io_method.standard_output(message);
	  message = "Before:" ;
	  io_method.standard_output(message);
	  region_system.show_region(type_index);
	  if(parameter_schedule[parameter_index].get_parameter_identifier()=="external_field")
	    {
	      region_system.get_region_field(type_index,work_field_vector);
	      work_field_vector[component_index]
		= incrimentation(parameter_index,sweep_step);
	      region_system.set_region_field(type_index,work_field_vector);
	    }
	  //
	  message = "After:" ;
	  io_method.standard_output(message);
	  region_system.show_region(type_index);
	  //
	  model_parameters[counter]=work_field_vector[component_index];
	  counter++;
	}
      else if( parameter_schedule[parameter_index].get_class_identifier()=="cell_type")
	{
	  type_index=parameter_schedule[parameter_index].get_type_identifier();
	  component_index=parameter_schedule[parameter_index].get_component_identifier();
	  message = "cell type [" + io_method.int_to_string(type_index) + "] is changed as follows";
	  io_method.standard_output(message);
	  message = "Before:" ;
	  io_method.standard_output(message);
	  cell_type_system.show_type(type_index);
	  // istropic coupling
	  if(parameter_schedule[parameter_index].get_parameter_identifier()=="istropic_coupling")
	    {
	      cell_type_system.get_adhesion_couplings(
						      type_index,
						      "isotropic",
						      work_type_vector
						      );
	      //
	      work_type_vector[component_index]= incrimentation(parameter_index,sweep_step);
	      //
	      model_parameters[counter]=work_type_vector[component_index];
	      counter++;
	      //
	      cell_type_system.set_parameter_double_vector(
							   type_index,
							   "Isotropic_Adhesion_Coupling_constant",
							   work_type_vector
							   );
	      //
	      cell_type_system.get_adhesion_couplings(
						      component_index,
						      "isotropic",
						      work_type_vector
						      );
	      work_type_vector[type_index]= incrimentation(parameter_index,sweep_step);
	      //
	      cell_type_system.set_parameter_double_vector(
							   component_index,
							   "Isotropic_Adhesion_Coupling_constant",
							   work_type_vector
							   );
	    }
	  else if(parameter_schedule[parameter_index].get_parameter_identifier()=="polar_coupling")
	    {
	      cell_type_system.get_adhesion_couplings(
						      type_index,
						      "dipolar",
						      work_type_vector
						      );
	      //
	      work_type_vector[component_index]= incrimentation(parameter_index,sweep_step);
	      //
	      model_parameters[counter]=work_type_vector[component_index];
	      counter++;
	      //
	      cell_type_system.set_parameter_double_vector(
							   type_index,
							   "Dipolar_Adhesion_Coupling_constant",
							   work_type_vector
							   );
	      //
	      cell_type_system.get_adhesion_couplings(
						      component_index,
						      "dipolar",
						      work_type_vector
						      );
	      work_type_vector[type_index]= incrimentation(parameter_index,sweep_step);
	      //
	      cell_type_system.set_parameter_double_vector(
							   component_index,
							   "Dipolar_Adhesion_Coupling_constant",
							   work_type_vector
							   );
	    }
	  else if(parameter_schedule[parameter_index].get_parameter_identifier()=="quadrapolar_coupling")
	    {
	      cell_type_system.get_adhesion_couplings(
						      type_index,
						      "quadrapolar",
						      work_type_vector
						      );
	      //
	      work_type_vector[component_index]= incrimentation(parameter_index,sweep_step);
	      //
	      model_parameters[counter]=work_type_vector[component_index];
	      counter++;
	      //
	      cell_type_system.set_parameter_double_vector(
							   type_index,
							   "Quadrapolar_Adhesion_Coupling_constant",
							   work_type_vector
							   );
	      //
	      cell_type_system.get_adhesion_couplings(
						      component_index,
						      "quadrapolar",
						      work_type_vector
						      );
	      work_type_vector[type_index]= incrimentation(parameter_index,sweep_step);
	      //
	      cell_type_system.set_parameter_double_vector(
							   component_index,
							   "Quadrapolar_Adhesion_Coupling_constant",
							   work_type_vector
							   );
	    }
	  else if(parameter_schedule[parameter_index].get_parameter_identifier()=="persistence_sensitivity")
	    {
	      work_double = incrimentation(parameter_index,sweep_step);
	      cell_type_system.set_parameter_double(
						    type_index,
						    "Polarity_sensitivity",
						    work_double
						    );
	      model_parameters[counter]=work_double;
	      counter++;
	    }
	  else if(parameter_schedule[parameter_index].get_parameter_identifier()=="persistence_quadrapolar_sensitivity")
	    {
	      work_double = incrimentation(parameter_index,sweep_step);
	      cell_type_system.set_parameter_double(
						    type_index,
						    "Quadrapolarity_sensitivity",
						    work_double
						    );
	      model_parameters[counter]=work_double;
	      counter++;
	    }
	  else if(parameter_schedule[parameter_index].get_parameter_identifier()=="persistent_time")
	    {
	      work_double = incrimentation(parameter_index,sweep_step);
	      cell_type_system.set_parameter_double(
						    type_index,
						    "Persistent_time",
						    work_double
						    );
	      model_parameters[counter]=work_double;
	      counter++;
	    }
	  else if(parameter_schedule[parameter_index].get_parameter_identifier()=="adhesion_sensitivity")
	    {
	      work_double = incrimentation(parameter_index,sweep_step);
	      cell_type_system.set_parameter_double(
						    type_index,
						    "Adhesion_sensitivity",
						    work_double
						    );
	      model_parameters[counter]=work_double;
	      counter++;
	    }
	  else if(parameter_schedule[parameter_index].get_parameter_identifier()=="field_sensitivity")
	    {
	      work_double = incrimentation(parameter_index,sweep_step);
	      cell_type_system.set_parameter_double(
						    type_index,
						    "Field_sensitivity",
						    work_double
						    );
	      model_parameters[counter]=work_double;
	      counter++;
	    };
	      // persistence
	      
	      //
	  io_method.standard_output(message);
	  message = "After:" ;
	  io_method.standard_output(message);
	  cell_type_system.show_type(type_index);
	}
      else
	{
	  io_method.error_output(
				 "schedule_system_class",
				 "input_parameters",
				 "due to undifiend class."
				 );
	};
    };
  //
};
//
void schedule_system_class::make_control_parameter_titles(
							  const type_system_class & cell_type_system
							  )
  {
  //
  int parameter_index;
  int type_index;
  int counter=0;
  int component_index;
  std::string message;
  for (
       parameter_index=0;
       parameter_index<number_of_control_parameters;
       parameter_index++
       )
    {
      if(parameter_schedule[parameter_index].get_class_identifier()=="simulation" )
	{
	  parameter_titles[counter]="simulation.";
	  if(parameter_schedule[parameter_index].get_parameter_identifier()=="temperature")
	    {
	      parameter_titles[counter]+="beta";
	      counter++;
	    }
	  if(parameter_schedule[parameter_index].get_parameter_identifier()=="beta")
	    {
	      parameter_titles[counter]+="beta";
	      counter++;
	    }
	  //
	}
      else if( parameter_schedule[parameter_index].get_class_identifier()=="adhesion")
	{
	  parameter_titles[counter]="adhesion.";
	  if(parameter_schedule[parameter_index].get_parameter_identifier()=="coupling_constant")
	    {
	      parameter_titles[counter]+="coupling_constant";
	      counter++;
	    };
	  //
	}
      else if( parameter_schedule[parameter_index].get_class_identifier()=="cell_type")
	{
	  type_index=parameter_schedule[parameter_index].get_type_identifier();
	  component_index=parameter_schedule[parameter_index].get_component_identifier();
	  parameter_titles[counter]="cell type [" 
	    + io_method.int_to_string(type_index) 
	    + "].("
	    + io_method.int_to_string(component_index) 
	    + ")";
	  parameter_titles[counter]+=
	    parameter_schedule[parameter_index].get_parameter_identifier();
	  counter++;
	}
      else
	{
	  io_method.error_output(
				 "schedule_system_class",
				 "input_parameters",
				 "due to undifiend class."
				 );
	};
    };
  //
};
//
inline double schedule_system_class::incrimentation(const int & parameter_index, const int & sweep_step)
{
  return 
    parameter_schedule[parameter_index].get_schedule_value("start_value")
    +(double)sweep_step
    *parameter_schedule[parameter_index].get_schedule_value("incriment");
};
//
void schedule_system_class::input_mid_plot_schedule()
{
  //
  number_of_midplot_time=io_method.get_input_longint(
						     "model_input",
						     "simulation.number_of_midplots"
						     );
  //
  if(number_of_midplot_time>0)
    {
      period_of_midplots=monte_carlo_steps_for_observation/number_of_midplot_time;
    }
  else
    {
      period_of_midplots=monte_carlo_steps_for_observation;
    };
  //
  if(period_of_midplots<=0)
    {
      period_of_midplots=0;
      /*
      io_method.error_output(
			     "schedule_system_class",
			     "input_mid_plot_schedule",
			     "due to low period for mid plot."
			     );
      */
    };
  //
};
//
void schedule_system_class::set_schedule()
{
  int parameter_index=0;
  model_parameter_schedule_class work_schedule;
  parameter_schedule.clear();
  std::string work_string;
  int work_int;
  double work_double;
  while(parameter_index<number_of_control_parameters)
    {
      // Class identifier
      set_structure_type_item(
			      parameter_schedule.size(),
			      "class"
			      );
      work_string=io_method.get_input_string(
					     "model_input",
					     structure_item
					     );
      work_schedule.set_class_identifier(
					 work_string
					 );
      // Type identifier
      set_structure_type_item(
			      parameter_schedule.size(),
			      "type"
			      );
      work_int=io_method.get_input_int(
				       "model_input" ,
				       structure_item
				       );
      work_schedule.set_type_identifier(
					work_int
					);
      // Parameter identifier
      set_structure_type_item(
			      parameter_schedule.size(),
			      "parameter"
			      );
      work_string=io_method.get_input_string(
					     "model_input",
					     structure_item
					     );
      work_schedule.set_parameter_identifier(
					     work_string
					     );
      // Component identifier
      set_structure_type_item(
			      parameter_schedule.size(),
			      "component"
			      );
      work_int=io_method.get_input_int(
				       "model_input",
				       structure_item
				       );
      work_schedule.set_component_identifier(
					     work_int
					     );
      // Start value
      set_structure_type_item(
			      parameter_schedule.size(),
			      "start_value"
			      );
      work_double=io_method.get_input_double(
					     "model_input",
					     structure_item
					     );
      work_schedule.set_schedule_value(
				       "start_value",
				       work_double
				       );
      // Incriment
      set_structure_type_item(
			      parameter_schedule.size(),
			      "incriment"
			      );
      work_double=io_method.get_input_double(
					     "model_input",
					     structure_item
					     );
      work_schedule.set_schedule_value(
				       "incriment",
				       work_double
				       );
      // Push back to Schedule Class
      parameter_schedule.push_back(work_schedule);
      parameter_schedule[parameter_index].output_schedule();
      parameter_index++;
    };
  // temperature & beta
  temperature=io_method.get_input_double(
					 "model_input",
					 "simulation.temperature"
					 );
  //
  if(temperature>0.0)
    {
      beta=1.0/temperature;
    }
  else
    {
      io_method.error_output(
			     "schedule",
			     "input_parameters",
			     "due to zero temperature."
			     );
    };
  //
  //beta=io_method.get_input_double(
  //				  "model_input",
  //				  "simulation.beta"
  //				  );
  //
};
//
void schedule_system_class::set_structure_type_item(
						    int parameter_index,
						    std::string child
						    )
{
  structure_item=io_method.generate_structure(
					      "simulation.schedule",
					      parameter_index,
					      child
					      );
};
//
void schedule_system_class::output_schedule()
{
  std::string message;
  int parameter_index;
  message = "=== Schedule list ===";
  io_method.standard_output(message);
  message = "Period of sweep:" 
    + io_method.int_to_string(number_of_sweep_steps);
  io_method.standard_output(message);
  for(
      parameter_index=0;
      parameter_index<number_of_control_parameters;
      parameter_index++
      )
    {
      parameter_schedule[parameter_index].output_schedule();
    };
};
//
schedule_system_class::schedule_system_class(
					     const model_parameters_cellular_potts_class & model,
					     const type_system_class & cell_type_system,
					     const simulation_system_class & simulation
					     )
{
  monte_carlo_steps_for_observation=simulation.get_number_of_monte_carlo_steps("observation");
  monte_carlo_steps_for_relaxzation=simulation.get_number_of_monte_carlo_steps("relaxation");
  number_of_sweep_steps=simulation.get_number_of_sweep_steps();
  number_of_observations=simulation.get_number_of_observations();
  space_dimension=model.get_space_dimension();
  input_mid_plot_schedule();
  period_of_observation=simulation.get_period_of_observation();
  number_of_control_parameters=simulation.get_number_of_control_parameters();
  number_of_types=model.get_number_of_cell_types();
  for(int type_index=0;type_index<number_of_types;type_index++)
  {
    work_type_vector.push_back(0.0);
  }
  for(int direction_index=0;direction_index<space_dimension;direction_index++)
  {
    work_field_vector.push_back(0.0);
  }
  set_schedule();
  output_schedule();
  int parameter_index;
  for(parameter_index=0;parameter_index<number_of_control_parameters;parameter_index++)
    {
      model_parameters.push_back(0.0);
      parameter_titles.push_back("none");
    };
  make_control_parameter_titles(
				cell_type_system
				);
};
