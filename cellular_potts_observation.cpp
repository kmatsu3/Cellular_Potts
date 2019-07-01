#include "cellular_potts_observation.hpp"
  /*======================
    Methods
   =======================*/
//
void observables_class::initialize(
				   const model_parameters_cellular_potts_class & model
				   )
{
  number_of_living_cells=-1;
  int direction_index;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      net_polarity[direction_index]=0.0;
      variance_of_polarity[direction_index]=0.0;
    };
};
//
observables_class::observables_class(
				     const model_parameters_cellular_potts_class & model
				     )
{
  number_of_living_cells=-1;
  int direction_index;
  space_dimension=model.get_space_dimension();
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      net_polarity.push_back(0.0);
      variance_of_polarity.push_back(0.0);
      net_displacements.push_back(0.0);
      variance_of_displacements.push_back(0.0);
    };
};
//
void observables_class::set_number_of_living_cells( 
						   const double & input_number_of_living_cells
						   )
{
  number_of_living_cells=input_number_of_living_cells;
};
//
double observables_class::get_number_of_living_cells()
const {
  return number_of_living_cells;
};
//
void observables_class::set_net_polarity( 
					 const std::vector<double> & input_average,
					 const std::vector<double> & input_variance
					  )
{
  int direction_index;
  double absolute_value=0.0;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      net_polarity[direction_index]=input_average[direction_index];
      variance_of_polarity[direction_index]=input_variance[direction_index];
      absolute_value+=input_average[direction_index]*input_average[direction_index];
    };
  absolute_value_of_polarity=pow(absolute_value,0.5);
};
//
void observables_class::get_net_polarity( 
					 std::vector<double> & output_average,
					 std::vector<double> & output_variance,
					 double & output_absolute
					  )
const {
  int direction_index;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      output_average[direction_index]=net_polarity[direction_index];
      output_variance[direction_index]=variance_of_polarity[direction_index];
    };
  output_absolute=absolute_value_of_polarity;
};
//
void observables_class::set_total_adhesion( 
					   const double & input_total_adhesion_energy,
					   const double & input_total_intercell_contact,
					   const double & input_total_intercell_head_contact_for_cell,
					   const double & input_total_intercell_tail_contact_for_cell,
					   const double & input_total_intercell_vartical_contact_for_cell,
					   const double & input_total_intercell_lateral_contact_for_cell,
					   const double & input_total_intercell_ordered_contact_for_cell,
					   const double & input_total_cellmedium_contact,
             const double & input_total_square_polar_product
					  )
{
  total_adhesion_energy_for_cell=input_total_adhesion_energy;
  total_intercell_contact_for_cell=input_total_intercell_contact;
  total_intercell_head_contact_for_cell=input_total_intercell_head_contact_for_cell;
  total_intercell_tail_contact_for_cell=input_total_intercell_tail_contact_for_cell;
  total_intercell_vartical_contact_for_cell=input_total_intercell_vartical_contact_for_cell;
  total_intercell_lateral_contact_for_cell=input_total_intercell_lateral_contact_for_cell;
  total_intercell_ordered_contact_for_cell=input_total_intercell_ordered_contact_for_cell;
  total_cellmedium_contact_for_cell=input_total_cellmedium_contact;
  total_square_polar_product=input_total_square_polar_product;
};
//
void observables_class::set_volume( 
				   const double & input_total_volume_energy,
				   const double & input_volume_fraction
				    )
{
  total_volume_energy_for_cell=input_total_volume_energy;
  volume_fraction=input_volume_fraction;
};
//
double observables_class::get_component( 
					const std::string & data_identifier,
					const int & component_index
					 )
  const {
  double return_value(0.0);
  if(data_identifier=="number_of_living_cells")
    {
      return_value=number_of_living_cells;
    }
  else if(data_identifier=="polarity.average")
    {
      return_value=net_polarity[component_index];
    }
  else if(data_identifier=="polarity.variance")
    {
      return_value=variance_of_polarity[component_index];
    }
  else if(data_identifier=="polarity.absolute")
    {
      return_value=absolute_value_of_polarity;
    }
  else if(data_identifier=="displacement.average")
    {
      return_value=net_displacements[component_index];
    }
  else if(data_identifier=="displacement.variance")
    {
      return_value=variance_of_displacements[component_index];
    }
  else if(data_identifier=="displacement.absolute")
    {
      return_value=absolute_value_of_displacement;
    }
  else if(data_identifier=="cell_squared_displacement")
    {
      return_value=cell_squared_displacement;
    }
  else if(data_identifier=="adhesion.energy")
    {
      return_value=total_adhesion_energy_for_cell;
    }
  else if(data_identifier=="adhesion.contact")
    {
      return_value=total_intercell_contact_for_cell;
    }
  else if(data_identifier=="adhesion.head_contact")
    {
      return_value=total_intercell_head_contact_for_cell;
    }
  else if(data_identifier=="adhesion.tail_contact")
    {
      return_value=total_intercell_tail_contact_for_cell;
    }
  else if(data_identifier=="adhesion.vartical_contact")
    {
      return_value=total_intercell_vartical_contact_for_cell;
    }
  else if(data_identifier=="adhesion.lateral_contact")
    {
      return_value=total_intercell_lateral_contact_for_cell;
    }
  else if(data_identifier=="adhesion.ordered_contact")
    {
      return_value=total_intercell_ordered_contact_for_cell;
    }
  else if(data_identifier=="adhesion.noncontact")
    {
      return_value=total_cellmedium_contact_for_cell;
    }
  else if(data_identifier=="cell_shape.quadrapolarity")
    {
      return_value=total_square_polar_product;
    }
  else if(data_identifier=="volume.energy")
    {
      return_value=total_volume_energy_for_cell;
    }
  else if(data_identifier=="volume.fraction")
    {
      return_value=volume_fraction;
    };
  return return_value;
};
//
void observables_class::set_net_displacement( 
					     const std::vector<double> & input_average,
					     const std::vector<double> & input_variance
					      )
{
  int direction_index;
  double absolute_value=0.0;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      net_displacements[direction_index]=input_average[direction_index];
      variance_of_displacements[direction_index]=input_variance[direction_index];
      absolute_value+=input_average[direction_index]*input_average[direction_index];
    };
  absolute_value_of_displacement=pow(absolute_value,0.5);
};
//
void observables_class::get_net_displacement( 
					     std::vector<double> & output_average,
					     std::vector<double> & output_variance
					      )
const {
  int direction_index;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      output_average[direction_index]=net_displacements[direction_index];
      output_variance[direction_index]=variance_of_displacements[direction_index];
    };
  fprintf(
	  stderr,
	  "(%f,%f)\n",
	  output_variance[0],
	  output_variance[1]
	  );
};
//
void observables_class::set_cell_squared_displacement( 
						      const double & input_value
						      )
{
  cell_squared_displacement=input_value;
};
//
double observables_class::get_cell_squared_displacement()
const {
  return cell_squared_displacement;
};
//
void observation_system_class::set_observable_settings()
{
};
//
void observation_system_class::initialize_load_flags()
{
  load_volume_flag=unload_value;
  load_polarity_flag=unload_value;
  load_displacement_flag=unload_value;
};
//
void observation_system_class::observation(
					   const model_parameters_cellular_potts_class & model,
					   const cell_system_class & cell_system,
					   const type_system_class & cell_type_system,
					   site_system_class & site_system,
					   state_system_class & state,
					   const long long int & monte_carlo_step
					   )
{
  initialize_load_flags();
  //
  load_volumes(state);
  observables[monte_carlo_step].set_number_of_living_cells(
							   calculate_number_of_living_cells(
											    model,
											    state  
											    )
							   );
  //
  load_polarities(state);
  calculate_net_polarity(
			 model,
			 state,
			 work_average
			 );
  //
  squared_vector(
		 cell_polarities,
		 work_producted_vectors,
		 number_of_cells*space_dimension
		 );
  //
  calculate_cell_variance_of_polarity(
				      model,
				      state,
				      work_average,
				      work_variance
				      );
  //
  observables[monte_carlo_step].set_net_polarity(
						 work_average,
						 work_variance
						 );
  //
  load_cell_displacements(state);
  //
  calculate_net_displacement(
			     model,
			     state,
			     work_average
			     );
  //
  squared_vector(
		 cell_displacements,
		 work_producted_vectors,
		 number_of_cells*space_dimension
		 );
  //
  calculate_cell_variance_of_displacement(
					  model,
					  state,
					  work_average,
					  work_variance
					  );
  //
  observables[monte_carlo_step].set_net_displacement(
						     work_average,
						     work_variance
						     );
  //
  double cell_squared_displacement=0.0;
  calculate_cell_average_squared_displacement(
					      model,
					      state,
					      cell_squared_displacement
					      );
  //
  observables[monte_carlo_step].set_cell_squared_displacement(
							      cell_squared_displacement
							      );
  //
  double total_adhesion_energy_for_cell=0.0;
  double total_intercell_contact_for_cell=0.0;
  double total_intercell_head_contact_for_cell=0.0;
  double total_intercell_tail_contact_for_cell=0.0;
  double total_intercell_vartical_contact_for_cell=0.0;
  double total_intercell_lateral_contact_for_cell=0.0;
  double total_intercell_ordered_contact_for_cell=0.0;
  double total_cellmedium_contact_for_cell=0.0;
  double total_square_polar_product=0.0;
  state.get_total_adhesion(
			   model,
			   cell_system,
			   site_system,
			   total_adhesion_energy_for_cell,
			   total_intercell_contact_for_cell,
			   total_intercell_head_contact_for_cell,
			   total_intercell_tail_contact_for_cell,
			   total_intercell_vartical_contact_for_cell,
			   total_intercell_lateral_contact_for_cell,
			   total_intercell_ordered_contact_for_cell,
			   total_cellmedium_contact_for_cell,
          total_square_polar_product
			   );
  //
  observables[monte_carlo_step].set_total_adhesion(
						   total_adhesion_energy_for_cell,
						   total_intercell_contact_for_cell,
						   total_intercell_head_contact_for_cell,
						   total_intercell_tail_contact_for_cell,
						   total_intercell_vartical_contact_for_cell,
						   total_intercell_lateral_contact_for_cell,
						   total_intercell_ordered_contact_for_cell,
						   total_cellmedium_contact_for_cell,
                total_square_polar_product
						   );
  //
  double total_volume_energy=0.0;
  double average_volume_fraction=0.0;
  calculate_volume(
		   cell_type_system,
		   cell_system,
		   total_volume_energy,
		   average_volume_fraction
		   );						   
  //
  observables[monte_carlo_step].set_volume(
					   total_volume_energy,
					   average_volume_fraction
					   );
  //
};
//
void observation_system_class::output_observation(
						  std::vector<std::string> & parameter_titles,
						  std::vector<double> & model_parameters
)
  {
  std::string message;
  long long int observation_index;
  int parameter_index;
  std::vector<double> sorted_observables(number_of_observables,0.0);
  std::vector<std::string> observable_titles(number_of_observables,"none");
  int data_index;
  if(observation_output_format=="simple")
    {
      message="";
      io_method.output_message(message,"observation_result.txt");
      io_method.output_message(message,"observation_result.txt");
      get_sorted_observables(0,sorted_observables,observable_titles,"title");
      message="#";
      for(parameter_index=0;parameter_index<(int)parameter_titles.size();parameter_index++)
	{
	  message+=parameter_titles[parameter_index];
	  if(parameter_index<(int)parameter_titles.size()-1) message+="/";
	};
      message+="::";
      for(data_index=0;data_index<number_of_observables;data_index++)
	{
	  message+=observable_titles[data_index];
	  if(data_index<(int)observable_titles.size()-1) message+="/";
	};
      io_method.output_message(message,"observation_result.txt");
      for(observation_index=0;observation_index<number_of_observations;observation_index++)
	{
	  message="";
	  for(parameter_index=0;parameter_index<(int)model_parameters.size();parameter_index++)
	    {
	      message+=io_method.double_to_string(model_parameters[parameter_index]);
	      message+=" ";
	    };
	  //
	  get_sorted_observables(observation_index,sorted_observables,observable_titles,"untitle");
	  //
	  for(data_index=0;data_index<number_of_observables;data_index++)
	    {
	      message+=io_method.double_to_string(sorted_observables[data_index]);
	      message+=" ";
	    };
	  //
	  io_method.output_message(message,"observation_result.txt");
	}
    }
  else if (observation_output_format=="xml")
    {
      io_method.error_output(
			     "observation_system_class",
			     "output_observations",
			     "due to under construction in output observation format."
			     );
    }
  else
    {
      io_method.error_output(
			     "observation_system_class",
			     "output_observations",
			     "due to undefined output observation format."
			     );
    }
};
//
void observation_system_class::output_average_observation(
							  std::vector<std::string> & parameter_titles,
							  std::vector<double> & model_parameters
							  )
{
  std::string message;
  long long int observation_index;
  int parameter_index;
  std::vector<double> sorted_observables(number_of_observables,0.0);
  std::vector<std::string> observable_titles(number_of_observables,"");
  std::vector<double> average(number_of_observables,0.0);
  std::vector<double> squared_average(number_of_observables,0.0);
  int data_index;
  if(observation_output_format=="simple")
    {
      message="#";
      for(parameter_index=0;parameter_index<(int)parameter_titles.size();parameter_index++)
	{
	  message+=parameter_titles[parameter_index];
	  if(parameter_index<(int)parameter_titles.size()-1) message+="/";
	};
      get_sorted_observables(0,sorted_observables,observable_titles,"title");
      message+="::";
      for(data_index=0;data_index<number_of_observables;data_index++)
	{
	  message+=observable_titles[data_index];
	  if(data_index<(int)observable_titles.size()-1) message+="/";
	};
      io_method.output_message(message,"observation_average.txt");
      for(observation_index=0;observation_index<number_of_observations;observation_index++)
	{
	  get_sorted_observables(observation_index,sorted_observables,observable_titles,"untitle");
	  //
	  for(data_index=0;data_index<number_of_observables;data_index++)
	    {
	      average[data_index]+=sorted_observables[data_index];
	      squared_average[data_index]+=sorted_observables[data_index]*sorted_observables[data_index];
	    };
	  //
	};
      for(data_index=0;data_index<number_of_observables;data_index++)
	{
	  average[data_index]=average[data_index]/(double)number_of_observations;
	  squared_average[data_index]=squared_average[data_index]/(double)number_of_observations;
	};
      message="";
      for(parameter_index=0;parameter_index<(int)model_parameters.size();parameter_index++)
	{
	  message+=io_method.double_to_string(model_parameters[parameter_index]);
	  message+=" ";
	};
      //
      //
      for(data_index=0;data_index<number_of_observables;data_index++)
	{
	  message+=io_method.double_to_string(average[data_index]);
	  message+=" ";
	  message+=io_method.double_to_string(squared_average[data_index]);
	  message+=" ";
	};
      //
      io_method.output_message(message,"observation_average.txt");
    }
  else if (observation_output_format=="xml")
    {
      io_method.error_output(
			     "observation_system_class",
			     "output_observations",
			     "due to under construction in output observation format."
			     );
    }
  else
    {
      io_method.error_output(
			     "observation_system_class",
			     "output_observations",
			     "due to undefined output observation format."
			     );
    }
};
//
void observation_system_class::get_sorted_observables(
						      const long long int & observation_index,
						      std::vector<double> & sorted_observables,
						      std::vector<std::string> & observable_titles,
						      const std::string & job_flag
						      )
const {
  int counter=0;
  int component_index;
  // inpput values
  sorted_observables[counter]
    =observables[observation_index].get_component("number_of_living_cells",0);
  counter=counter+1;
  //
  for(component_index=0;component_index<space_dimension;component_index++)
    {
    sorted_observables[counter+component_index]
      =observables[observation_index].get_component("polarity.average",component_index);
    };
  counter=counter+space_dimension;
  //
  for(component_index=0;component_index<space_dimension;component_index++)
    {
    sorted_observables[counter+component_index]
      =observables[observation_index].get_component("polarity.variance",component_index);
    };
  counter=counter+space_dimension;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("polarity.absolute",0);
  counter=counter+1;
  //
  for(component_index=0;component_index<space_dimension;component_index++)
    {
    sorted_observables[counter+component_index]
      =observables[observation_index].get_component("displacement.average",component_index);
    };
  counter=counter+space_dimension;
  //
  for(component_index=0;component_index<space_dimension;component_index++)
    {
    sorted_observables[counter+component_index]
      =observables[observation_index].get_component("displacement.variance",component_index);
    };
  counter=counter+space_dimension;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("displacement.absolute",0);
  counter=counter+1;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("cell_squared_displacement",0);
  counter=counter+1;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("adhesion.energy",0);
  counter=counter+1;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("adhesion.contact",0);
  counter=counter+1;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("adhesion.head_contact",0);
  counter=counter+1;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("adhesion.tail_contact",0);
  counter=counter+1;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("adhesion.vartical_contact",0);
  counter=counter+1;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("adhesion.lateral_contact",0);
  counter=counter+1;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("adhesion.ordered_contact",0);
  counter=counter+1;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("adhesion.noncontact",0);
  counter=counter+1;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("cell_shape.quadrapolarity",0);
  counter=counter+1;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("volume.energy",0);
  counter=counter+1;
  //
  sorted_observables[counter]
    =observables[observation_index].get_component("volume.fraction",0);
  counter=counter+1;
  // make titles
  if(job_flag=="title")
    {
      counter=0;
      observable_titles[counter]="number of living cells";
      counter=counter+1;
      //
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  observable_titles[counter+component_index]
	    ="polarity.average[" + io_method.int_to_string(component_index)  + "]";
	};
      counter=counter+space_dimension;
  //
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  observable_titles[counter+component_index]
	    ="polarity.variance[" + io_method.int_to_string(component_index)  + "]";
	};
      counter=counter+space_dimension;
  //
      observable_titles[counter]="polarity.absolute";
      counter=counter+1;
  //
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  observable_titles[counter+component_index]
	    ="displacement.average[" + io_method.int_to_string(component_index)  + "]";
	};
      counter=counter+space_dimension;
  //
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  observable_titles[counter+component_index]
	    ="displacement.variance[" + io_method.int_to_string(component_index)  + "]";
	};
      counter=counter+space_dimension;
  //
      observable_titles[counter]="displacement.absolute";
      counter=counter+1;
  //
      observable_titles[counter]="cell_squared_displacement";
      counter=counter+1;
  //
      observable_titles[counter]="adhesion.energy";
      counter=counter+1;
  //
      observable_titles[counter]="adhesion.contact";
      counter=counter+1;
  //
      observable_titles[counter]="adhesion.head_contact";
      counter=counter+1;
  //
      observable_titles[counter]="adhesion.tail_contact";
      counter=counter+1;
  //
      observable_titles[counter]="adhesion.vartical_contact";
      counter=counter+1;
  //
      observable_titles[counter]="adhesion.lateral_contact";
      counter=counter+1;
  //
      observable_titles[counter]="adhesion.ordered_contact";
      counter=counter+1;
  //
      observable_titles[counter]="adhesion.noncontact";
      counter=counter+1;
  //
      observable_titles[counter]="cell_shape.quadrapolarity";
      counter=counter+1;
  //
      observable_titles[counter]="volume.energy";
      counter=counter+1;
  //
      observable_titles[counter]="volume.fraction";
      counter=counter+1;
    };
};
//
void observation_system_class::load_volumes(
					    const state_system_class & state  
					    )
{
  state.get_cell_volumes(cell_volumes);
  load_volume_flag=load_value;
};
//
double observation_system_class::calculate_number_of_living_cells(
								  const model_parameters_cellular_potts_class & model,
								  const state_system_class & state  
								  ) 
 const {
  double return_value=0.0;
  if(load_volume_flag==load_value)
    {
      std::vector<long long int>::const_iterator cell_index=cell_volumes.begin();
      while(cell_index!=cell_volumes.end())
	{
	  if(*cell_index>0)
	    {
	      return_value=return_value+1.0;
	    };
	  cell_index++;
	};
    } 
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "observation_system_class",
			     "calculate_number_of_living_cells",
			     "due to unloaded cell_volumes."
			     );
    };
  return return_value;
};
//
void observation_system_class::load_polarities(
					       const state_system_class & state  
					       )
{
  state.get_sorted_cell_polarities(
				   cell_polarities
				   );
  //
  load_polarity_flag=load_value;
};
//
void observation_system_class::calculate_net_polarity(
						      const model_parameters_cellular_potts_class & model,
						      const state_system_class & state,
						      std::vector<double> & net_polarity
						      ) 
  const {
  //
  int component_index;
  //]
  if(load_polarity_flag==load_value)
    {
      sorted_vector_sum(
			net_polarity,
			cell_polarities,
			number_of_cells,
			space_dimension
			);
  //
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  net_polarity[component_index]
	    =net_polarity[component_index]/((double)number_of_cells);
	};
      //
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "observation_system_class",
			     "calculate_net_polarity",
			     "due to unloaded cell polarities."
			     );
    };
};
//
inline void observation_system_class:: sorted_vector_sum(
							 std::vector<double> & sum_vector,
							 const std::vector<double> & sorted_vectors,
							 const long int & number_of_vectors,
							 const int & vector_size
							 )
const {
  int component_index;
  for(component_index=0;component_index<vector_size;component_index++)
    {
      sum_vector[component_index]
	=std::accumulate(
			 sorted_vectors.begin()+(number_of_vectors* component_index),
			 sorted_vectors.begin()+(number_of_vectors*(component_index+1)),
			 0.0
			 );
    };
};
//
void observation_system_class::calculate_cell_variance_of_polarity(
								   const model_parameters_cellular_potts_class & model,
								   const state_system_class & state,
								   const std::vector<double> & net_polarity,
								   std::vector<double> & polarity_variance
								   ) 
  const{
  int component_index;
  if(load_polarity_flag==load_value)
    {
      sorted_vector_sum(
			polarity_variance,
			work_producted_vectors,
			number_of_cells,
			space_dimension
			);
  //
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  polarity_variance[component_index]
	    =(polarity_variance[component_index]/((double)number_of_cells)
	      -net_polarity[component_index]*net_polarity[component_index]);
	};
      //
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "observation_system_class",
			     "calculate_cell_variance_of_polarity",
			     "due to unloaded cell polarities."
			     );
    };
};
//
inline void observation_system_class:: squared_vector(
						      const std::vector<double> & input_vector,
						      std::vector<double> & output_vector,
						      const long int & vector_size
						      )
{
  long int vector_index;
  for(vector_index=0;vector_index<vector_size;vector_index++)
    {
      output_vector[vector_index]=input_vector[vector_index]*input_vector[vector_index];
    };
};
//
/*
inline void observation_system_class:: sorted_vector_producted_sum(
								   std::vector<double> & sum_vector,
								   const std::vector<double> & sorted_vectors,
								   const long int & number_of_vectors,
								   const int & vector_size
								   )
const {
  int component_index;
  for(component_index=0;component_index<vector_size;component_index++)
    {
      sum_vector[component_index]
	=std::accumulate(
			 sorted_vectors.begin()+(number_of_vectors* component_index),
			 sorted_vectors.begin()+(number_of_vectors*(component_index+1)),
			 0.0
			 );
    };
};
*/
//
void observation_system_class::load_cell_displacements(
						       const state_system_class & state  
						       )
{
  state.get_sorted_cell_displacements(
				      cell_displacements
				      );
  //
  load_displacement_flag=load_value;
};
//
void observation_system_class::calculate_net_displacement(
							  const model_parameters_cellular_potts_class & model,
							  const state_system_class & state,
							  std::vector<double> & net_cell_displacements		       
							  ) 
  const {
  int component_index;
  if(load_displacement_flag==load_value)
    {
      sorted_vector_sum(
			net_cell_displacements,
			cell_displacements,
			number_of_cells,
			space_dimension
			);
      //
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  net_cell_displacements[component_index]
	    =net_cell_displacements[component_index]/((double)number_of_cells);
	};
      //
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "observation_system_class",
			     "calculate_net_displacement",
			     "due to unloaded cell_displacements."
			     );
    };
};
//
void observation_system_class::calculate_cell_variance_of_displacement(
								       const model_parameters_cellular_potts_class & model,
								       const state_system_class & state,
								       const std::vector<double> & net_cell_displacement,
								       std::vector<double> & displacement_variance
								       ) 
  const{
  int component_index;
  if(load_displacement_flag==load_value)
    {
      sorted_vector_sum(
			displacement_variance,
			work_producted_vectors,
			number_of_cells,
			space_dimension
			);
      //
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  displacement_variance[component_index]
	    =(displacement_variance[component_index]/((double)number_of_cells)
	      -net_cell_displacement[component_index]*net_cell_displacement[component_index]);
	};
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "observation_system_class",
			     "calculate_cell_variance_of_displacement",
			     "due to unloaded cell displacements."
			     );
    };
};
//
void observation_system_class::calculate_cell_average_squared_displacement(
									   const model_parameters_cellular_potts_class & model,
									   const state_system_class & state,
									   double & average_squared_displacement
									   ) 
  const{
  //  int component_index;
  if(load_displacement_flag==load_value)
    {
      average_squared_displacement
	=std::inner_product(
			    cell_displacements.begin(),
			    cell_displacements.end(),
			    cell_displacements.begin(),
			    0.0
			    );
      average_squared_displacement=average_squared_displacement/(double)number_of_cells;
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "observation_system_class",
			     "calculate_cell_average_displacement",
			     "due to unloaded cell displacements."
			     );
    };
};
//
double observation_system_class::calculate_adhesion_hamiltonian(
								const cell_system_class & cell_system,
								const type_system_class & cell_type_system,
								const state_system_class & state  
								) 
  const{
  double return_value;
  fprintf(stderr,"Not implementd (calculate_cell_variance_of_displacement)");
  abort();
  
  return return_value;
};
//
void observation_system_class::calculate_volume(
						const type_system_class & cell_type_system,
						const cell_system_class & cell_system,
						double & volume_energy,  
						double & volume_fraction 
						) 
{
  long int cell_index;
  int type_index;
  for(type_index=0;type_index<number_of_cell_types;type_index++)
    {
      balk_modulus[type_index]=cell_type_system.get_balk_modulus(type_index);
      natural_volumes[type_index]=cell_type_system.get_natural_volume(type_index);
    };
  volume_energy=0.0;
  volume_fraction=0.0;
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      type_index=cell_system.get_type(cell_index);
      if(type_index!=buffer_type&&cell_volumes[cell_index]>0)
	{
	  volume_energy+=balk_modulus[type_index]
	    *(double)
	    (
	     (cell_volumes[cell_index]-natural_volumes[type_index])
	     *(cell_volumes[cell_index]-natural_volumes[type_index])
	     );
	  volume_fraction+=(double)cell_volumes[cell_index]/(double)natural_volumes[type_index];
	};
    };
  volume_energy=volume_energy/(double)number_of_cells;
  volume_fraction=volume_fraction/(double)number_of_cells;
};
//
double observation_system_class::calculate_contact_number(
							  const std::string anisotropy
							  ) 
  const{
  double return_value;
  fprintf(stderr,"Not implementd (calculate_cell_variance_of_displacement)");
  abort();
  return return_value;
};
  /*======================
    Constructor
   =======================*/
observation_system_class::observation_system_class(
						   const model_parameters_cellular_potts_class & model,
						   const cell_system_class & cell_system,
						   const type_system_class & cell_type_system,
						   const simulation_system_class & simulation
						   )
{
  std::string file_identifier="observation_result";
  buffer_cell=cell_system.get_buffer_cell();
  buffer_type=cell_type_system.get_buffer_type();
  io_method.file_initialize(file_identifier);
  number_of_control_parameters=simulation.get_number_of_control_parameters();
  space_dimension=model.get_space_dimension();
  number_of_cells=model.get_number_of_cells();
  number_of_cell_types=model.get_number_of_cell_types();
  number_of_observables=1+space_dimension*4+2+2+1+1+4+2+1+1;
  observation_output_format=simulation.get_observation_output_format();
  long int cell_index;
  int direction_index;
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      cell_volumes.push_back(0);
      //
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  cell_polarities.push_back(0.0);
	  cell_displacements.push_back(0.0);
	  work_producted_vectors.push_back(0.0);
	};
      //
    };
  //
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      work_average.push_back(0.0);
      work_variance.push_back(0.0);
    };
  int type_index;
  for(type_index=0;type_index<number_of_cell_types;type_index++)
    {
      natural_volumes.push_back(0);
      balk_modulus.push_back(0.0);
    };
  //
  load_value=1;
  unload_value=0;
  steps_for_observation=simulation.get_number_of_monte_carlo_steps("observation");
  number_of_observations=simulation.get_number_of_observations();
  period_of_observation=simulation.get_period_of_observation();
  observables_class initialization_observables(model);
  initialization_observables.initialize(model);
  long long int observation_step;
  for(observation_step=0;observation_step<number_of_observations;observation_step++)
    {
      observables.push_back(initialization_observables);
    };
  io_method.file_initialize("observation_result");
  io_method.file_initialize("observation_average");
};
