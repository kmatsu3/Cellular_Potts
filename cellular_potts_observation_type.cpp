

#include "cellular_potts_observation_type.hpp"
observables_type_class::observables_type_class()
{
  net_polarity.clear();
  variance_of_polarity.clear();
  net_displacement.clear();
  variance_of_displacement.clear();
  total_displacement.clear();
  //
  variance_longitudinal_polarity.clear();
  variance_lateral_polarity.clear();
  longitudinal_displacement.clear();
  lateral_displacement.clear(); 
  variance_longitudinal_displacement.clear();
  variance_lateral_displacement.clear(); 
};
//
observables_type_system_class::observables_type_system_class()
{
};
//
void observables_type_system_class::initialize(
					       const long long int & input_number_of_observations,
					       const int & input_sweep_step,
					       const model_parameters_cellular_potts_class & model,
					       const cell_system_class & cell_system,	
					       const type_system_class & cell_type_system
					       )
{
  std::string message;
  number_of_cells=model.get_number_of_cells();
  number_of_cell_types=model.get_number_of_cell_types();
  number_of_observations=input_number_of_observations;
  space_dimension=model.get_space_dimension();
  observables_type_class work_observables_type;
  long int type_pointer=0;
  sweep_step=input_sweep_step;
  iteration_counter=0;
  matrix_dimension=space_dimension*(space_dimension+1)/2;
  long int number_of_cells_for_type=0;
  track_steps=0;
  track_index=0;
  observables_types.clear();
  for( 
      int type_index = 0;
      type_index < number_of_cell_types;
      type_index++
       )
    {
      number_of_cells_for_type
	=cell_type_system.get_number_of_cells(type_index);
      work_observables_type.type=type_index;
      work_observables_type.space_dimension=space_dimension;
      work_observables_type.start_pointer=type_pointer;
      if(number_of_cells_for_type>0)
	{
	  type_pointer+=number_of_cells_for_type;
	};
      work_observables_type.end_pointer=type_pointer;
      work_observables_type.number_of_living_cells.clear();
      work_observables_type.net_polarity.clear();
      work_observables_type.variance_of_polarity.clear();
      work_observables_type.variance_longitudinal_polarity.clear();
      work_observables_type.variance_lateral_polarity.clear();
      work_observables_type.net_displacement.clear();
      work_observables_type.variance_of_displacement.clear();
      work_observables_type.longitudinal_displacement.clear();
      work_observables_type.lateral_displacement.clear(); 
      work_observables_type.variance_longitudinal_displacement.clear();
      work_observables_type.variance_lateral_displacement.clear(); 
      work_observables_type.total_displacement.clear();
      work_observables_type.represent_position.clear();
      work_observables_type.averaged_position.clear();
      for(
	  long int cell_index=0;
	  cell_index<number_of_cells_for_type;
	  cell_index++
	  )
	{
	  for(
	      int direction_index=0;
	      direction_index<space_dimension;
	      direction_index++
	      )
	    {
	      work_observables_type.total_displacement.push_back(0.0);
	    };
	};
      work_observables_type.total_displacement.push_back(0.0);
      work_observables_type.represent_position.push_back(0.0);
      work_observables_type.represent_position.push_back(0.0);
      work_observables_type.averaged_position.push_back(0.0);
      work_observables_type.averaged_position.push_back(0.0);
      for(
	  long long int iteration_index=0;
	  iteration_index<number_of_observations;
	  iteration_index++
	  )
	{
	  work_observables_type.number_of_living_cells.push_back(0);
	  work_observables_type.total_displacement.push_back(0.0);
	  work_observables_type.longitudinal_displacement.push_back(0.0);
	  work_observables_type.lateral_displacement.push_back(0.0);
	  work_observables_type.variance_longitudinal_polarity.push_back(0.0);
	  work_observables_type.variance_lateral_polarity.push_back(0.0);
	  work_observables_type.variance_longitudinal_displacement.push_back(0.0);
	  work_observables_type.variance_lateral_displacement.push_back(0.0);
	  for(
	      int direction_index=0;
	      direction_index<space_dimension;
	      direction_index++
	      )
	    {
	      work_observables_type.net_polarity.push_back(0.0);
	      work_observables_type.net_displacement.push_back(0.0);
	      work_observables_type.represent_position.push_back(0.0);
	      work_observables_type.averaged_position.push_back(0.0);
	      for(
		  int component_index=direction_index;
		  component_index<space_dimension;
		  component_index++
		  )
		{
		  work_observables_type.variance_of_polarity.push_back(0.0);
		  work_observables_type.variance_of_displacement.push_back(0.0);
		  // debug 3
		  // message ="Debug 3:"+io_method.longlongint_to_string(iteration_index);
		  // message+=io_method.int_to_string(direction_index);
		  // message+=io_method.int_to_string(component_index);
		  // message+=io_method.longint_to_string((long int)work_observables_type.variance_of_displacement.size());
		  //
		};
	    };
	};
      observables_types.push_back(work_observables_type);
    };
  //
  cell_volumes.clear();
  cell_displacements.clear();
  for(
      long int cell_index=0;
      cell_index<number_of_cells;
      cell_index++
      )
    {
      cell_volumes.push_back(0);
      for(
	  int direction_index=0;
	  direction_index<space_dimension;
	  direction_index++
	  )
	{
	  cell_polarities.push_back(0.0);
	  cell_displacements.push_back(0.0);
	  total_cell_displacements.push_back(0.0);
	};
    };
  //
  projection_vector.clear();
  for(
      int direction_index=0;
      direction_index<space_dimension;
      direction_index++
      )
    {
      projection_vector.push_back(0.0);
    };
  // Debug
  //  std::string message="Debug0:";
  // message+=io_method.longint_to_string((long int)number_of_cells)+",";
  //message+=io_method.longint_to_string((long int)cell_volumes.size());
  //io_method.standard_output(message);
  //
  work_polarity_vector.clear();
  work_polarity_variance.clear();
  work_polarity_pair.clear();
  work_polarity_pair_variance.clear();
  work_displacement_vector.clear();
  work_displacement_variance.clear();
  work_displacement_pair.clear();
  work_displacement_pair_variance.clear();
  work_norm.clear();
  work_square_norm.clear();
  work_position.clear();
  //
  work_polarity_pair.push_back(0.0);
  work_polarity_pair.push_back(0.0);
  work_polarity_pair_variance.push_back(0.0);
  work_polarity_pair_variance.push_back(0.0);
  work_displacement_pair.push_back(0.0);
  work_displacement_pair.push_back(0.0);
  work_displacement_pair_variance.push_back(0.0);
  work_displacement_pair_variance.push_back(0.0);
  //
  for(
      int direction_index=0;
      direction_index<space_dimension;
      direction_index++
      )
    {
      work_polarity_vector.push_back(0.0);
      work_displacement_vector.push_back(0.0);
      work_position.push_back(0.0);
      for(
	  int component_index=direction_index;
	  component_index<space_dimension;
	  component_index++
	  )
	{
	  work_polarity_variance.push_back(0.0);
	  work_displacement_variance.push_back(0.0);
	};
    };
  for(
      long int cell_index=0;
      cell_index<number_of_cells;
      cell_index++
      )
    {
      work_norm.push_back(0.0);
      work_square_norm.push_back(0.0);
      work_longitudinal.push_back(0.0);
      work_lateral.push_back(0.0);
    };
  //
  averaged_net_polarity.clear();
  variance_for_net_polarity.clear();
  averaged_variance_of_polarity.clear();
  averaged_net_displacement.clear();
  variance_for_displacement.clear();
  averaged_variance_of_displacement.clear();
  for(
      int direction_index=0;
      direction_index<space_dimension;
      direction_index++
      )
    {
      averaged_net_polarity.push_back(0.0);
      variance_for_net_polarity.push_back(0.0);
      averaged_net_displacement.push_back(0.0);
      variance_for_displacement.push_back(0.0);
      for(
	  int component_index=direction_index;
	  component_index<space_dimension;
	  component_index++
	  )
	{
	  averaged_variance_of_polarity.push_back(0.0);
	  averaged_variance_of_displacement.push_back(0.0);
	};
    };
  steps.clear();
  steps.push_back(0);
  for(
      long long int iteration_index=0;
      iteration_index<number_of_observations;
      iteration_index++
      )
    {
      steps.push_back(0);
    };
  averaged_number_of_cells=0.0;
};
//
void observables_type_system_class::load_data(
					      const state_system_class & state
					      )
{
  state.get_cell_volumes(cell_volumes);
  state.get_sorted_cell_polarities(
				   cell_polarities
				   );
  state.get_sorted_cell_displacements(
				      cell_displacements
				      );
  /*
  for(int direction_index=0;
      direction_index<space_dimension;
      direction_index++)
    {
    for(long int cell_index=0;
	cell_index<number_of_cells;
	cell_index++)
      {
	total_cell_displacements[direction_index+cell_index*space_dimension]
	  +=cell_displacements[number_of_cells*direction_index+cell_index];
      };
    };
  */
};
//
void observables_type_system_class::load_unsorted_data(
						       const state_system_class & state
						       )
{
  state.get_cell_volumes(cell_volumes);
  state.get_cell_polarities(
			    cell_polarities
			    );
  state.get_cell_displacements(
			       cell_displacements
			       );
  /*
  for(int direction_index=0;
      direction_index<space_dimension;
      direction_index++)
    {
    for(long int cell_index=0;
	cell_index<number_of_cells;
	cell_index++)
      {
	total_cell_displacements[direction_index+cell_index*space_dimension]
	  +=cell_displacements[number_of_cells*direction_index+cell_index];
      };
    };
  */
};
//
void observables_type_system_class::calculation(
						const state_system_class & state
						)
{
  if(
     iteration_counter<number_of_observations
     )
    {
      load_data(state);
      store_living_cell_number();
      for(int type_index=0;
	  type_index<number_of_cell_types;
	  type_index++)
	{
	  calculate_net_value(
			      type_index,
			      cell_polarities,
			      work_polarity_vector
			      );
	  calculate_variance(
			     type_index,
			     cell_polarities,
			     work_polarity_variance
			     );
	  calculate_projections(
				type_index,
				cell_polarities,
				work_polarity_vector,
				work_polarity_pair,
				work_polarity_pair_variance
				);
	  /*
	    io_method.standard_output(
	    "debug polarity :" + io_method.double_to_string(work_polarity_pair[0])
	    + '/' + io_method.double_to_string(work_polarity_pair[1])
	    + '/' + io_method.double_to_string(work_polarity_pair_variance[0])
	    + '/' + io_method.double_to_string(work_polarity_pair_variance[1])
	    + '/' + io_method.double_to_string(work_polarity_vector[0])
	    + '/' + io_method.double_to_string(work_polarity_vector[1])
	    + '/' + io_method.double_to_string(work_polarity_pair_variance[0])
	    );
	  */
	  calculate_net_value(
			      type_index,
			      cell_displacements,
			      work_displacement_vector
			      );
	  calculate_variance(
			     type_index,
			     cell_displacements,
			     work_displacement_variance
			     );
	  calculate_projections(
				type_index,
				cell_displacements,
				work_polarity_vector,
				work_displacement_pair,
				work_displacement_pair_variance
				);
	  /*
	    io_method.standard_output(
	    "debug displacement: "
	    +       io_method.double_to_string(work_displacement_pair[0])
	    + '/' + io_method.double_to_string(work_displacement_pair[1])
	    + '/' + io_method.double_to_string(work_displacement_pair_variance[0])
	    + '/' + io_method.double_to_string(work_displacement_pair_variance[1])
	    + '/' + io_method.double_to_string(work_polarity_vector[0])
	    + '/' + io_method.double_to_string(work_polarity_vector[1])
	    );
	  */
	  value_set(type_index);
	};
      track_push();
      iteration_counter++;
    };
};
//
void  observables_type_system_class::track(
					   const state_system_class & state
					   )
{
  state.get_sorted_cell_displacements(
				      cell_displacements
				      );
  for(
      int direction_index=0;
      direction_index<space_dimension;
      direction_index++
      )
    {
      for(
	  long int cell_index=0;
	  cell_index<number_of_cells;
	  cell_index++
	  )
	{
	  total_cell_displacements[
				   direction_index*number_of_cells
				   +cell_index
				   ]
	    +=cell_displacements[
				 direction_index*number_of_cells
				 +cell_index
				 ];
	};
    };
  track_steps++;
};
//
void  observables_type_system_class::track_push()
{
  //
  //  std::string message;
  //
  for(
      int type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      if(
	 observables_types[
			   type_index
			   ].end_pointer
	 >
	 observables_types[
			   type_index
			   ].start_pointer
	 &&
	 track_index<number_of_observations+1
	 )
	{
	  work_average_square_norm
	    =calculate_square_norm_average(
					   type_index,
					   total_cell_displacements,
					   observables_types[
							     type_index
							     ].number_of_living_cells[
										      iteration_counter
										      ]
					   );
	  calculate_net_value(
			      type_index,
			      total_cell_displacements,
			      work_position
			      );
	  observables_types[
			    type_index
			    ].total_displacement[
						 track_index
						 ]
	    =work_average_square_norm;
	  // debug
	  //	  message ="(";
	  //message+=io_method.int_to_string(type_index)+',';
	  //message+=io_method.double_to_string(work_average_norm)+',';
	  //message+=")";
	  //
      for(
	  int direction_index=0;
	  direction_index<space_dimension;
	  direction_index++
	  )
	{
	  observables_types[
			    type_index
			    ].averaged_position[
						(track_index)*space_dimension+direction_index
						]
	    =work_position[
			   direction_index
			   ];
	  observables_types[
			    type_index
			    ].represent_position[
						 (track_index)*space_dimension+direction_index
						 ]
	    =total_cell_displacements[
				      direction_index*number_of_cells
				      +observables_types[type_index].start_pointer
				      ];
	};
      //      io_method.standard_output(message);
	};
    };
  steps[track_index]=track_steps;
  track_index++;
};
//
void  observables_type_system_class::store_living_cell_number()
{
  for(int type_index=0;
      type_index<number_of_cell_types;
      type_index++)
    {
      for(
	  long int cell_index=observables_types[type_index].start_pointer;
	  cell_index<observables_types[type_index].end_pointer;
	  cell_index++
	  )
	{
	  if(cell_volumes[cell_index]>0)
	    {
	      observables_types[type_index].number_of_living_cells[iteration_counter]++;
	    };
	};
    };
};
//
void  observables_type_system_class::value_set(
					       const int & type_index
					       )
{
  for(int component_index=0;component_index<space_dimension;component_index++)
    {
      observables_types[
			type_index
			].net_polarity[
				       iteration_counter
				       *space_dimension
				       +component_index
				       ]
	=work_polarity_vector[
			      component_index
			      ];
      observables_types[
			type_index
			].net_displacement[
					   iteration_counter
					   *space_dimension
					   +component_index
					   ]
	=work_displacement_vector[
				  component_index
				  ];
    };
  for(
      int direction_index=0;
      direction_index<space_dimension;
      direction_index++
      )
    {
      for(
	  int component_index=direction_index;
	  component_index<space_dimension;
	  component_index++
	  )
	{
	  observables_types[
			    type_index
			    ].variance_of_polarity[
						   iteration_counter
						   *matrix_dimension
						   +direction_index
						   *(space_dimension-1)
						   +component_index
						   ]
	    =work_polarity_variance[
				    direction_index
				    *(space_dimension-1)
				    +component_index
				    ];
	  observables_types[
			    type_index
			    ].variance_of_displacement[
						       iteration_counter
						       *matrix_dimension
						       +direction_index
						       *(space_dimension-1)
						       +component_index
						       ]
		    =work_displacement_variance[
						direction_index
						*(space_dimension-1)
						+component_index
						];
	};
    };
  observables_types[
		    type_index
		    ].variance_longitudinal_polarity[
						     iteration_counter
						     ]
    =work_polarity_pair_variance[0];
  observables_types[
		    type_index
		    ].variance_lateral_polarity[
						iteration_counter
						]
    =work_polarity_pair_variance[1];
  observables_types[
		    type_index
		    ].longitudinal_displacement[
						iteration_counter
						]
    =work_displacement_pair[0];
  observables_types[
		    type_index
		    ].lateral_displacement[
					   iteration_counter
					   ]
    =work_displacement_pair[1];
  observables_types[
		    type_index
		    ].variance_longitudinal_displacement[
							 iteration_counter
							 ]
    =work_displacement_pair_variance[0];
  observables_types[
		    type_index
		    ].variance_lateral_displacement[
						    iteration_counter
						    ]
    =work_displacement_pair_variance[1];

};
//
void  observables_type_system_class::calculate_net_value(
							 const int & type_index,
							 const std::vector<double> & input_vectors,
							 std::vector<double> & output_vectors
							 )
{
  if(
     observables_types[type_index].end_pointer
     >observables_types[type_index].start_pointer
     )
    {
      sorted_vector_sum(
			output_vectors,
			input_vectors,
			observables_types[type_index].start_pointer,
			observables_types[type_index].end_pointer,
			number_of_cells,
			space_dimension
			);
      //
      for(int component_index=0;component_index<space_dimension;component_index++)
	{
	  output_vectors[component_index]
	    =output_vectors[component_index]
	    /((double)observables_types[type_index].number_of_living_cells[iteration_counter]);
	};
    };
};
//
void  observables_type_system_class::calculate_variance(
							 const int & type_index,
							 const std::vector<double> & input_vectors,
							 std::vector<double> & output_vectors
							 )
{
  /*
  for(
      int direction_index=0;
      direction_index<space_dimension;
      direction_index++
      )
    {
      for(
	  int component_index=direction_index;
	  component_index<space_dimension;
	  component_index++
	  )
	{
	  observables_types[type_index].variance_of_polarity.push_back(0.0);
	};
    };
  */
  if(observables_types[type_index].end_pointer>observables_types[type_index].start_pointer)
    {
      sorted_vector_inner_product(
				  output_vectors,
				  input_vectors,
				  input_vectors,
				  observables_types[type_index].start_pointer,
				  observables_types[type_index].end_pointer,
				  number_of_cells,
				  space_dimension
				  );
      //
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  for(int component_index=direction_index;component_index<space_dimension;component_index++)
	    {
	      output_vectors[
			     direction_index*(space_dimension-1)
			     +component_index
				 ]
		=output_vectors[
				direction_index*(space_dimension-1)
				+component_index
				]
		/((double)observables_types[type_index].number_of_living_cells[iteration_counter]);
	    };
	};
    };
};
//
void observables_type_system_class::calculate_projections(
							  const int & type_index,
							  const std::vector<double> & input_vectors,
							  const std::vector<double> & input_projection_direction,
							  std::vector<double> & output_pair,
							  std::vector<double> & output_pair_variance
							  )
{
  if(
     observables_types[type_index].end_pointer
     >observables_types[type_index].start_pointer
     )
    {
      work_double=norm(
		       input_projection_direction,
		       0,
		       space_dimension
		       );
      for(
	  int direction_index=0;
	  direction_index<space_dimension;
	  direction_index++
	  )
	{
	  projection_vector[direction_index]
	    =input_projection_direction[direction_index]/work_double;
	};
      //
      for(
	  long int cell_index=observables_types[type_index].start_pointer;
	  cell_index<observables_types[type_index].end_pointer;
	  cell_index++
	  )
	{
	  work_square_norm[cell_index]=0.0;
	  work_longitudinal[cell_index]=0.0;
	  work_lateral[cell_index]=0.0;
	};
      for(
	  int direction_index=0;
	  direction_index<space_dimension;
	  direction_index++
	  )
	{
	  for(
	      long int cell_index=observables_types[type_index].start_pointer;
	      cell_index<observables_types[type_index].end_pointer;
	      cell_index++
	      )
	    {
	      work_longitudinal[cell_index]
		+= input_vectors[direction_index*number_of_cells+cell_index]
		* projection_vector[direction_index];
	      work_square_norm[cell_index]
		+= input_vectors[direction_index*number_of_cells+cell_index]
		* input_vectors[direction_index*number_of_cells+cell_index];
	    };
	};
      for(
	  long int cell_index=observables_types[type_index].start_pointer;
	  cell_index<observables_types[type_index].end_pointer;
	  cell_index++
	  )
	{
	  work_lateral[cell_index]
	    +=std::sqrt(
			work_square_norm[cell_index]
			-work_longitudinal[cell_index]
			*work_longitudinal[cell_index]
			);
	};
      output_pair[0]=0.0;
      output_pair[1]=0.0;
      output_pair_variance[0]=0.0;
      output_pair_variance[1]=0.0;
      for(
	  long int cell_index=observables_types[type_index].start_pointer;
	  cell_index<observables_types[type_index].end_pointer;
	  cell_index++
	  )
	{
	  output_pair[0]+=work_longitudinal[cell_index];
	  output_pair[1]+=work_lateral[cell_index];
	  output_pair_variance[0]+=work_longitudinal[cell_index]
	    *work_longitudinal[cell_index];
	  output_pair_variance[1]+=work_lateral[cell_index]
	    *work_lateral[cell_index];
	};
      output_pair[0]
	=output_pair[0]
	/double(
		observables_types[type_index].end_pointer
		-observables_types[type_index].start_pointer
		);
      output_pair[1]
	=output_pair[1]
	/double(
		observables_types[type_index].end_pointer
		-observables_types[type_index].start_pointer
		);
      output_pair_variance[0]
	=output_pair_variance[0]
	/double(
		observables_types[type_index].end_pointer
		-observables_types[type_index].start_pointer
		);
      output_pair_variance[1]
	=output_pair_variance[1]
	/double(
		observables_types[type_index].end_pointer
		-observables_types[type_index].start_pointer
		);
    };
};
//
double observables_type_system_class::calculate_square_norm_average(
								    const int & type_index,
								    const std::vector<double> & input_vectors,
								    const double & diviser
								    )
{
  for(
      long int cell_index=observables_types[type_index].start_pointer;
      cell_index<observables_types[type_index].end_pointer;
      cell_index++
      )
    {
      work_square_norm[cell_index]=square_norm(
					       input_vectors,
					       cell_index,
					       space_dimension
					       );
    };
  return vector_sum(
		    work_square_norm,
		    observables_types[type_index].start_pointer,
		    observables_types[type_index].end_pointer
		    )
    /diviser;
};
//
inline double observables_type_system_class:: square_norm(
							  const std::vector<double> & input_vector,
							  const long int & index,
							  const long int & vector_size
							  ) const
{
  return std::inner_product(
			    input_vector.begin()+( index   *vector_size),
			    input_vector.begin()+((index+1)*vector_size),
			    input_vector.begin()+( index   *vector_size),
			    0.0
			    );
};
//
inline double observables_type_system_class:: norm(
						    const std::vector<double> & input_vector,
						    const long int & index,
						    const long int & vector_size
						    ) const
{
  return sqrt(
	      square_norm(input_vector,
			  index,
			  vector_size)
	      );
};
//
inline double observables_type_system_class:: vector_sum(
							 const std::vector<double> & input_vector,
							 const long int & start_pointer,
							 const long int & end_pointer
							 ) const
{
  return std::accumulate(
			 input_vector.begin()+start_pointer,
			 input_vector.begin()+end_pointer,
			 0.0
			 );
};
//
inline void observables_type_system_class:: sorted_vector_sum(
							       std::vector<double> & sum_vector,
							       const std::vector<double> & sorted_vectors,
							       const long int & start_pointer,
							       const long int & end_pointer,
							       const long int & number_of_vectors,
							       const int & vector_size
							       )
  const {
  for(int component_index=0;component_index<vector_size;component_index++)
    {
      sum_vector[component_index]
	=std::accumulate(
			 sorted_vectors.begin()
			 +(number_of_vectors*component_index)
			 +start_pointer,
			 sorted_vectors.begin()
			 +(number_of_vectors*component_index)
			 +end_pointer,
			 0.0
			 );
    };
};
//
inline void observables_type_system_class:: sorted_vector_inner_product(
									 std::vector<double> & sum_vector,
									 const std::vector<double> & sorted_vectors_1,
									 const std::vector<double> & sorted_vectors_2,
									 const long int & start_pointer,
									 const long int & end_pointer,
									 const long int & number_of_vectors,
									 const int & vector_size
									 )
  const {
  for(int direction_index=0;direction_index<vector_size;direction_index++)
    {
      for(int component_index=direction_index;component_index<vector_size;component_index++)
	{
	  sum_vector[direction_index*(vector_size-1)+component_index]
	    =std::inner_product(
				sorted_vectors_1.begin()+(number_of_vectors*direction_index)+start_pointer,
				sorted_vectors_1.begin()+(number_of_vectors*direction_index)+end_pointer,
				sorted_vectors_2.begin()+(number_of_vectors*component_index)+start_pointer,
				0.0
				);
	};
    };
};
//

//
void observables_type_system_class::output(
					   const std::vector<std::string> & parameter_titles,
					   const std::vector<double> & model_parameters
					   )
{
  //
  output_initialize();
  for(
      int type_index=0;
      type_index< number_of_cell_types;
      type_index++
      )
    {
      output_labels(
		    parameter_titles,
		    "obseravbles_for_type_",
		    type_index
		    );
  //
      if(observables_types[type_index].end_pointer
	 >observables_types[type_index].start_pointer)
	{
	  for(
	      long long int index=0;
	      index<number_of_observations;
	      index++
	      )
	    {
	      message="";
	      for(
		  int parameter_index=0;
		  parameter_index<(int)model_parameters.size();
		  parameter_index++
		  )
		{
		  message+=io_method.double_to_string(
						      model_parameters[parameter_index]
						      );
		  message+=" ";
		};
	      message 
		+= io_method.longint_to_string(
					       observables_types[type_index].number_of_living_cells[index]
					       );
	      message+= " ";
	      for(
		  int direction_index=0;
		  direction_index<space_dimension;
		  direction_index++
		  )
		{
		  message
		    +=io_method.double_to_string(
						  observables_types[type_index].net_polarity[
											     index*space_dimension
											     +direction_index
											     ]
						  );
		  message+=" ";
		};
	      //
	      for(
		  int direction_index=0;
		  direction_index<space_dimension;
		  direction_index++
		  )
		{
		  for(
		      int component_index=direction_index;
		      component_index<space_dimension;
		      component_index++
		      )
		    {
		      message
			+=io_method.double_to_string(
						      observables_types[type_index].variance_of_polarity[
													 index*matrix_dimension
													 +direction_index*(space_dimension-1)
													 +component_index
													 ]
						      );
		      message+=" ";
		    };
		};
	      // longitudinal & lateral components
	      message
		+=io_method.double_to_string(
					     observables_types[type_index].variance_longitudinal_polarity[
													  index
													  ]
					     );
	      message+=" ";
	      message
		+=io_method.double_to_string(
					     observables_types[type_index].variance_lateral_polarity[
												     index
												     ]
					     );
	      message+=" ";
	      //
	      for(
		  int direction_index=0;
		  direction_index<space_dimension;
		  direction_index++
		  )
		{
		  message
		    +=io_method.double_to_string(
						  observables_types[type_index].net_displacement[
												 index*space_dimension
												 +direction_index
												 ]
						      );
		  message+=" ";
		};
	      for(
		  int direction_index=0;
		  direction_index<space_dimension;
		  direction_index++
		  )
		{
		  for(
		      int component_index=direction_index;
		      component_index<space_dimension;
		      component_index++
		      )
		    {
		      message
			+=io_method.double_to_string(
						      observables_types[type_index].variance_of_displacement[
													     index*matrix_dimension
													     +direction_index*(space_dimension-1)
													     +component_index
													     ]
						      );
		      message+=" ";
		    };
		};
		  //		  message+=io_method.double_to_string(
		  //				      observables_types[type_index].total_displacement[index]
		  //				      );
	      // longitudinal & lateral
	      message
		+=io_method.double_to_string(
					     observables_types[type_index].longitudinal_displacement[
												     index
												     ]
					     );
	      message+=" ";
	      message
		+=io_method.double_to_string(
					     observables_types[type_index].variance_longitudinal_displacement[
													      index
													      ]
					     );
	      message+=" ";
	      message
		+=io_method.double_to_string(
					     observables_types[type_index].lateral_displacement[
												index
												]
					     );
	      message+=" ";
	      message
		+=io_method.double_to_string(
					     observables_types[type_index].variance_lateral_displacement[
													 index
													 ]
					     );
	      message+=" ";
	      //
	      io_method.output_message(message,filename);
	    };   
	};
	      io_method.output_message("",filename);
	      io_method.output_message("",filename);
    };
	  //
  for(
      int type_index=0;
      type_index< number_of_cell_types;
      type_index++
      )
    {
      output_labels(
		    parameter_titles,
		    "cell_track_for_type_",
		    type_index
		    );
      for(
	  long long int index=0;
	  index < number_of_observations+1;
	  index++
	  )
	{
	  message=io_method.longlongint_to_string(steps[index]);
	  message
	    +=" " 
	    +io_method.double_to_string(
					observables_types[type_index].total_displacement[index]
					);
	  for(
	      int direction_index=0;
	      direction_index<space_dimension;
	      direction_index++
	      )
	    {
	      message
		+=" " 
		+io_method.double_to_string(
					    observables_types[
							      type_index
							      ].averaged_position[
										  index*space_dimension
										  +direction_index
										  ]
					    );
	    };		       
	  message+=
	    " "+io_method.longint_to_string(
					    observables_types[
							      type_index
							      ].start_pointer
					    );
      for(
	  int direction_index=0;
	  direction_index<space_dimension;
	  direction_index++
	  )
	{
	  message
	    +=" " 
	    +io_method.double_to_string(
					observables_types[type_index].represent_position[
											 index*space_dimension
											 +direction_index
											 ]
					);
	};
      io_method.output_message(message,filename);
    };
      io_method.output_message("",filename);
      io_method.output_message("",filename);
    };
  //
  for(
      int type_index=0;
      type_index< number_of_cell_types;
      type_index++
      )
    {
      if(
	 observables_types[type_index].end_pointer
	 >observables_types[type_index].start_pointer
	 )
	{
	  output_initialize();
	  output_labels(
			parameter_titles,
			"averaged_obseravbles_for_type_",
			type_index
			);
	  calculation_average(type_index);
	  //
	  message="";
	  for(int parameter_index=0;parameter_index<(int)model_parameters.size();parameter_index++)
	    {
	      message+=io_method.double_to_string(model_parameters[parameter_index]);
	      message+=" ";
	    };
	  message 
	    += io_method.longint_to_string(
					   (long int)((double)averaged_number_of_cells/(double)number_of_observations)
					   );
	  //polarity component + variance # =space_dimension*2
	  message+= " ";
	  for(
	      int direction_index=0;
	      direction_index<space_dimension;
	      direction_index++
	      )
	    {
	      message 
		+= io_method.double_to_string(
					      averaged_net_polarity[direction_index]/(double)number_of_observations
					      );
	      message+= " ";
	      message 
		+= io_method.double_to_string(
					      variance_for_net_polarity[direction_index]/(double)number_of_observations
					      -std::pow(averaged_net_polarity[direction_index]/(double)number_of_observations,2.0)
					      );
	      message+= " ";
	    };
	  // absolute + its variance # = 2
	  message 
	    += io_method.double_to_string(
					  absolute_value_of_polarity/(double)number_of_observations
					  );
	  message+= " ";
	  message 
	    += io_method.double_to_string(
					  variance_for_absolute_value_of_polarity/(double)number_of_observations
					  -std::pow(absolute_value_of_polarity/(double)number_of_observations,2.0)
					  );
	  message+= " ";
	  // covariance space_dimension*(space_dimension+1)/2
	  for(
	      int direction_index=0;
	      direction_index<space_dimension;
	      direction_index++
	      )
	    {
	      for(
		  int component_index=direction_index;
		  component_index<space_dimension;
		  component_index++
		  )
		{
		  message 
		    += io_method.double_to_string(
						  averaged_variance_of_polarity[
										direction_index
										*(space_dimension-1)
										+component_index
										]
						  /(double)number_of_observations
						  );
		  message+= " ";
		};
	    };
	  message+= " ";
	  // longitudinal & lateral
	  message 
	    += io_method.double_to_string(
					  averaged_variance_of_longitudinal_polarity
					  /(double)number_of_observations
					  );
	  message+= " ";
	  message 
	    += io_method.double_to_string(
					  averaged_variance_of_lateral_polarity
					  /(double)number_of_observations
					  );
	  message+= " ";
	  //
	  for(
	      int direction_index=0;
	      direction_index<space_dimension;
	      direction_index++
	      )
	    {
	      message 
		+= io_method.double_to_string(
					      averaged_net_displacement[direction_index]
					      /(double)number_of_observations
					      );
	      message+= " ";
	      message 
		+= io_method.double_to_string(
					      variance_for_displacement[direction_index]
					      /(double)number_of_observations
					      -
					      std::pow(
						       averaged_net_displacement[direction_index]
						       /(double)number_of_observations
						       ,2.0)
					      );
	      message+= " ";
	    };
	  // absolute
	  message 
	    += io_method.double_to_string(
					  absolute_value_of_displacement/(double)number_of_observations
					  );
	  message+= " ";
	  message 
	    += io_method.double_to_string(
					  variance_for_absolute_value_of_displacement/(double)number_of_observations
					  -std::pow(absolute_value_of_displacement/(double)number_of_observations,2.0)
					  );
	  message+= " ";
	  //
	  for(
	      int direction_index=0;
	      direction_index<space_dimension;
	      direction_index++
	      )
	    {
	      for(
		  int component_index=direction_index;
		  component_index<space_dimension;
		  component_index++
		  )
		{
		  message 
		+= io_method.double_to_string(
					      averaged_variance_of_displacement[
										direction_index
										*(space_dimension-1)
										+component_index
										]
					      /(double)number_of_observations
					      );
		  message+= " ";
		};
	    }
	  message+= " ";
	  // longitudinal & lateral
	  message 
	    += io_method.double_to_string(
					  averaged_variance_of_longitudinal_displacement
					  /(double)number_of_observations
					  );
	  message+= " ";
	  message 
	    += io_method.double_to_string(
					  variance_for_longitudinal_displacement
					  /(double)number_of_observations
					  -
					  std::pow(
						   averaged_variance_of_longitudinal_displacement
						   /(double)number_of_observations,
						   2.0)
					  );
	  message+= " ";
	  message 
	    += io_method.double_to_string(
					  averaged_variance_of_lateral_displacement
					  /(double)number_of_observations
					  );
	  message+= " ";
	  message 
	    += io_method.double_to_string(
					  variance_for_lateral_displacement
					  /(double)number_of_observations
					  -
					  std::pow(
						   averaged_variance_of_lateral_displacement
						   /(double)number_of_observations,
						   2.0)
					  );
	  message+= " ";
	  message 
	    += io_method.double_to_string(
					  averaged_variance_of_longitudinal_displacement
					  /(double)number_of_observations
					  );
	  message+= " ";
	  message 
	    += io_method.double_to_string(
					  averaged_variance_of_lateral_displacement
					  /(double)number_of_observations
					  );
	  message+= " ";
	  io_method.output_message(message,filename);
	  io_method.output_message("",filename);
	  io_method.output_message("",filename);
	};
    };
};
//
void observables_type_system_class::output_initialize()
{
  for(
      int direction_index=0;
      direction_index<space_dimension;
      direction_index++
      )
    {
      averaged_net_polarity[direction_index]=0.0;
      averaged_net_displacement[direction_index]=0.0;
      for(
	  int component_index=direction_index;
	  component_index<space_dimension;
	  component_index++
	  )
	{
	  averaged_variance_of_polarity[direction_index*space_dimension+component_index]=0.0;
	  averaged_variance_of_displacement[direction_index*space_dimension+component_index]=0.0;
	};
    };
  averaged_number_of_cells=0.0;
};
//
void observables_type_system_class::output_labels(
						  const std::vector<std::string> & parameter_titles,
						  const std::string & filename_head,
						  const int & type_index
						  )
{
      filename=filename_head
	+io_method.int_to_string(type_index)
	+"_"
	+io_method.longint_to_format_string(sweep_step,"%04d")
	+".dat";
      message ="# ";
      for(int parameter_index=0;parameter_index<(int)parameter_titles.size();parameter_index++)
	{
	  message+=parameter_titles[parameter_index]+"/";
	};
      message+="Ncells/net_polarity("
	+io_method.int_to_string(space_dimension)
	+")/variance("
	+io_method.int_to_string(matrix_dimension)+")/net_displacement/variance";
      io_method.output_message(message,filename);
};
//
void  observables_type_system_class::calculation_average(const int & type_index)
 {
   std::string message;
   averaged_number_of_cells=0;
   absolute_value_of_polarity=0.0;
   absolute_value_of_displacement=0.0;
   variance_for_absolute_value_of_polarity=0.0;
   variance_for_absolute_value_of_displacement=0.0;
   //
   averaged_variance_of_longitudinal_polarity=0.0;
   averaged_variance_of_lateral_polarity=0.0;
   //
   averaged_longitudinal_displacement=0.0;
   variance_for_longitudinal_displacement=0.0;
   averaged_variance_of_longitudinal_displacement=0.0;
   averaged_lateral_displacement=0.0;
   variance_for_lateral_displacement=0.0;
   averaged_variance_of_lateral_displacement=0.0;
   //
   for(
       int direction_index=0;
       direction_index<space_dimension;
       direction_index++
       )
     {
       averaged_net_polarity[direction_index]=0.0;
       averaged_net_displacement[direction_index]=0.0;
       variance_for_net_polarity[direction_index]=0.0;
       variance_for_displacement[direction_index]=0.0;
	   for(
	       int component_index=direction_index;
	       component_index<space_dimension;
	       component_index++
	       )
	     {
	       averaged_variance_of_polarity[direction_index*space_dimension+component_index]=0.0;
	       averaged_variance_of_displacement[direction_index*space_dimension+component_index]=0.0;
	     };
     };
   for(
       long long int index=0;
       index<number_of_observations;
       index++
       )
     {
       averaged_number_of_cells
	 +=observables_types[type_index].number_of_living_cells[index];
       // for polarity
       work_double=0.0;
       for(
	   int direction_index=0;
	   direction_index<space_dimension;
	   direction_index++
	   )
	 {
	   work_double+=std::pow(observables_types[type_index].net_polarity[index*space_dimension+direction_index],2.0);
	 };
       absolute_value_of_polarity+=sqrt(work_double);
       variance_for_absolute_value_of_polarity+=std::pow(sqrt(work_double),2.0);
       
       //
       averaged_variance_of_longitudinal_polarity
	 +=sqrt(
		observables_types[type_index].variance_longitudinal_polarity[index]
		-
		work_double*work_double
		);
       averaged_variance_of_lateral_polarity
	 +=sqrt(observables_types[type_index].variance_lateral_polarity[index]);
       // for displacement
       work_double=0.0;
       for(
	   int direction_index=0;
	   direction_index<space_dimension;
	   direction_index++
	   )
	 {
	   work_double+=std::pow(observables_types[type_index].net_displacement[index*space_dimension+direction_index],2.0);
	 };
       absolute_value_of_displacement+=sqrt(work_double);
       variance_for_absolute_value_of_displacement+=std::pow(sqrt(work_double),2.0);
       //
       averaged_longitudinal_displacement
	 +=observables_types[type_index].longitudinal_displacement[index];
       variance_for_longitudinal_displacement
	 +=observables_types[type_index].longitudinal_displacement[index]
	 * observables_types[type_index].longitudinal_displacement[index];
       averaged_variance_of_longitudinal_displacement
	 +=sqrt(observables_types[type_index].variance_longitudinal_displacement[index]
		- observables_types[type_index].longitudinal_displacement[index]
		* observables_types[type_index].longitudinal_displacement[index]);
       averaged_lateral_displacement
	 +=observables_types[type_index].lateral_displacement[index];
       variance_for_lateral_displacement
	 +=observables_types[type_index].lateral_displacement[index]
	 * observables_types[type_index].lateral_displacement[index];
       averaged_variance_of_lateral_displacement
	 +=sqrt(observables_types[type_index].variance_lateral_displacement[index]
		- observables_types[type_index].lateral_displacement[index]
		* observables_types[type_index].lateral_displacement[index]);
       //
       for(
	   int direction_index=0;
	   direction_index<space_dimension;
	   direction_index++
	   )
	 {
	   averaged_net_polarity[direction_index]
	     +=observables_types[type_index].net_polarity[index*space_dimension+direction_index];
	   variance_for_net_polarity[direction_index]
	     +=std::pow(observables_types[type_index].net_polarity[index*space_dimension+direction_index],2.0);
	   averaged_net_displacement[direction_index]
	     +=observables_types[type_index].net_displacement[index*space_dimension+direction_index];
	   variance_for_displacement[direction_index]
	     +=std::pow(observables_types[type_index].net_displacement[index*space_dimension+direction_index],2.0);
	   for(
	       int component_index=direction_index;
	       component_index<space_dimension;
	       component_index++
	       )
	     {
	       averaged_variance_of_polarity[direction_index*(space_dimension-1)+component_index]
		 +=observables_types[type_index].variance_of_polarity[
								      index*matrix_dimension
								      +direction_index*(space_dimension-1)
								      +component_index
								      ];
	       averaged_variance_of_displacement[
						 direction_index*(space_dimension-1)
						 +component_index
						 ]
		 +=observables_types[type_index].variance_of_displacement[
									  index*matrix_dimension
									  +direction_index*(space_dimension-1)
									  +component_index
									  ];
	     };
	 };
     };
 };
