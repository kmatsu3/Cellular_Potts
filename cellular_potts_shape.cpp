#include "cellular_potts_shape.hpp"
shape_class::shape_class()
{
  polar_moments.clear();
};
//
shape_tensor_class::shape_tensor_class()
{
  shape_moment_tensors.clear();
};
//
shape_system_class::shape_system_class()
{
};
//
shape_observable_for_type_class::shape_observable_for_type_class()
{
  shape_observables.clear();
};
//
void shape_class::initialize(
			     long int & number_of_cells
			     )
{
  polar_moments.clear();
  for(
      long int moment_index=0;
      moment_index<number_of_cells;
      moment_index++
      )
    {
      polar_moments.push_back(-1.0);
    };
};
//
void shape_tensor_class::initialize(
				    long int & number_of_cells,
				    int & space_dimension
				    )
{
  shape_moment_tensors.clear();
  for(
      long int moment_index=0;
      moment_index<number_of_cells*space_dimension*space_dimension;
      moment_index++
      )
    {
      shape_moment_tensors.push_back(-1.0);
    };
};
//
void shape_system_class::initialize(
				    const model_parameters_cellular_potts_class & model,
				    const cell_system_class & cell_system,
				    const type_system_class & cell_type_system,
				    site_system_class & site_system,
				    int & input_sweep_step
				    )
{
  number_of_polarities=1;
  number_of_components=3;
  counter=0;
  polar_gen_flag=true;
  sweep_step=input_sweep_step;
  number_of_cells=model.get_number_of_cells();
  number_of_cell_types=model.get_number_of_cell_types();
  number_of_sites=model.get_number_of_sites();
  space_dimension=model.get_space_dimension();
  number_of_cell_types=model.get_number_of_cell_types();
  buffer_cell=cell_system.get_buffer_cell();
  for(
      int direction_index=0;
      direction_index<space_dimension;
      direction_index++
      )
    {
      work_long_int.push_back(0);
      relative_coordinate.push_back(0.0);
    };
  shape_class temporary_shapes;
  temporary_shapes.initialize(
			      number_of_cells
			      );
  for(
      int polarity_index=0;
      polarity_index<number_of_polarities*number_of_components;
      polarity_index++
      )
    {
      shapes.push_back(temporary_shapes);
    };
  shape_tensors.initialize(
			   number_of_cells,
			   space_dimension
			   );
  //
  shape_observable_for_type_class temporary_observables;
  temporary_observables.shape_observables.clear();
  for(
      int data_index=0;
      data_index<number_of_polarities*number_of_components;
      data_index++
	  )
    {
      temporary_observables.shape_observables.push_back(0.0);
    };
  long int number_of_cells_for_type_pointer;
  //
  polar_observables.clear();
  polar_observable_variances.clear();
  squared_polar_observables.clear();
  shape_traces.clear();
  shape_trace_variances.clear();
  squared_shape_traces.clear();
  shape_traceless_determinants.clear();
  shape_traceless_determinant_variances.clear();
  squared_shape_traceless_determinants.clear();
  //
  for(
      int type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      polar_observables.push_back(temporary_observables);
      //
      polar_observable_variances.push_back(temporary_observables);
      //
      squared_polar_observables.push_back(temporary_observables);
      //
    };
  number_of_cells_for_type_pointer=0;
  // debug
  std::string message;
  //
  for(
      int type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      polar_observables[type_index].start_pointer
	=number_of_cells_for_type_pointer;
      //
      polar_observable_variances[type_index].start_pointer
	=number_of_cells_for_type_pointer;
      //
      squared_polar_observables[type_index].start_pointer
	=number_of_cells_for_type_pointer;
      //
      // debug
      //      message="debug:mem2:" + io_method.longint_to_string(number_of_cells_for_type_pointer);
      // debug
      if(cell_type_system.get_number_of_cells(type_index)>0)
	{
	  number_of_cells_for_type_pointer
	    =number_of_cells_for_type_pointer
	    +cell_type_system.get_number_of_cells(type_index);
	};
      // debug
      //      message += ","+ io_method.longint_to_string(number_of_cells_for_type_pointer);
      //io_method.standard_output(message);
      // debug
      polar_observables[type_index].end_pointer
	=number_of_cells_for_type_pointer;
      polar_observable_variances[type_index].end_pointer
	=number_of_cells_for_type_pointer;
      squared_polar_observables[type_index].end_pointer
	=number_of_cells_for_type_pointer;
      //
    };
  //
  temporary_observables.shape_observables.clear();
  temporary_observables.shape_observables.push_back(0.0);
  for(
      int type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      shape_traces.push_back(temporary_observables);
      //
      shape_trace_variances.push_back(temporary_observables);
      //
      squared_shape_traces.push_back(temporary_observables);
      //
      shape_traceless_determinants.push_back(temporary_observables);
      //
      shape_traceless_determinant_variances.push_back(temporary_observables);
      //
      squared_shape_traceless_determinants.push_back(temporary_observables);
    };
  number_of_cells_for_type_pointer=0;
  for(
      int type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      shape_traces[type_index].start_pointer
	=number_of_cells_for_type_pointer;
      //
      shape_trace_variances[type_index].start_pointer
	=number_of_cells_for_type_pointer;
      //
      squared_shape_traces[type_index].start_pointer
	=number_of_cells_for_type_pointer;
      //
      shape_traceless_determinants[type_index].start_pointer
	=number_of_cells_for_type_pointer;
      //
      shape_traceless_determinant_variances[type_index].start_pointer
	=number_of_cells_for_type_pointer;
      //
      squared_shape_traceless_determinants[type_index].start_pointer
	=number_of_cells_for_type_pointer;
      //
      //
      if(cell_type_system.get_number_of_cells(type_index)>0)
	{
	  number_of_cells_for_type_pointer
	    =number_of_cells_for_type_pointer
	    +cell_type_system.get_number_of_cells(type_index);
	};
      //
      //
      shape_traces[type_index].end_pointer
	=number_of_cells_for_type_pointer;
      //
      shape_trace_variances[type_index].end_pointer
	=number_of_cells_for_type_pointer;
      //
      squared_shape_traces[type_index].end_pointer
	=number_of_cells_for_type_pointer;
      //
      shape_traceless_determinants[type_index].end_pointer
	=number_of_cells_for_type_pointer;
      //
      shape_traceless_determinant_variances[type_index].end_pointer
	=number_of_cells_for_type_pointer;
      //
      squared_shape_traceless_determinants[type_index].end_pointer
	=number_of_cells_for_type_pointer;
    };
  tmp_matrix.clear();
  tmp_vector.clear();
  for(
      int direction_index=0;
      direction_index<space_dimension;
      direction_index++
      )
    {
      tmp_vector.push_back(0.0);
      for(
	  int sub_index=0;
	  sub_index<space_dimension;
	  sub_index++
	  )
	{
	  tmp_matrix.push_back(0.0);
	};
    };
  //
  traces.clear();
  traceless_determinants.clear();
  for(
      long int cell_index=0;
      cell_index < number_of_cells;
      cell_index++
      )
    {
      traces.push_back(0.0);
      traceless_determinants.push_back(0.0);
    };
};
//
void shape_system_class::calculate_polar_correlation(
						     const state_system_class & state,
						     const cell_system_class & cell_system
						     )
{
  if(polar_gen_flag)
    {
      for(
	  int polarity_index=0;
	  polarity_index<number_of_polarities;
	  polarity_index++
	  )
	{
	  for(
	      long int cell_index=0;
	      cell_index<number_of_cells;
	      cell_index++
	      )
	    {
	      shapes[polarity_index].polar_moments[cell_index]
		=0.0;
	    };
	  for(
	      long long int site_index=0;
	      site_index<number_of_sites;
	      site_index++
	      )
	    {
	      if(state.configuration[site_index]!=cell_system.get_buffer_cell())
		{
		  shapes[number_of_components*polarity_index ].polar_moments[
									     state.configuration[site_index]
									     ]
		    =+state.polar_product[site_index];
		  shapes[number_of_components*polarity_index+1].polar_moments[
									      state.configuration[site_index]
									      ]
		    =+pow(state.polar_product[site_index],2.0);
		  shapes[number_of_components*polarity_index+2].polar_moments[
									      state.configuration[site_index]
									      ]
		    =+(1.0-pow(state.polar_product[site_index],2.0));
		};
	    };
	  //	  io_method.standard_output(io_method.longint_to_string(number_of_cells));
	  for(
	      long int cell_index=0;
	      cell_index<number_of_cells;
	      cell_index++
	      )
	    {
	      shapes[number_of_components*polarity_index].polar_moments[
									cell_index
									]
		=shapes[number_of_components*polarity_index].polar_moments[
									   cell_index
									   ]
		/state.cell_volumes[cell_index];
	      shapes[number_of_components*polarity_index+1].polar_moments[
									  cell_index
									  ]
		=shapes[number_of_components*polarity_index+1].polar_moments[
									     cell_index
									     ]
		/state.cell_volumes[cell_index];
	      shapes[number_of_components*polarity_index+2].polar_moments[
									  cell_index
									  ]
		=shapes[number_of_components*polarity_index+2].polar_moments[
									     cell_index
									     ]
		/state.cell_volumes[cell_index];
	    };      
	};
      //   polar_gen_flag=false;
    };
  accumlate_polar_correlation();
};
//
void shape_system_class::accumlate_polar_correlation()
{
  double observable;
  double squared_observable;
  counter=counter+1;
  for(
      int type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      for(
	  int data_index=0;
	  data_index < number_of_polarities*number_of_components;
	  data_index++
	  )
	{
	  observable=0.0;
	  squared_observable=0.0;
	  for(
	      long int cell_index=polar_observables[type_index].start_pointer;
	      cell_index < polar_observables[type_index].end_pointer;
	      cell_index++
	      )
	    {
	      //   io_method.standard_output("debug mem:"+ io_method.int_to_string(data_index)+","
	      //			+io_method.int_to_string(cell_index));
	      observable
		=observable
		+shapes[data_index].polar_moments[
						  cell_index
						  ];
	      squared_observable
		=squared_observable
		+pow(
		     shapes[data_index].polar_moments[
						      cell_index
						      ],
		     2.0
		     );
	    };
	  observable=observable/((double)number_of_cells);
	  squared_observable=squared_observable/((double)number_of_cells);
	  polar_observables[type_index].shape_observables[data_index]
	    =polar_observables[type_index].shape_observables[data_index]
	    +observable;
	  polar_observable_variances[type_index].shape_observables[data_index]
	    =polar_observable_variances[type_index].shape_observables[data_index]
	    +squared_observable-pow(observable,2.0);
	  squared_polar_observables[type_index].shape_observables[data_index]
	    =squared_polar_observables[type_index].shape_observables[data_index]
	    +pow(observable,2.0);
	};
    };
};
//
void shape_system_class::output_shape_characters(
						 const std::vector<std::string> & parameter_titles,
						 const std::vector<double> & model_parameters
						 )
{
  std::string message;
  for(
      int type_index=0;
      type_index< number_of_cell_types;
      type_index++
      )
    {
      output_labels(
		    parameter_titles,
		    "shape_character_for_type_",
		    type_index
		    );
      message = "";
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
      message+=io_method.double_to_string(
					  shape_traces[type_index].shape_observables[0]
					  /(double)counter
					  );
      message+=" ";
      message+=io_method.double_to_string(
					  shape_trace_variances[type_index].shape_observables[0]
					  /(double)counter
					  );
      message+=" ";
      if(counter>0)
	{
	  message+=io_method.double_to_string(
					      (
					       squared_shape_traces[type_index].shape_observables[0]
						   -shape_traces[type_index].shape_observables[0]
						   *shape_traces[type_index].shape_observables[0]
						   /(double)counter
						   )/(double)counter
						  );
	      message+=" ";
	};
      message+=io_method.double_to_string(
					  shape_traceless_determinants[type_index].shape_observables[0]
					  /(double)counter
					  );
      message+=" ";
      message+=io_method.double_to_string(
					  shape_traceless_determinant_variances[type_index].shape_observables[0]
					  /(double)counter
					  );
      message+=" ";
      if(counter>0)
	{
	  message+=io_method.double_to_string(
					      (
					       squared_shape_traceless_determinants[type_index].shape_observables[0]
						   -shape_traceless_determinants[type_index].shape_observables[0]
						   *shape_traceless_determinants[type_index].shape_observables[0]
						   /(double)counter
						   )/(double)counter
						  );
	      message+=" ";
	};
      io_method.output_message(message, filename);
      io_method.output_message("",filename);
      io_method.output_message("",filename);
    };
};
//
void shape_system_class::output_polar_correlation(
						  const std::vector<std::string> & parameter_titles,
						  const std::vector<double> & model_parameters
						  )
{
  std::string message;
  for(
      int type_index=0;
      type_index< number_of_cell_types;
      type_index++
      )
    {
      output_labels(
		    parameter_titles,
		    "shape_polar_correlation_for_type_",
		    type_index
		    );
      message = "";
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
      for(
	  int data_index=0;
	  data_index<number_of_polarities*number_of_components;
	  data_index++
	  )
	{
	  message+=io_method.double_to_string(
					      polar_observables[type_index].shape_observables[data_index]
					      /(double)counter
					      );
	  message+=" ";
	  message+=io_method.double_to_string(
					      polar_observable_variances[type_index].shape_observables[data_index]
					      /(double)counter
					      );
	  message+=" ";
	  if(counter>0)
	    {
	      message+=io_method.double_to_string(
						  (
						   squared_polar_observables[type_index].shape_observables[data_index]
						   -polar_observables[type_index].shape_observables[data_index]
						   *polar_observables[type_index].shape_observables[data_index]
						   /(double)counter
						   )/(double)counter
						  );
	      message+=" ";
	    };
	};
      io_method.output_message(message, filename);
      io_method.output_message("",filename);
      io_method.output_message("",filename);
    };
};
//
void shape_system_class::output_labels(
				       const std::vector<std::string> & parameter_titles,
				       const std::string & filename_head,
				       const int & type_index
				       )
{
  std::string message;
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
  message+="raw/variance/fluctuation/parallel/variance/fluctuation/lateral/variance/fluctuation";
  io_method.output_message(message,filename);
};
//
void shape_system_class::calculate_shape_tensor(
						const state_system_class & state,
						const model_parameters_cellular_potts_class & model,
						site_system_class & site_system,
						const cell_system_class & cell_system
						
						)
{
  for(
      long int cell_index=0;
      cell_index<number_of_cells;
      cell_index++
      )
    {
      for(
	  int direction_index=0;
	  direction_index<space_dimension;
	  direction_index++
	  )
	{
	  for(
	      int sub_index=direction_index;
	      sub_index<space_dimension;
	      sub_index++
	      )
	    {
	      shape_tensors.shape_moment_tensors[
						 cell_index*space_dimension*space_dimension
						 +direction_index*space_dimension+sub_index
						 ]=0.0;
		};
	};
    };
  long int cell_index;
  for(
      long long int site_index=0;
      site_index<number_of_sites;
      site_index++
      )
    {
      cell_index = state.configuration[site_index];
      if(cell_index!=buffer_cell)
	{
	  calculate_relative_coordinate(
					state,
					model,
					site_system,
					site_index,	
					cell_index
				    );
	  for(
	      int direction_index=0;
	      direction_index<space_dimension;
	      direction_index++
	      )
	    {
	      for(
	      int sub_index=direction_index;
	      sub_index<space_dimension;
	      sub_index++
		  )
		{
		  shape_tensors.shape_moment_tensors[
						     cell_index*space_dimension*space_dimension
						     +direction_index*space_dimension+sub_index
						     ]
		    = shape_tensors.shape_moment_tensors[
							 cell_index*space_dimension*space_dimension
							 +direction_index*space_dimension+sub_index
							 ]
		    +relative_coordinate[direction_index]
		    *relative_coordinate[sub_index];
		};
	    };
	};
    };
  for(
      cell_index = 0;
      cell_index < number_of_cells;
      cell_index++
      )
    {
      for(
	  int direction_index=0;
	  direction_index<space_dimension;
	  direction_index++
	  )
	{
	  for(
	      int sub_index=direction_index;
	      sub_index<space_dimension;
	      sub_index++
	      )
	    {
	      shape_tensors.shape_moment_tensors[
						 cell_index
						 *space_dimension*space_dimension
						 +direction_index*space_dimension+sub_index
						 ]
		= shape_tensors.shape_moment_tensors[
						     cell_index
						     *space_dimension*space_dimension
						     +direction_index*space_dimension+sub_index
						     ]
		/ (double)(state.cell_volumes[cell_index]*state.cell_volumes[cell_index])*4.0*M_PI;
	    };
	};
    };
  accumlate_shape_tensor();
};
//
void shape_system_class::accumlate_shape_tensor() 
{
  calculate_characters(
		       shape_tensors
		       );
  accumlate_tensor(
		   traces,
		   shape_traces,
		   shape_trace_variances,
		   squared_shape_traces
		   );
  //debug
  //std::vector<double>::iterator debug_index=traceless_determinants.begin();
  //size_t index_num;
  //while(
  //	debug_index!=traceless_determinants.end()
  //	)
  // {
  //   index_num= std::distance(traceless_determinants.begin(),debug_index);
  //  io_method.standard_output(
  //				"debug :" 
  //				+ io_method.longint_to_string(index_num) + ','
  //				+ io_method.double_to_string(*debug_index)
  //				);
  //  debug_index++;
  //};
  // debug
  accumlate_tensor(
		   traceless_determinants,
		   shape_traceless_determinants,
		   shape_traceless_determinant_variances,
		   squared_shape_traceless_determinants
		   );
};
//
void shape_system_class::calculate_characters(
					      shape_tensor_class & tensors
					      )
{
  bool zero_determinant;
  double ratio;
  for(
      long int cell_index=0;
      cell_index < number_of_cells;
      cell_index++
      )
    {
      traces[cell_index]=0.0;
      for(
	  int direction_index=0;
	  direction_index<space_dimension;
	  direction_index++
	  )
	{
	  traces[cell_index]
	    +=tensors.shape_moment_tensors[
					   cell_index
					   *space_dimension*space_dimension
					   +direction_index
					   *space_dimension
					   +direction_index
					   ];
	};
      traceless_determinants[cell_index]=0.0;
      for(
	  int direction_index=0;
	  direction_index<space_dimension;
	  direction_index++
	  )
	{
	  for(
	      int sub_index=0;
	      sub_index<space_dimension;
	      sub_index++
	      )
	    {
	      if(direction_index<=sub_index)
		{
		  tmp_matrix[direction_index*space_dimension+sub_index]
		    =tensors.shape_moment_tensors[
						  cell_index
						  *space_dimension*space_dimension
						  +direction_index
						  *space_dimension
						  +sub_index
						  ];
		}
	      else
		{
		  tmp_matrix[direction_index*space_dimension+sub_index]
		    =tensors.shape_moment_tensors[
						  cell_index
						  *space_dimension*space_dimension
						  +sub_index
						  *space_dimension
						  +direction_index
						  ];
		};
	      if(direction_index==sub_index)
		{
		  tmp_matrix[direction_index*space_dimension+sub_index]
		    -=traces[cell_index]/(double)space_dimension;
		};
	    };
	};
      // debug
      //      io_method.standard_output("debug:("+io_method.longint_to_string(cell_index)
      //				+"):"+io_method.double_to_string(traces[cell_index]));
      //   io_method.standard_output("n_debug: matrix ("+io_method.double_to_string(tmp_matrix[0])
      //				+","+io_method.double_to_string(tmp_matrix[1])+")");
      // io_method.standard_output("n_debug: matrix ("+io_method.double_to_string(tmp_matrix[2])
      //			+","+io_method.double_to_string(tmp_matrix[3])+")");	  
      // debug
      for(
	  int direction_index=0;
	  direction_index<space_dimension;
	  direction_index++
	  )
	{
	  zero_determinant=false;
	  if(
	     tool.finite_abs_check(
				   tmp_matrix[direction_index*space_dimension+direction_index]
				   )
	     !=
	     tool.get_value("true")
	     )
	    {
	      for(
		  int sub_index=direction_index+1;
		  sub_index<space_dimension;
		  sub_index++
		  )
		{
		  if(
		     tool.finite_abs_check(
					   tmp_matrix[direction_index*space_dimension+sub_index]
					   )
		     ==
		     tool.get_value("true")
		     )
		    {
		      for(
			  int component_index=0;
			  component_index<space_dimension;
			  component_index++
			  )
			{
			  tmp_vector[component_index]
			    =tmp_matrix[direction_index*space_dimension+component_index];
			};
		      for(
			  int component_index=0;
			  component_index<space_dimension;
			  component_index++
			  )
			{
			  tmp_matrix[direction_index*space_dimension+component_index]
			    =tmp_matrix[sub_index*space_dimension+component_index];
			};
		      for(
			  int component_index=0;
			  component_index<space_dimension;
			  component_index++
			  )
			{
			  tmp_matrix[sub_index*space_dimension+component_index]
			    =tmp_vector[component_index];
			};		      
		      break;
		    };
		};
	    };
	  //
	  if(
	     tool.finite_abs_check(
				   tmp_matrix[direction_index*space_dimension+direction_index]
				   )
	     !=
	     tool.get_value("true")
	     )
	    {
	      traceless_determinants[cell_index]=0.0;
	      zero_determinant=true;
	      break;
	    };
	  for(
	      int sub_index=direction_index+1;
	      sub_index<space_dimension;
	      sub_index++
	      )
	    {
	      ratio=tmp_matrix[sub_index*space_dimension+direction_index]
		/tmp_matrix[direction_index*space_dimension+direction_index];
	      for(
		  int component_index=direction_index;
		  component_index<space_dimension;
		  component_index++
		  )
		{
		  tmp_matrix[sub_index*space_dimension+component_index]
		    -=tmp_matrix[direction_index*space_dimension+component_index]
		    *ratio;
		};
	    };
	};
      if(!zero_determinant)
	{
	  //
	  traceless_determinants[cell_index]=1.0;
	  for(
	      int direction_index=1.0;
	      direction_index<space_dimension;
	      direction_index++
	      )
	    {
	      traceless_determinants[cell_index]
		*=tmp_matrix[
			     direction_index*space_dimension
			     +direction_index
			     ];
	    };
	};
      //debug
      //      io_method.standard_output("ok" + io_method.longint_to_string(cell_index));
      //io_method.standard_output("debug: T_matrix ("+io_method.double_to_string(tmp_matrix[0])
      //				+","+io_method.double_to_string(tmp_matrix[1])+")");
      //io_method.standard_output("debug: T_matrix ("+io_method.double_to_string(tmp_matrix[2])
      // 				+","+io_method.double_to_string(tmp_matrix[3])+")");
      //io_method.standard_output("debug: det=" + io_method.double_to_string(traceless_determinants[cell_index]));
      // debug
    };
};
//
void shape_system_class::accumlate_tensor(
					  std::vector<double> & characters,
					  std::vector<shape_observable_for_type_class> & observables,
					  std::vector<shape_observable_for_type_class> & observable_variances,
					  std::vector<shape_observable_for_type_class> & squared_observables
					  )
{
  double tmp_observable=0.0;
  double tmp_observable_variance=0.0;
  long int counter=0;
  for(
      int type_index =0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      counter=0;
      for(
	  long int cell_index = observables[type_index].start_pointer;
	  cell_index < observables[type_index].end_pointer;
	  cell_index++
	  )
	{
	  tmp_observable
	    +=characters[
			 cell_index
			 ];
	  tmp_observable_variance
	    +=pow(
		  characters[
			     cell_index
			     ]
		  , 2.0
		  );
	  counter++;
	};
      if(counter!=0)
	{
	  observables[type_index].shape_observables[0]
	    +=tmp_observable/(double)counter;
	  observable_variances[type_index].shape_observables[0]
	    +=tmp_observable_variance/(double)counter
	    -pow(tmp_observable,2.0)/(double)(counter*counter);
	  squared_observables[type_index].shape_observables[0]
	    +=pow(tmp_observable/(double)counter,2.0);
	};
    };
};
