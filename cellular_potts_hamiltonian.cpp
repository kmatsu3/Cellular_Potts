#include "cellular_potts_hamiltonian.hpp"
#include "cellular_potts_state.hpp"
//
double local_state_class::get_local_isotropic_adhesion_energy() 
{
  double return_value=0.0;
  int neighbor_index;
  for(neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
    {
      if(cell!=neighbor_cells[neighbor_index])
	{
	  work_vector_double_for_neighbors[neighbor_index]
	    =isotropic_adhesion_couplings[cell_type][neighbor_types[neighbor_index]];
	}
      else
	{
	  work_vector_double_for_neighbors[neighbor_index]=0.0;
	};
    };
  //
  return_value
    =accumulate(
		work_vector_double_for_neighbors.begin(),
		work_vector_double_for_neighbors.end(),
		0.0
		);
  //
  return return_value;
};
//
void local_state_class::set_product_polarity()
{
  product_polarity=0.0;
  //  int debug_index;
  if(
     (cell_type!=buffer_type)
     &&
     (iszero_norm(relative_coordinates)==false_value)
     )
    {
      product_polarity=normalized_product(polarities,relative_coordinates);
//      std::string message="debug:";
//      message+=io_method.double_to_string(relative_coordinates[0])+',';
//      message+=io_method.double_to_string(relative_coordinates[0]);
//      io_method.standard_output(message);
    };
};
//
double local_state_class::get_local_polarity_driving_energy()
{
  return product_polarity*polarity_sensitivity[cell_type]
 +(product_polarity*product_polarity-0.5)*quadrapolarity_sensitivity[cell_type];
};
//
void local_state_class::set_product_external_field()
{
  product_external_field=0.0;
  if(
     (cell_type!=buffer_type&&field_flag)
     &&
     (iszero_norm(relative_coordinates)==false_value)
     )
    {
      product_external_field=normalized_product(external_field,relative_coordinates);
    };
};
//
double local_state_class::get_local_field_driving_energy()
{
  return product_external_field*field_sensitivity[cell_type];
};
//
double local_state_class::get_local_dipolar_adhesion_energy()
{
  double return_value=0.0;
  double polar_product_neighbor=0.0;
  double work_product_polar=0.0;
  double work_product_polar_neighbor=0.0;
  double work_product_quadrapolar=0.0;
  double work_product_quadrapolar_neighbor=0.0;
  int neighbor_index;
  //  int debug_index;
  if(
     (cell_type!=buffer_type)
     &&
     (iszero_norm(relative_coordinates)==false_value)
     )
    {
      work_product_polar=product_polarity+dipolar_adhesion_basals[cell_type];
      work_product_quadrapolar=product_polarity*product_polarity+quadrapolar_adhesion_basals[cell_type];
      //
      for(neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
	{
	  if(
	     (
	      neighbor_types[neighbor_index]
	      ==
	      buffer_type
	      )
	     ||
	     (
	      cell
	      ==
	      neighbor_cells[neighbor_index]
	      )
	     ||
	     (iszero_norm(neighbor_relative_coordinates[neighbor_index])
	      ==
	      true_value
	      )
	     )
	    {
	      work_vector_double_for_neighbors[neighbor_index]=0.0;
	    }
	  else
	    {
	      //
	      polar_product_neighbor=normalized_product(
							neighbor_polarities[neighbor_index],
							neighbor_relative_coordinates[neighbor_index]
							);
	      work_product_polar_neighbor=polar_product_neighbor
		+dipolar_adhesion_basals[neighbor_types[neighbor_index]];
	      work_product_quadrapolar_neighbor=polar_product_neighbor*polar_product_neighbor
		+quadrapolar_adhesion_basals[neighbor_types[neighbor_index]];
	      work_vector_double_for_neighbors[neighbor_index]
		=dipolar_adhesion_couplings[cell_type][neighbor_types[neighbor_index]]
		*work_product_polar*work_product_polar_neighbor
		+quadrapolar_adhesion_couplings[cell_type][neighbor_types[neighbor_index]]
		*work_product_quadrapolar*work_product_quadrapolar_neighbor;
	      //
	    };
	};
      return_value
	=accumulate(
		    work_vector_double_for_neighbors.begin(),
		    work_vector_double_for_neighbors.end(),
		    0.0
		    );
    };
  //
  return return_value;
};
//
/*
double local_state_class::get_local_dipolar_adhesion_energy()
{
  double return_value=0.0;
  double work_product=0.0;
  int neighbor_index;
  //  int debug_index;
  if(
     (cell_type!=buffer_type)
     &&
     (iszero_norm(relative_coordinates)==false_value)
     )
    {
      work_product=product_polarity+dipolar_adhesion_basals[cell_type];
      //
      for(neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
	{
	  if(
	     (
	      neighbor_types[neighbor_index]
	      ==
	      buffer_type
	      )
	     ||
	     (
	      cell
	      ==
	      neighbor_cells[neighbor_index]
	      )
	     ||
	     (iszero_norm(neighbor_relative_coordinates[neighbor_index])
	      ==
	      true_value
	      )
	     )
	    {
	      work_vector_double_for_neighbors[neighbor_index]=0.0;
	    }
	  else
	    {
	      work_vector_double_for_neighbors[neighbor_index]
		=dipolar_adhesion_couplings[cell_type][neighbor_types[neighbor_index]]
		*work_product
	        *(
		  normalized_product(
				     neighbor_polarities[neighbor_index],
				     neighbor_relative_coordinates[neighbor_index]
				     )
		  +dipolar_adhesion_basals[neighbor_types[neighbor_index]]
		  );
	      //
	    };
	};
      return_value
	=accumulate(
		    work_vector_double_for_neighbors.begin(),
		    work_vector_double_for_neighbors.end(),
		    0.0
		    );
    };
  //
  return return_value;
};
*/
//
inline int local_state_class::iszero_norm(
					  const std::vector<double> & vector_data
					  )
  const {
  int return_value=false_value;
  double norm = std::inner_product(vector_data.begin(),vector_data.end(),vector_data.begin(),0.0);
  if(norm<0.00000000001)
    {
      return_value=true_value;
    };
  return return_value;
};
//
inline double local_state_class::normalized_product(
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
//
void local_state_class::get_cpu_time(
				     clock_t & cumulative_time,
				     const std::string & job_flag 
				     )
{
  if(job_flag=="init")
    {
      cumulative_time_local=clock()-clock();
    }
  else if(job_flag=="start")
    {
      sub_start_time=clock();
    }
  else if(job_flag=="add")
    {
      sub_end_time=clock();
      cumulative_time_local=cumulative_time_local+(sub_end_time-sub_start_time);
    }
  else if(job_flag=="get")
    {
      cumulative_time=cumulative_time_local;
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "local state class",
			     "get_cpu_time",
			     "due to irrigal job_flag."
			     );
    };
};		
//
double local_state_class::get_volume_energy_difference(
						       const std::vector<long int> & configuration
						       )
const {
  double return_value=0.0;
//    fprintf(stderr,"volume:%ld,%ld\n", volume,natural_volumes[cell_type]);
  if((configuration[local_site]==cell)&&(cell_type!=buffer_type))
    {
        return_value = -balk_moduluses[cell_type]*(2.0*(double)(volume-natural_volumes[cell_type])-1);
    }
  else if(cell_type!=buffer_type)
    {
        return_value = balk_moduluses[cell_type]*(2.0*(double)(volume-natural_volumes[cell_type])+1);
    };
  return return_value;
};
//
double state_system_class::get_energy_difference() 
  const {
//  	fprintf(stderr,"%f,%f,%f,%f\n",candidate_adhesion_energy,present_adhesion_energy,present_volume_energy_difference,candidate_volume_energy_difference);
  return 
      candidate_adhesion_energy
    - present_adhesion_energy
    + present_volume_energy_difference
    + candidate_volume_energy_difference
    + polarity_driving_difference
    + field_driving_difference ;
};
//
void state_system_class::get_total_adhesion(
					    const model_parameters_cellular_potts_class & model,
					    const cell_system_class & cell_system,
					    site_system_class & site_system,
					    double & total_adhesion_energy_for_cell,
					    double & total_intercell_contact_for_cell,
					    double & total_intercell_head_contact_for_cell,
					    double & total_intercell_tail_contact_for_cell,
					    double & total_intercell_vartical_contact_for_cell,
					    double & total_intercell_lateral_contact_for_cell,
					    double & total_intercell_ordered_contact_for_cell,
					    double & total_cellmedium_contact_for_cell,
					    double & total_square_polar_product
					    )
{
  long long int site_index;
  long long int neighbor_site;
  int neighbor_index;
  total_adhesion_energy_for_cell=0.0;
  total_intercell_contact_for_cell=0.0;
  total_intercell_head_contact_for_cell=0.0;
  total_intercell_tail_contact_for_cell=0.0;
  total_intercell_lateral_contact_for_cell=0.0;
  total_intercell_ordered_contact_for_cell=0.0;
  total_cellmedium_contact_for_cell=0.0;
  total_square_polar_product=0.0;
  //
  gen_polar_product(
		    model,
		    site_system,
		    polar_product
		    );
  //
  //  int debug_flag=0;
  //
  for(site_index=0;
      site_index<number_of_sites;
      site_index++)
    {
      if(configuration[site_index]!=buffer_cell)
	{
	  if(cell_volumes[configuration[site_index]]>0)
	    {
	      total_square_polar_product=total_square_polar_product
		+polar_product[site_index]*polar_product[site_index]
		/cell_volumes[configuration[site_index]];
	    };
	  for(neighbor_index=0;
	      neighbor_index<number_of_neighbor_sites;
	      neighbor_index++)
	    {
	      neighbor_site=site_system.get_neighbor_site(site_index,neighbor_index);
	      if(
		 configuration[neighbor_site]!=configuration[site_index]
		 )
		{
		  if(
		     configuration[neighbor_site]!=buffer_cell
		     &&
		     neighbor_site<site_index
		     )
		    {
		      total_adhesion_energy_for_cell
			+=isotropic_adhesion_energy(
						    cell_types[configuration[site_index]],
						    cell_types[configuration[neighbor_site]]
						    );
		      
		      total_adhesion_energy_for_cell
			+=dipolar_adhesion_energy(
						  cell_types[configuration[site_index]],
						  polar_product[site_index]
						  +dipolar_adhesion_basals[cell_types[configuration[site_index]]],
						  cell_types[configuration[neighbor_site]],
						  polar_product[neighbor_site]
						  +dipolar_adhesion_basals[cell_types[configuration[neighbor_site]]]
						  );
		      total_intercell_contact_for_cell+=1.0;
		      total_intercell_head_contact_for_cell+=(1.0+polar_product[site_index])/2.0;
		      total_intercell_tail_contact_for_cell+=(1.0-polar_product[site_index])/2.0;
		      total_intercell_vartical_contact_for_cell+=polar_product[site_index]*polar_product[site_index]*polar_product[neighbor_site]*polar_product[neighbor_site];
		      total_intercell_lateral_contact_for_cell+=(1.0-polar_product[site_index]*polar_product[site_index])*(1.0-polar_product[neighbor_site]*polar_product[neighbor_site]);
		      total_intercell_ordered_contact_for_cell+=polar_product[site_index]*polar_product[neighbor_index];
		    }
		  else if(configuration[neighbor_site]==buffer_cell)
		    {
		      total_adhesion_energy_for_cell
			+=isotropic_adhesion_energy(
						    cell_types[configuration[site_index]],
						    cell_types[configuration[neighbor_site]]
						    );
		      total_cellmedium_contact_for_cell+=1.0;
		    };
		}
	      else
		{
		};
	    };
	};
    };
  total_adhesion_energy_for_cell=total_adhesion_energy_for_cell/(double)number_of_cells;
  total_intercell_contact_for_cell=2.0*total_intercell_contact_for_cell/(double)number_of_cells; 
  total_intercell_head_contact_for_cell=2.0*total_intercell_head_contact_for_cell/(double)number_of_cells; 
  total_intercell_tail_contact_for_cell=2.0*total_intercell_tail_contact_for_cell/(double)number_of_cells; 
  total_intercell_vartical_contact_for_cell=2.0*total_intercell_vartical_contact_for_cell/(double)number_of_cells; 
  total_intercell_lateral_contact_for_cell=2.0*total_intercell_lateral_contact_for_cell/(double)number_of_cells; 
  total_intercell_ordered_contact_for_cell=2.0*total_intercell_ordered_contact_for_cell/(double)number_of_cells; 
  total_cellmedium_contact_for_cell=total_cellmedium_contact_for_cell/(double)number_of_cells;
  total_square_polar_product=total_square_polar_product/(double)number_of_cells;
};
//
void state_system_class::get_adhesion_field(
					    const model_parameters_cellular_potts_class & model,
					    const cell_system_class & cell_system,
					    site_system_class & site_system,
					    const adhesion_system_class & adhesion_system
					    )
{
  int cell_type;
  int neighbor_type;
  long long int neighbor_site;
  int pair_index;
  double adhesion_strength;
  double adhesion_product;
  double adhesion_density;
  double work_norm;
  //debug
  ///*
  std::string debug_message;
  //io_method.standard_output("enter:");
  //*/
  //debug
  gen_polar_product(
		    model,
		    site_system,
		    polar_product
		    );
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
	  adhesion_field[cell_index*space_dimension+direction_index]=0.0;
	};
    };
  for(
      long long int site_index=0;
      site_index<number_of_sites;
      site_index++
      )
    {
      if(configuration[site_index]!=buffer_cell)
	{
	  // debug
	  //if(configuration[site_index]==90) message="pass 0:";
	  //
	  cell_type=cell_system.get_type(configuration[site_index]);
	  calculate_relative_coordinate(
					model,
					site_system,
					site_index,
					work_relative_coordinate
					);
	  work_norm=norm(work_relative_coordinate);
	  if(tool.finite_abs_check(work_norm)==tool.get_value("true"))
	    {
	      for(
		  int neighbor_index=0;
		  neighbor_index<number_of_neighbor_sites;
		  neighbor_index++
		  )
		{
		  neighbor_site=site_system.get_neighbor_site(
							      site_index,
							      neighbor_index
							      );
		  neighbor_type=cell_system.get_type(configuration[neighbor_site]);
		  pair_index=cell_type
		    *adhesion_system.number_of_cell_types
		    +neighbor_type;
		  if(
		     adhesion_system.typepair_to_adhesion_flags[pair_index]
		     &&
		     neighbor_type!=buffer_type
		     &&
		     configuration[site_index]!=configuration[neighbor_site]
		     )
		    {
		      // debug
		      ///*
		      //if(configuration[site_index]==90) message="pass 1:" 
		      //				      +io_method.int_to_string(neighbor_index);
		      //*/
		      // debug
		      for(
			  int map_index=adhesion_system.typepair_to_adhesion_pointers[pair_index];
			  map_index<adhesion_system.typepair_to_adhesion_pointers[pair_index+1];
			  map_index++
			  )
			{
			  // debug
			  /*
			  if(configuration[site_index]==90) message="pass 2:" 
							      +io_method.int_to_string(map_index);
			  */
			  // debug
			  // should be debugging
			  // debug
			  //	  debug_message= adhesion_system.get_adhesion_type(map_index)
			  //	    + ":";
			  // debug
			  if(
			     adhesion_system.get_adhesion_type(map_index)
			     ==
			     adhesion_system.get_adhesion_type_identifier("normal")
			     )
			    {
			      adhesion_strength
				=adhesion_system.coupling_constants[map_index];
			      // debug
			      //debug_message += "debug list A: strength =" 
			      //	+ io_method.double_to_string(adhesion_strength)
			      //	+ "map no." + io_method.longint_to_string(map_index)
			      //	+ "cell no." + io_method.longint_to_string(configuration[site_index])
			      //	+ "neig no." + io_method.longint_to_string(configuration[neighbor_site]);
			      //  io_method.standard_output(debug_message);
			      // debug
			    }
			  else if (
				   adhesion_system.get_adhesion_type(map_index)
				   ==
				   adhesion_system.get_adhesion_type_identifier("tight")
				   )
			    {
			      adhesion_strength
				=adhesion_system.get_adhesion_tight_junction_constant(
										      map_index,
										      configuration[site_index],
										      configuration[neighbor_site]
										      ) ;
			      // debug
			      //     debug_message += "debug list B: strength =" 
			      //	+ io_method.double_to_string(adhesion_strength)
			      //	+ "map no." + io_method.longint_to_string(map_index)
			      //	+ "cell no." + io_method.longint_to_string(configuration[site_index])
			      //	+ "neig no." + io_method.longint_to_string(configuration[neighbor_site]);
			      //    io_method.standard_output(debug_message);
			      // debug
				};
			  // debug
			  //	  if(configuration[site_index]==135) 
			  //{
			  //  message="debug 2:" 
			  //	+io_method.double_to_string(adhesion_system.coupling_constants[map_index]);
			  //  io_method.standard_output(message);
			  //};
			  // debug
			  adhesion_product=1.0;
			  adhesion_density=0.0;
			  for(
			      int component_index=1;
			      component_index<adhesion_system.number_of_adhesion_components;
			      component_index++
			      )
			    {
			      adhesion_density
				+=((double)component_index)
				*adhesion_product
				*adhesion_system.polar_components_1[
								    map_index
								    *adhesion_system.number_of_adhesion_components
								    +component_index
								    ];
			      adhesion_product
				=adhesion_product
				*polar_product[site_index];
			    };
			  adhesion_strength
			    =adhesion_strength*adhesion_density;
			  // debug
			  /*
			  if(configuration[site_index]==90) 
			    message+="pass 3a:"+io_method.double_to_string(adhesion_strength)+',';
			  */
		      // debug
			  adhesion_product=1.0;
			  adhesion_density=0.0;
			  for(
			      int component_index=0;
			      component_index<adhesion_system.number_of_adhesion_components;
			      component_index++
			      )
			    {
			      adhesion_density
			    +=adhesion_product
				*adhesion_system.polar_components_2[
								    map_index
								    *adhesion_system.number_of_adhesion_components
								    +component_index
								    ];
			      adhesion_product
				=adhesion_product
				*polar_product[neighbor_index];
			      // debug
			      /*
			      if(configuration[site_index]==90) 
				{
				  message+="pass 3b01:"+io_method.int_to_string(component_index)+',';
				  message+="pass 3b02:"+io_method.double_to_string(adhesion_system.polar_components_2[
														      map_index
														      *adhesion_system.number_of_adhesion_components
														      +component_index
														      ])+',';
				};
			      */
			      // debug
			    };
			  adhesion_strength
			    =adhesion_strength*adhesion_density;
			  // debug
			  /*
			  if(configuration[site_index]==90)
			    {
			      message+="pass 3b:"+io_method.double_to_string(adhesion_strength)+',';
			      message+="pass 3b2:"+io_method.double_to_string(adhesion_product)+',';
			    };
			  */
			  //debug
			  // debug
			  //    if(
			  //	   std::isnan(
			  //		      adhesion_density
			  //		      )
			  //	   )
			  //	  {
			  //	    std::string message;
			  //	    message =io_method.double_to_string(adhesion_strength);
			  //	    message+=io_method.double_to_string(adhesion_product);
			  //	    io_method.error_output("debug 1","",message);
			  //	  };
			  //
			  for(
			      int direction_index=0;
			      direction_index<space_dimension;
			      direction_index++
			      )
			    {
			      adhesion_field[
					     configuration[site_index]
					     *space_dimension
					     +direction_index
					     ]
				+=adhesion_strength
				*work_relative_coordinate[direction_index]
				/work_norm;
			      // debug
			      /*
			      if(configuration[site_index]==90)
				{
				  message+="pass 3b3:"+io_method.double_to_string(adhesion_strength)+',';
				  message+="pass 3b4:"+io_method.double_to_string(work_norm)+',';
				  message+="pass 3b5:"+io_method.double_to_string(work_relative_coordinate[direction_index])+',';
				  message+="pass 3b6:"+io_method.double_to_string(adhesion_field[
												 configuration[site_index]
												 *space_dimension
												 +direction_index
												 ])+',';
				}
			      */
			      // debug
			      // debug
			      //if(configuration[site_index]==50&&map_index==0)
			      //  {
			      //    io_method.standard_output("pass 3b:"+io_method.double_to_string(adhesion_strength));
			      //  };
			      // debug
			    };
			};
		    };
		};
	    };
	  // debug
	  /*
	   if(configuration[site_index]==90) 
	     io_method.standard_output(message);
	  */
	  // debug
	};
    };
  // debug
  /*
  message="";
  for(
      int direction_index=0;
      direction_index<space_dimension;
      direction_index++
      )
    {
      message+=io_method.double_to_string(
	adhesion_field[
		       50
		       *space_dimension
		       +direction_index
		       ])+',';
    };
  io_method.standard_output(message);
  io_method.standard_output("exit:");
  */
  //debug
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
	  adhesion_field[cell_index*space_dimension+direction_index]
	    =-adhesion_field[cell_index*space_dimension+direction_index];
	};
      /*
      if(cell_index==90)
	{
	  message="AF:";
	    for(
		int direction_index=0;
		direction_index<space_dimension;
		direction_index++
		)
	      {
		message+=io_method.double_to_string(adhesion_field[cell_index*space_dimension+direction_index])+',';
	      };
	  io_method.standard_output(message);
	};
      */
    };
};
//
inline double state_system_class::norm(
				       const std::vector<double> & input_vector
				       )
{
  return sqrt(
	      std::inner_product(
				 input_vector.begin(),
				 input_vector.end(),
				 input_vector.begin(),
				 0.0
				 )
	      );
};
