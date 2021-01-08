#include "cellular_potts_hamiltonian.hpp"
#include "cellular_potts_state.hpp"
//
double local_state_class::get_local_adhesion_energy(
						    const adhesion_system_class & adhesion_system
						    )
{
  double return_value=0.0;
  double work_density=0.0;
  double work_density_neighbor=0.0;
  double work_polarity_value=1.0;
  double polar_product_neighbor=0.0;
  int pair_index;
  int neighbor_index;
  std::string debug_message;
  if(
     (cell_type!=buffer_type)
     &&
     (iszero_norm_in_adhesive_hamiltonian(relative_coordinates)==false_value)
     )
    {
      for(neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
	{
	  pair_index
	    =cell_type
	    *adhesion_system.number_of_cell_types
	    +neighbor_types[neighbor_index];
	  work_vector_double_for_neighbors[neighbor_index]=0.0;
	  if(
	     (
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
	      (iszero_norm_in_adhesive_hamiltonian(neighbor_relative_coordinates[neighbor_index])
	       ==
	       true_value
	       )
	      )
	     )
	    {
	    }
	  else if(adhesion_system.typepair_to_adhesion_flags[pair_index])
	    {
	      for(
		  int map_index
		    = adhesion_system.typepair_to_adhesion_pointers[pair_index];
		  map_index
		    < adhesion_system.typepair_to_adhesion_pointers[pair_index+1];
		  map_index++
		  )
		{
		  work_polarity_value=1.0;
		  work_density=0.0;
		  for(
		      int component_index=0;
		      component_index<adhesion_system.number_of_adhesion_components;
		      component_index++
		      )
		    {
		      work_density 
			+= work_polarity_value
			*adhesion_system.polar_components_1[
							    map_index*adhesion_system.number_of_adhesion_components
							    +component_index
							    ];
		      work_polarity_value
			=work_polarity_value
			*product_polarity;
		    };
		  polar_product_neighbor
		    =normalized_product_in_adhesive_hamiltonian(
								neighbor_polarities[neighbor_index],
								neighbor_relative_coordinates[neighbor_index]
								);
		  work_polarity_value=1.0;
		  work_density_neighbor=0.0;
		  for(
		      int component_index=0;
		      component_index<adhesion_system.number_of_adhesion_components;
		      component_index++
		      )
		    {
		      work_density_neighbor
			+=work_polarity_value
			*adhesion_system.polar_components_2[
							    map_index
							    *adhesion_system.number_of_adhesion_components
							    +component_index
							    ];    
		      work_polarity_value
			=work_polarity_value
			*polar_product_neighbor;
		    };
		  // should be debugging
		  //		  debug_message = adhesion_system.get_adhesion_type(map_index)
		  //  + ":";
		  if(
		     adhesion_system.get_adhesion_type(map_index)
		     ==
		     adhesion_system.get_adhesion_type_identifier("normal")
		     )
		    {
		      work_vector_double_for_neighbors[neighbor_index]
			+=adhesion_system.coupling_constants[map_index]
			*work_density*work_density_neighbor;
		      //debug
		      //  debug_message += "debug list A: strength =" 
		      // 	+ io_method.double_to_string(adhesion_system.coupling_constants[map_index])
		      // 	+ "map no." + io_method.longint_to_string(map_index)
		      //	+ "cell no." + io_method.longint_to_string(cell)
		      //	+ "neig no." + io_method.longint_to_string(neighbor_cells[neighbor_index]);
		      //io_method.standard_output(debug_message);
		      //debug
		    }
		  else if(
			  adhesion_system.get_adhesion_type(map_index)
			  ==
			  adhesion_system.get_adhesion_type_identifier("tight")
			  )
		    {
		      work_vector_double_for_neighbors[neighbor_index]
			+=adhesion_system.get_adhesion_tight_junction_constant(
									      map_index,
									      cell,
									      neighbor_cells[neighbor_index]
									      )
			*work_density*work_density_neighbor;
		      //debug
		      //      debug_message += "debug list B: strength =" 
		      //	+ io_method.double_to_string(
		      //				     adhesion_system.get_adhesion_tight_junction_constant(
		      //											  map_index,
		      //											  cell,
		      //											  neighbor_cells[neighbor_index]
		      //										  )
		      //				     )
			//	+ "map no." + io_method.longint_to_string(map_index)
			//	+ "cell no." + io_method.longint_to_string(cell)
			//	+ "neig no." + io_method.longint_to_string(neighbor_cells[neighbor_index]);
		      //		      io_method.standard_output(debug_message);
		      //debug
		    };
		  //
		};
	    };
	  // ;
	};
      return_value
	=accumulate(
		    work_vector_double_for_neighbors.begin(),
		    work_vector_double_for_neighbors.end(),
		    0.0
		    );
    };
  // debug
  /*
  std::string messages;
  messages = "debug:";
  messages+= "cell_type =" + io_method.int_to_string(cell_type) + '_';
  io_method.standard_output(messages);
  for( int neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
    {
      messages = "energy["+ io_method.int_to_string(neighbor_index)+ ","
	+  io_method.int_to_string(neighbor_types[neighbor_index]) + "]"
	+ "=" + io_method.double_to_string(work_vector_double_for_neighbors[neighbor_index]);
      io_method.standard_output(messages);
    };
  //
  */
  return return_value;
};
//
inline int local_state_class::iszero_norm_in_adhesive_hamiltonian(
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
inline double local_state_class::normalized_product_in_adhesive_hamiltonian(
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
