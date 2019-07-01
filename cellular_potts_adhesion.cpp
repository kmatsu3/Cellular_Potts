#include "cellular_potts_adhesion.hpp"
#include <cstdlib>
  /*======================
    Constructor
   =======================*/
adhesion_binding_table::adhesion_binding_table()
{
  bind_partner_index.clear();
  bind_partner_pointer.clear();
};
/*======================
    Method
 =======================*/
void adhesion_binding_table::set_table(
				       const int & input_adhesion_index,
				       const long int & number_of_cells,
				       const double & coupling_constant
				       )
{
  adhesion_index=input_adhesion_index;
  std::string structure_item;
  std::vector<long int> work_vector_longint;
  bind_coupling.clear();
  for(
      long int cell_index=0;
      cell_index<number_of_cells;
      cell_index++
      )
    {
      for(
	  long int partner_index=0;
	  partner_index < number_of_cells;
	  partner_index++
	  )
	{
	  bind_coupling.push_back(0.0);
	}
    };
  bind_coupling.push_back(0.0);
  bind_partner_pointer.push_back(0);
  for(
      long int cell_index=0;
      cell_index<number_of_cells;
      cell_index++
      )
    {
      structure_item=set_structure_bind_table_item(
						   adhesion_index,
						   cell_index,
						   "bind_table_components"
						   );
      work_vector_longint.clear();
      io_method.get_input_longint_array(
					"bind_table",
					structure_item,
					work_vector_longint
					);
      bind_partner_pointer.push_back((long int)work_vector_longint.size());
      for(
	  long int partner_index=work_vector_longint[cell_index];
	  partner_index<work_vector_longint[cell_index+1];
	  partner_index++
	  )
	{
	  bind_partner_index.push_back(work_vector_longint[partner_index]);
	  bind_coupling[cell_index*number_of_cells+partner_index]
	    =coupling_constant;
	} 
    };
  show_binding_table(number_of_cells);
};
//
void adhesion_binding_table::show_binding_table(
						const long int & number_of_cells
						)
{
  io_method.standard_output("/// bind table No." 
			    + io_method.int_to_string(adhesion_index)
			    + " ///");
  std::string message;
  for(long int cell_index=0;
      cell_index<number_of_cells;
      cell_index++)
    {
      message = "  Cell [" + io_method.longint_to_string(cell_index)
	+ "]";
      io_method.standard_output(message);
      for(long int partner_index=0;
	  partner_index<number_of_cells;
	  partner_index++)
	{
	  if(tools.finite_abs_check(
			      bind_coupling[
					    cell_index*number_of_cells
					    +partner_index
					    ]
			      )
	     )
	    {
	      message = "No." + io_method.longint_to_string(partner_index)
		+ ": J = " +io_method.double_to_string(
						       bind_coupling[
								     cell_index*number_of_cells
								     +partner_index]);
	      io_method.standard_output(message);
	    };
	};
    };
};
//
std::string adhesion_binding_table::set_structure_bind_table_item(
								  const int & adhesion_index,
								  const long int & cell_index,
								  const std::string & child
								  )
{
  return io_method.generate_structure(
				      "bind_table["+io_method.int_to_string(adhesion_index)+"].cell",
				      cell_index,
				      child
				      );
};
  /*======================
    Constructor
   =======================*/
adhesion_cellular_potts_class::adhesion_cellular_potts_class()
  {
    coupling_constant=0.0;
    cell_type_1=-1;
    polarity_type_1=-1;
    adhesion_polarity_1.clear();
    cell_type_2=-1;
    polarity_type_2=-1;
    adhesion_polarity_2.clear();
    interaction_type="empty";
  };
/*======================
    Method
 =======================*/

void adhesion_cellular_potts_class::set_adhesion_int_value(
							   const std::string & type,
							   const int & value
							   )
{
  if(type=="cell_type_1")
    {
      cell_type_1=value;
    }
  else if(type=="polarity_type_1")
    {
      polarity_type_1=value;
    }
  else if(type=="cell_type_2")
    {
      cell_type_2=value;
    }
  else if(type=="polarity_type_2")
    {
      polarity_type_2=value;
    }
  else
    {
      io_cellular_potts io_method;
      std::string message;
      message = "In input, ";
      message+= type;
      message+= " is undefined(adhesion_cellular_potts_class::set_adhesion_int_value).";
      io_method.standard_output(message);
      std::abort();
    };
};
//
void adhesion_cellular_potts_class::set_adhesion_double_value(
							      const std::string & type,
							      const double & value
							      )
{
  if(type=="coupling_constant")
    {
      coupling_constant=value;
    }
  else
    {
      io_cellular_potts io_method;
      std::string message;
      message = "In input, ";
      message+= type;
      message+= " is undefined(adhesion_cellular_potts_class::set_adhesion_double_value).";
      io_method.standard_output(message);
      std::abort();
    };
};
//
void adhesion_cellular_potts_class::set_adhesion_double_vector(
							       const int & input_cell_index,
							       const std::vector<double> & adhesion_polarity
							       )
{
  if(input_cell_index==1)
    {
      std::vector<double>::const_iterator index=adhesion_polarity.begin();
      adhesion_polarity_1.clear();
      while(index!=adhesion_polarity.end())
	{
	  adhesion_polarity_1.push_back(*index);
	  index++;
	};
    }
  else if(input_cell_index==2)
    {
      std::vector<double>::const_iterator index=adhesion_polarity.begin();
      adhesion_polarity_2.clear();
      while(index!=adhesion_polarity.end())
	{
	  adhesion_polarity_2.push_back(*index);
	  index++;
	};
    }
  else
    {
      io_cellular_potts io_method;
      std::string message;
      message = "In input, cell type ";
      message+= io_method.int_to_string(input_cell_index);
      message+= " is undefined(adhesion_cellular_potts_class::set_adhesion_double_vector).";
      io_method.standard_output(message);
      std::abort();
    };
};
//
void adhesion_cellular_potts_class::set_adhesion_string(
							const std::string & type,
							const std::string & value
							      )
{
  if(type=="interaction_type")
    {
      interaction_type=value;
    }
  else
    {
      io_cellular_potts io_method;
      std::string message;
      message = "In input, ";
      message+= type;
      message+= " is undefined(adhesion_cellular_potts_class::set_adhesion_double_value).";
      io_method.standard_output(message);
      std::abort();
    };
};
//
double adhesion_cellular_potts_class::get_coupling_constant() const 
{
  return coupling_constant;
};
//
int adhesion_cellular_potts_class::get_adhesion_cell_type(
							  const int & input_cell_index
							  ) const 
{
  int return_value;
  if(input_cell_index==1)
    {
      return_value=cell_type_1;
    }
  else if(input_cell_index==2)
    {
      return_value=cell_type_2;
    }
  else
    {
      io_cellular_potts io_method;
      std::string message;
      message = "In input, cell type";
      message+= '[' + io_method.int_to_string(input_cell_index) + ']';
      message+= " is undefined(adhesion_cellular_potts_class::get_adhesion_cell).";
      io_method.standard_output(message);
      std::abort();
    };
  return return_value;
};
//
int adhesion_cellular_potts_class::get_adhesion_polarity_type(
							      const int & input_cell_index
							      ) const 
{
  int return_value;
  if(input_cell_index==1)
    {
      return_value=polarity_type_1;
    }
  else if(input_cell_index==2)
    {
      return_value=polarity_type_2;
    }
  else
    {
      io_cellular_potts io_method;
      std::string message;
      message = "In input, cell type";
      message+= '[' + io_method.int_to_string(input_cell_index) + ']';
      message+= " is undefined(adhesion_cellular_potts_class::get_adhesion_cell).";
      io_method.standard_output(message);
      std::abort();
    };
  return return_value;
};
//
void adhesion_cellular_potts_class::get_adhesion_polarity(
							  const int & input_cell_index,
							  std::vector<double> & polarity
							  ) const 
{
  if((int)polarity.size()!=(int)adhesion_polarity_1.size())
    {
      io_cellular_potts io_method;
      std::string message;
      message = "In input, ";
      message+= "polarity[" + io_method.int_to_string(input_cell_index) + "]";
      message+= "does not have the compatible size (adhesion_cellular_potts_class::get_adhesion_polarity).";
      io_method.standard_output(message);
      std::abort();
    };
  if(input_cell_index==1)
    {
      std::vector<double>::const_iterator index=adhesion_polarity_1.begin();
      int counter=0;
      while(index!=adhesion_polarity_1.end())
	{
	  polarity[counter]=*index;
	  index++;
	  counter++;
	};
    }
  else if(input_cell_index==2)
    {
      std::vector<double>::const_iterator index=adhesion_polarity_2.begin();
      int counter=0;
      while(index!=adhesion_polarity_2.end())
	{
	  polarity[counter]=*index;
	  index++;
	  counter++;
	};
    }
  else
    {
      io_cellular_potts io_method;
      std::string message;
      message = "In input, ";
      message+= "input type = " + io_method.int_to_string(input_cell_index) + "";
      message+= " is undefined(adhesion_cellular_potts_class::get_adhesion_polarity).";
      io_method.standard_output(message);
      std::abort();
    };
};
//
double adhesion_cellular_potts_class::get_adhesion_cell1_component(
								    const int & component_index
								    ) const 
{
  return adhesion_polarity_1[component_index];
};
//
double adhesion_cellular_potts_class::get_adhesion_cell2_component(
								    const int & component_index
								    ) const 
{
  return adhesion_polarity_2[component_index];
};
//
std::string adhesion_cellular_potts_class::get_interaction_type() const 
{
  return interaction_type;
};
//
void adhesion_system_class::initialize_adhesion()
{
  adhesion_cellular_potts_class work_adhesion;
  int work_int;
  double work_double;
  std::string work_string;
  std::vector<double> work_vector_double;
  adhesions.clear();
  std::string structure_item;
  for(
      int adhesion_index=0;
      adhesion_index<number_of_adhesion;
      adhesion_index++
      )
    {
      // cell_type_1
      structure_item=set_structure_adhesion_item(
						 adhesion_index,
						 "cell[1].cell_type"
						 );
      work_int=io_method.get_input_int(
				       "adhesion_input",
				       structure_item
				       );
      work_adhesion.set_adhesion_int_value(
					   "cell_type_1",
					   work_int
					   );
      // polarity_type_1
      structure_item=set_structure_adhesion_item(
						 adhesion_index,
						 "cell[1].polarity_type"
						 );
      work_int=io_method.get_input_int(
				       "adhesion_input",
				       structure_item
				       );
      work_adhesion.set_adhesion_int_value(
					   "polarity_type_1",
					   work_int
					   );
      //cell 1 polarity
      structure_item=set_structure_adhesion_item(
						 adhesion_index,
						 "cell[1].adhesion_components"
						 );
      work_vector_double.clear();
      io_method.get_input_double_array(
				       "adhesion_input",
				       structure_item,
				       work_vector_double
				       );
      work_adhesion.set_adhesion_double_vector(
					       1,
					       work_vector_double
					       );
      // cell_type_2
      structure_item=set_structure_adhesion_item(
						 adhesion_index,
						 "cell[2].cell_type"
						 );
      work_int=io_method.get_input_int(
				       "adhesion_input",
				       structure_item
				       );
      work_adhesion.set_adhesion_int_value(
					   "cell_type_2",
					   work_int
					   );
      // polarity_type_2
      structure_item=set_structure_adhesion_item(
						 adhesion_index,
						 "cell[2].polarity_type"
						 );
      work_int=io_method.get_input_int(
				       "adhesion_input",
				       structure_item
				       );
      work_adhesion.set_adhesion_int_value(
					   "polarity_type_2",
					   work_int
					   );
      //cell 2 polarity
      structure_item=set_structure_adhesion_item(
						 adhesion_index,
						 "cell[2].adhesion_components"
						 );
      work_vector_double.clear();
      io_method.get_input_double_array(
				       "adhesion_input",
				       structure_item,
				       work_vector_double
				       );
      work_adhesion.set_adhesion_double_vector(
					       2,
					       work_vector_double
					       );
      //coupling_constant
      structure_item=set_structure_adhesion_item(
						 adhesion_index,
						 "coupling_constant"
						 );
      work_double=io_method.get_input_double(
					     "adhesion_input",
					     structure_item
					     );
      work_adhesion.set_adhesion_double_value(
					      "coupling_constant",
					      work_double
					      );
      //interaction_type
      structure_item=set_structure_adhesion_item(
						 adhesion_index,
						 "interaction_type"
						 );
      work_double=io_method.get_input_double(
					     "adhesion_input",
					     structure_item
					     );
      interaction_key_check(work_string);
      work_adhesion.set_adhesion_string(
					"interaction_type",
					work_string
					);
      // generation table
      if(work_adhesion.get_interaction_type()==adhesion_type_tight)
	work_adhesion.bind_table.set_table(
					   adhesion_index,
					   number_of_cells,
					   work_adhesion.get_coupling_constant()
					   );
	//
      adhesions.push_back(work_adhesion);
    };
  make_typepair_to_adhesion_maps();
};
//
 void adhesion_system_class::show_adhesion()
 {
   std::string messages;
   std::vector<double> work_vector(number_of_adhesion_components,0.0);
   messages= "==== Adhesion data started ====";
   io_method.standard_output(messages);
   for(
       int adhesion_index=0;
       adhesion_index<number_of_adhesion;
       adhesion_index++
       )
     {
       messages = "****adhesion [" 
	 + io_method.int_to_string(adhesion_index) + "]"
	 + "****";
       io_method.standard_output(messages);
       messages = "+ cell_type[1] = " 
	 + io_method.int_to_string(adhesions[adhesion_index].get_adhesion_cell_type(1)) + "]"
	 + "****";
       io_method.standard_output(messages);
       messages = "+ cell polarity[1] = " 
	 + io_method.int_to_string(adhesions[adhesion_index].get_adhesion_polarity_type(1)) + "]";
       io_method.standard_output(messages);
       //
       adhesions[adhesion_index].get_adhesion_polarity(
						       1,
						       work_vector
						       );
       messages = "+ cell[1] adhesion componet = (";
       for(
	   int component_index=0;
	   component_index<(int)work_vector.size();
	   component_index++
	   )
	 {
	   messages += io_method.double_to_string(work_vector[component_index]);
	   if(component_index<(int)work_vector.size()-1)
	     {messages += ", ";}
	 };
       messages += ")";
       io_method.standard_output(messages);
       //
       messages = "****cell type [2] = " 
	 + io_method.int_to_string(adhesions[adhesion_index].get_adhesion_cell_type(2)) + "]"
	 + "****";
       io_method.standard_output(messages);
       messages = "****cell polarity[2] = " 
	 + io_method.int_to_string(adhesions[adhesion_index].get_adhesion_polarity_type(2)) + "]";
       io_method.standard_output(messages);
       //
       adhesions[adhesion_index].get_adhesion_polarity(
						       2,
						       work_vector
						       );
       messages = "+ cell[2] adhesion componet = (";
       for(
	   int component_index=0;
	   component_index<(int)work_vector.size();
	   component_index++
	   )
	 {
	   messages += io_method.double_to_string(work_vector[component_index]);
	   if(component_index<(int)work_vector.size()-1)
	     {messages += ", ";}
	 };
       messages += ")";
       io_method.standard_output(messages);
       //
       messages = "**** coupling constant  = " 
	 + io_method.double_to_string(
				      adhesions[adhesion_index].get_coupling_constant()
				      ) ;
       io_method.standard_output(messages);
     };
   show_map_table();
   messages= "==== Adhesion data finished ====";
 };
//
void adhesion_system_class::show_map_table()
{
  int pair_index;
  int map_index;
  std::string messages;
  messages = "**** map table: ";
    io_method.standard_output(messages);
  for(
      int type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      //      io_method.standard_output("type=" + io_method.int_to_string(type_index));
      for(
	  int pair_type_index=0;
	  pair_type_index<number_of_cell_types;
	  pair_type_index++
	  )
	{
	  pair_index=type_index*number_of_cell_types+pair_type_index;
	  //	  io_method.standard_output("pair_index="+ io_method.int_to_string(pair_index));
	  //	  io_method.standard_output("pair_type="+ io_method.int_to_string(pair_type_index));
	  if(typepair_to_adhesion_flags[pair_index])
	    {
	      for(
		  map_index=typepair_to_adhesion_pointers[pair_index];
		  map_index<typepair_to_adhesion_pointers[pair_index+1];
		  map_index++
		  )
		{
		  //		  io_method.standard_output("map_index="+ io_method.int_to_string(map_index));
		  messages = "+ map [" 
		    + io_method.int_to_string(type_index) + ","
		    + io_method.int_to_string(pair_type_index) + "] ="
		    + io_method.int_to_string(
					      typepair_to_adhesion_maps[map_index]
					      );
		  io_method.standard_output(messages);
		  messages = "* coupling const = "
		    + io_method.double_to_string(coupling_constants[map_index]);
		  io_method.standard_output(messages);
		  messages = "* cell 1 component = (";
		  for(
		      int component_index=0;
		      component_index<number_of_adhesion_components;
		      component_index++
		      )
		    {
		      messages += io_method.double_to_string(
							     polar_components_1[
										map_index*number_of_adhesion_components+component_index
										]
							     ) ;
		      if(component_index<number_of_adhesion_components-1) messages += ',';
		    };
		  messages += ")";
		  io_method.standard_output(messages);
		  messages = "* cell 2 component = (";
		  for(
		      int component_index=0;
		      component_index<number_of_adhesion_components;
		      component_index++
		      )
		    {
		      messages += io_method.double_to_string(
							     polar_components_2[
										map_index*number_of_adhesion_components+component_index
										]
							     ) ;
		      if(component_index<number_of_adhesion_components-1) messages += ',';
		    };
		  messages += ")";
		  io_method.standard_output(messages);
		};
	    };
	};
    };
};
//
int adhesion_system_class::get_adhesion_cell_type(
						  const int & adhesion_index,
						  const int & cell_index
						  ) const 
{
  return adhesions[adhesion_index].get_adhesion_cell_type(
							  cell_index
							  );
};
//
int adhesion_system_class::get_adhesion_polarity_type(
						      const int & adhesion_index,
						      const int & cell_index
						      ) const 
{
  return adhesions[adhesion_index].get_adhesion_polarity_type(
							      cell_index
							      );
};
//
void adhesion_system_class::get_adhesion_polarity(
						  const int & adhesion_index,
						  const int & cell_index,
						  std::vector<double> & polarity
						  ) const 
{
  adhesions[adhesion_index].get_adhesion_polarity(
						  cell_index,
						  polarity
						  );
};
//
double  adhesion_system_class::get_adhesion_cell1_component(
							    const int & adhesion_index,
							    const int & component_index
							    ) const
{
  return adhesions[adhesion_index].get_adhesion_cell1_component(component_index);
};
//
double  adhesion_system_class::get_adhesion_cell2_component(
							    const int & adhesion_index,
							    const int & component_index
							    ) const
{
  return adhesions[adhesion_index].get_adhesion_cell2_component(component_index);
};
//
double adhesion_system_class::get_coupling_constant(
						    const int & adhesion_index
						    ) const 
{
  return adhesions[adhesion_index].get_coupling_constant();
};
//
void adhesion_system_class::set_coupling_constant(
						    const int & adhesion_index,
						    const double & value
						    )  
{
  adhesions[adhesion_index].set_adhesion_double_value(
						      "coupling_constant",
						      value
						      );
};
//
std::string adhesion_system_class::get_interaction_type(
							 const int & adhesion_index
							 ) const 
{
  return adhesions[adhesion_index].get_interaction_type();
};
//
void adhesion_system_class::make_typepair_to_adhesion_maps()
{
  int pair_index;
  std::vector<int> counter; 
  typepair_to_adhesion_flags.clear();
  typepair_to_adhesion_pointers.clear();
  typepair_to_adhesion_pointers.push_back(0);
  typepair_to_adhesion_maps.clear();
  for (
       int type_index=0;
       type_index<number_of_cell_types;
       type_index++
       )
    {
      for(
	  int pairing_type_index=0;
	  pairing_type_index<number_of_cell_types;
	  pairing_type_index++
	  )
	{
	  typepair_to_adhesion_flags.push_back(false);
	  typepair_to_adhesion_pointers.push_back(0);
	  counter.push_back(0);
	};
    };
  for (
       int adhesion_index=0;
       adhesion_index< number_of_adhesion;
       adhesion_index++
       )
    {
      pair_index=number_of_cell_types
	*adhesions[adhesion_index].get_adhesion_cell_type(1)
	+adhesions[adhesion_index].get_adhesion_cell_type(2);
      typepair_to_adhesion_flags[pair_index]=true;
      typepair_to_adhesion_pointers[pair_index+1]++;
      if(adhesions[adhesion_index].get_adhesion_cell_type(2)
	 !=adhesions[adhesion_index].get_adhesion_cell_type(1))
	{
	  pair_index=number_of_cell_types
	    *adhesions[adhesion_index].get_adhesion_cell_type(2)
	    +adhesions[adhesion_index].get_adhesion_cell_type(1);
	  typepair_to_adhesion_pointers[pair_index+1]++;
	  typepair_to_adhesion_flags[pair_index]=true;
	};
    };
  //
  for(
      int type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      for(
	  int pairing_type_index=0;
	  pairing_type_index<number_of_cell_types;
	  pairing_type_index++
	  )
	{
	  pair_index=number_of_cell_types*type_index+pairing_type_index;
	  typepair_to_adhesion_pointers[pair_index+1]
	    += typepair_to_adhesion_pointers[pair_index];
	};
    };
  //
  typepair_to_adhesion_maps.reserve(
				    typepair_to_adhesion_pointers[
								  number_of_cell_types
								  *number_of_cell_types
								  ]
				    );
  for(
      int tmp_pair_index=0;
      tmp_pair_index<typepair_to_adhesion_pointers[number_of_cell_types*number_of_cell_types];
      tmp_pair_index++
      )
    {
      typepair_to_adhesion_maps.push_back(-1); 
    };
  //
  for(
      int adhesion_index=0;
      adhesion_index< number_of_adhesion;
      adhesion_index++
      )
    {
      pair_index=number_of_cell_types*adhesions[adhesion_index].get_adhesion_cell_type(1)
	+adhesions[adhesion_index].get_adhesion_cell_type(2);
      typepair_to_adhesion_maps[
				typepair_to_adhesion_pointers[
							      pair_index
							      ]
				+counter[pair_index]
				]=adhesion_index;
      counter[pair_index]++;
      if(adhesions[adhesion_index].get_adhesion_cell_type(1)
	 !=adhesions[adhesion_index].get_adhesion_cell_type(2))
	{
	  pair_index=number_of_cell_types*adhesions[adhesion_index].get_adhesion_cell_type(2)
	    +adhesions[adhesion_index].get_adhesion_cell_type(1);
	  typepair_to_adhesion_maps[
				    typepair_to_adhesion_pointers[
								  pair_index
								  ]
				    +counter[pair_index]
				    ]=adhesion_index;
	  counter[pair_index]++;
	};
    };
  //
  for(
      int map_index=0;
      map_index<typepair_to_adhesion_pointers[
					      number_of_cell_types*number_of_cell_types
					      ];
      map_index++
      )
    {
      coupling_constants.push_back(0.0);
      for(
	  int component_index=0;
	  component_index<number_of_adhesion_components;
	  component_index++
	  )
	{
	  polar_components_1.push_back(0.0);
	  polar_components_2.push_back(0.0);
	}
    }
  //
  std::vector<double> work_vector;
  for(
      int component_index=0;
      component_index<number_of_adhesion_components;
      component_index++
      )
    {
      work_vector.push_back(0.0);
    };
  //
  for (
       int type_index=0;
       type_index<number_of_cell_types;
       type_index++
       )
    {
      for(
	  int pairing_type_index=0;
	  pairing_type_index<number_of_cell_types;
	  pairing_type_index++
	  )
	{
	  pair_index=number_of_cell_types*type_index+pairing_type_index;
	  for(
	      int map_index=typepair_to_adhesion_pointers[
							  pair_index
							  ];
	      map_index< typepair_to_adhesion_pointers[
						       pair_index+1
						       ];
	      map_index++
	      )
	    {
	      coupling_constants[
				 map_index
				 ]
		=adhesions[typepair_to_adhesion_maps[map_index]].get_coupling_constant();
	      adhesions[typepair_to_adhesion_maps[map_index]].get_adhesion_polarity(
										    1,
										    work_vector
										    );
	      //
	      if(
		 adhesions[typepair_to_adhesion_maps[map_index]].get_adhesion_cell_type(1)
		 ==type_index
		 )
		{
		  for(
		      int component_index=0;
		      component_index<number_of_adhesion_components;
		      component_index++
		      )
		    {
		      polar_components_1[
					 map_index*number_of_adhesion_components+component_index
					 ]
			=work_vector[component_index];
		    };
		}
	      else if(
		      adhesions[typepair_to_adhesion_maps[map_index]].get_adhesion_cell_type(1)
		      ==pairing_type_index
		      )
		{
		  for(
		      int component_index=0;
		      component_index<number_of_adhesion_components;
		      component_index++
		      )
		    {
		      polar_components_2[
					 map_index*number_of_adhesion_components+component_index
					 ]
			=work_vector[component_index];
		    };
		}
	      else
		{
		  io_method.error_output(
					 "cellular_potts_adhesion",
					 "make_typepair_to_adhesion_maps",
					 "due to irrigal cell_type 1 in adhesion input."
					 );
		};
	      //
	      adhesions[typepair_to_adhesion_maps[map_index]].get_adhesion_polarity(
										    2,
										    work_vector
										    );
	      //
	      if(
		 adhesions[typepair_to_adhesion_maps[map_index]].get_adhesion_cell_type(2)
		 ==type_index
		 &&
		 type_index
		 !=
		 pairing_type_index
		 )
		{
		  for(
		      int component_index=0;
		      component_index<number_of_adhesion_components;
		      component_index++
		      )
		    {
		      polar_components_1[
					 map_index*number_of_adhesion_components+component_index
					 ]
			=work_vector[component_index];
		    };
		}
	      else if(
		      adhesions[typepair_to_adhesion_maps[map_index]].get_adhesion_cell_type(2)
		      ==pairing_type_index
		      )
		{
		  for(
		      int component_index=0;
		      component_index<number_of_adhesion_components;
		      component_index++
		      )
		    {
		      polar_components_2[
					 map_index*number_of_adhesion_components+component_index
					 ]
			=work_vector[component_index];
		    };
		}
	      else
		{
		  io_method.standard_output(
					    io_method.int_to_string(type_index)+','
					    +io_method.int_to_string(pairing_type_index)+','
					    +io_method.int_to_string(adhesions[typepair_to_adhesion_maps[map_index]].get_adhesion_cell_type(1))+','
					    +io_method.int_to_string(adhesions[typepair_to_adhesion_maps[map_index]].get_adhesion_cell_type(2))
					    );
		  io_method.error_output(
					 "cellular_potts_adhesion",
					 "make_typepair_to_adhesion_maps",
					 "due to irrigal cell_type 2 in adhesion input."
					 );
		};
	    };
	};
    };
};
//
void adhesion_system_class::interaction_key_check(const std::string & value)
{
  if(
     value !=adhesion_type_normal
     &&
     value !=adhesion_type_tight
     )
    {
      io_method.error_output(
			     "cellular_potts_adhesion",
			     "interaction_key_check",
			     "due to irrigal interaction type"+ value + "in adhesion input."
			     );
    };
};
//
adhesion_system_class::adhesion_system_class(
					     const model_parameters_cellular_potts_class & model,
					     const type_system_class & type_system
					     )
{
  space_dimension=model.get_space_dimension();
  number_of_cells=model.get_number_of_cells();
  number_of_cell_types=model.get_number_of_cell_types();
  number_of_adhesion=model.get_number_of_adhesion();
  number_of_adhesion_components=model.get_number_of_adhesion_components();
  adhesion_type_normal="normal";
  adhesion_type_tight="tight";
};
//
std::string adhesion_system_class::set_structure_adhesion_item(
							       const int & adhesion_index,
							       const std::string & child
							       )
{
  return io_method.generate_structure(
				      "adhesion",
				      adhesion_index,
				      child
				      );
};
//
