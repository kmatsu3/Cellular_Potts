#include "cellular_potts_cell.hpp"
  /*======================
    Constructor
   =======================*/
cell_cellular_potts_class::cell_cellular_potts_class()
{
  Type=-1;
};
//
void cell_cellular_potts_class::set_type(const int & input_type)
{
  Type=input_type;
};
//
cell_system_class::cell_system_class(
				     const model_parameters_cellular_potts_class & model,
				     const type_system_class & type_system
				     )
{
  int_error_number=-1;
  long int cell_index;
  int additional=0;
  space_dimension=model.get_space_dimension();
  number_of_sites=model.get_number_of_sites();
  number_of_cell_types=model.get_number_of_cell_types();
  number_of_cells=model.get_number_of_cells();
  if(type_system.get_buffer_type()!=-1)
    {
      additional=1;
      set_buffer_cell(
		      model,
		      type_system
		      );
    }
  cell_cellular_potts_class work_cell;
  for(cell_index=0;cell_index<number_of_cells+additional;cell_index++)
    {
      work_cell.set_type(
			 type_define(
				     model,
				     type_system,
				     cell_index
				     )
			 );
      //
      cells.push_back(work_cell);
    };
  //	
  int direction_index;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      work_positions.push_back(0);
    };
  //
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      system_dimensions.push_back(
				  model.get_system_dimension(direction_index)	  
				  );
    };
  //
  for(long int cell_index=0;cell_index<number_of_cells+1;cell_index++)
    {
      mobile_table.push_back(true);
    };
  set_mobile_table();
};
//
int cell_system_class::type_define(
				   const model_parameters_cellular_potts_class & model,
				   const type_system_class & type_system,
				   const long int & cell_index
				   ) const
{
  std::vector<long int> bound_list;
  int type_index;
  int type=-1;
  bound_list.push_back(0); 
  for (type_index=0;type_index<number_of_cell_types;type_index++)
    {
      if(type_system.get_number_of_cells(type_index)>0)
	{
	  bound_list.push_back(bound_list[type_index]+type_system.get_number_of_cells(type_index));
	} 
      else 
	{
	  bound_list.push_back(bound_list[type_index]);
	}
    };
  for (type_index=0;type_index<number_of_cell_types;type_index++)
    {
      if(bound_list[type_index]<=cell_index&&cell_index<bound_list[type_index+1]) type=type_index;
    };
  if(cell_index==buffer_cell) type=type_system.get_buffer_type();
  if(type==-1)
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "cell_system_class",
			     "type_define",
			     "due to inconsistency between type define and the number of cells"
			     );
    };
  return type;
};
//
void cell_system_class::show_cells(
				   const model_parameters_cellular_potts_class & model
				   ) const
{
  io_cellular_potts io_method;
  std::string message;
  long int cell_index;
  io_method.standard_output("=== Cell Data ===");
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      message ="Cell[";
      message+=io_method.longint_to_string(cell_index);
      message+="] = ";
      message+=io_method.int_to_string(cells[cell_index].get_type());
      io_method.standard_output(message);
    };
};
int cell_cellular_potts_class::get_type()
const {
  return Type;
};
//
long int cell_system_class:: get_buffer_cell()
  const{
  return buffer_cell;
};
//
void cell_system_class:: set_buffer_cell(
					 const model_parameters_cellular_potts_class & model,
					 const type_system_class & type_system
					 )
{
  if(type_system.get_buffer_type()!=-1)
    {
      buffer_cell=model.get_number_of_cells();
    };
};
//
int cell_system_class::get_type(const long int & cell_index)
  const { 
  if((long int)cells.size()>cell_index&&cell_index>=0)
    {
      return cells[cell_index].get_type();
    }else{
    return int_error_number;
  };
};
//
void cell_system_class::set_type(
				 const long int & cell_index,
				 const int & type_index
				 )
{ 
  cells[cell_index].set_type(type_index);
};
//
void cell_system_class::calculate_cell_volumes(
					       const model_parameters_cellular_potts_class & model,
					       const std::vector<long int> & configuration,
					       std::vector<long long int> & cell_volumes
					       ) 
  const{
  long long int site_index;
  for(int cell_index=0;cell_index<number_of_cells;cell_index++)
  {
    cell_volumes[cell_index]=0;
  };
  for(site_index=0;site_index<number_of_sites;site_index++)
    {
      if(configuration[site_index]>=0&&configuration[site_index]<number_of_cells)
       cell_volumes[configuration[site_index]]++;
    };
  //
};
//
void cell_system_class::assign_cell_origins(
					    const model_parameters_cellular_potts_class & model,
					    const std::vector<long int> & configuration,
					    std::vector<long long int> & cell_origins
					    ) 
{
  for(long int cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      cell_origins[cell_index]=-1;
    };
  for(long long int site_index=0;site_index<number_of_sites;site_index++)
    {
      if(cell_origins[configuration[site_index]]==-1) 
	{
	  cell_origins[configuration[site_index]]=site_index;
	};
    };
  for(long int cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      if(cell_origins[cell_index]==-1)
	{
	  std::string error_message="due to cell" + io_method.longint_to_string(cell_index) + "disappear !";
	  io_method.error_output(
				 "cell_system_class",
				 "type_define",
				 error_message
				 );
	};
    };
};
//
void cell_system_class::calculate_cell_total_positions(
						       const model_parameters_cellular_potts_class & model,
						       site_system_class & site_system,
						       const std::vector<long int> & configuration,
						       const std::vector<long long int> & cell_origins,
						       std::vector<long long int> & cell_positions
						       ) 
const {
  long long int site_index;
  long int cell_index;
  int direction_index;
  std::vector<long int> work_position(space_dimension);
  //
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      for(cell_index=0;cell_index<number_of_cells;cell_index++)
	{
	  cell_positions[cell_index*space_dimension+direction_index]=0;
	};
    };
  //
    /*
    // debug begin
    std::vector<long int> debug_coordinates(space_dimension,0);
    int counter=0;
    io_cellular_potts io_method;
    std::string message;
    // debug end
    */
    /*
    //debug begin
    message="debug totalcord[origin]: ";
    site_system.get_site_coordinates(cell_origins[configuration[0]],debug_coordinates);
    for(int direction_index=0;direction_index<space_dimension;direction_index++)
    {
      message+=io_method.longint_to_string(cell_positions[direction_index])+",";
    };
    io_method.standard_output(message);
    // debug end
    */
    //
  for(site_index=0;site_index<number_of_sites;site_index++)
    {
      if(configuration[site_index]!=buffer_cell)
	{
    /*
    //debug begin
    message="debug totalcord["+ 
    io_method.int_to_string(counter)+ "]: ";
    counter++;
    site_system.get_site_coordinates(site_index,debug_coordinates);
    for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      message+=io_method.longint_to_string(debug_coordinates[direction_index])+",";
    };
    io_method.standard_output(message);
    //ddebug end
    */
	  site_system.origin_shift(
				   model,
				   cell_origins[configuration[site_index]],
				   site_index,
				   work_position
				   );
	  //
	  for(direction_index=0;direction_index<space_dimension;direction_index++)
	    {
	      cell_positions[space_dimension*configuration[site_index]+direction_index]
	       +=work_position[direction_index];
	    };
	  //
	};
      //
    };
  //
    /*
    //debug
    message="debug totalcord[total]: ";
    for(int direction_index=0;direction_index<space_dimension;direction_index++)
    {
      message+=io_method.longint_to_string(cell_positions[direction_index])+",";
    };
    io_method.standard_output(message);
    //
    */
};
//
/*
const std::vector<long long int> cell_system_class::calculate_cell_total_positions(
										   const model_parameters_cellular_potts_class & model,
										   site_system_class & site_system,
										   const std::vector<long int> & configuration,
										   const std::vector<long long int> & cell_origins
										   ) 
const {
  int space_dimension=model.get_space_dimension();
  long long int number_of_sites=model.get_number_of_sites();
  long long int site_index;
  long int number_of_cells=model.get_number_of_cells();
  long int cell_index;
  int direction_index;
  std::vector<long long int> cell_positions(space_dimension*number_of_cells);
  std::vector<long int> positions(number_of_sites*space_dimension,0);
  std::vector<long int> work_position(space_dimension);
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      cell_positions.push_back(0);
      cell_positions.push_back(0);
    };
  //
  for(site_index=0;site_index<number_of_sites;site_index++)
    {
      //   fprintf(stderr,"%ld,%ld,%lld\n",site_index,configuration[site_index],cell_origins[configuration[site_index]]);
      if(configuration[site_index]!=buffer_cell)
	{
	  //
	  //
	  //work_position=site_system.origin_shift(
	//					 model,
	//					 cell_origins[configuration[site_index]],
	//					 site_index
	//					 );
	  //
	  site_system.origin_shift(
				   model,
				   cell_origins[configuration[site_index]],
				   site_index,
				   work_position
				   );
	  //
	  for(direction_index=0;direction_index<space_dimension;direction_index++)
	    {
	      positions[space_dimension*site_index+direction_index]=work_position[direction_index];
	    };
	  //
	};
      //
    };
  //
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      for(site_index=0;site_index<number_of_sites;site_index++)
	{
	  if(configuration[site_index]!=buffer_cell)
	    {
	      cell_positions[space_dimension*configuration[site_index]+direction_index]
		+=positions[space_dimension*site_index+direction_index];
	    };
	};
    };
  //
  return cell_positions;
};
*/
//		    
std::vector<double> cell_system_class::calculate_cell_position(
							       const model_parameters_cellular_potts_class & model,
							       const site_system_class & site_system,
							       const std::vector<long long int> & cell_origins,
							       const std::vector<long long int> & cell_total_positions,
							       const std::vector<long long int> & cell_volumes,
							       const long int & cell_index
							       )
const {
  int space_dimension=model.get_space_dimension();
  std::vector<double> return_value(space_dimension);
  int component_index;
  std::vector<long int> origin_position(space_dimension);
  site_system.get_site_coordinates(
				   cell_origins[cell_index],
				   origin_position
				   );
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      return_value[component_index]
	=(double)cell_total_positions[space_dimension*cell_index+component_index]
	/(double)cell_volumes[cell_index]
	+(double)origin_position[component_index];
    };
  return return_value;
};
//
//		    
void cell_system_class::calculate_cell_position_from_cell_origin(
								 const model_parameters_cellular_potts_class & model,
								 const site_system_class & site_system,
								 const std::vector<long long int> & cell_origins,
								 const std::vector<long long int> & cell_total_positions,
								 const std::vector<long long int> & cell_volumes,
								 const long int & cell_index,
								 std::vector<double> & relative_cell_position 
								 )
const {
  int space_dimension=model.get_space_dimension();
  int component_index;
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      relative_cell_position[component_index]
	=(double)cell_total_positions[space_dimension*cell_index+component_index]
	/(double)cell_volumes[cell_index];
    };
};
//		    
void cell_system_class::calculate_cell_position_in_system(
							  const model_parameters_cellular_potts_class & model,
							  const site_system_class & site_system,
							  const std::vector<long long int> & cell_origins,
							  const std::vector<long long int> & cell_total_positions,
							  const std::vector<long long int> & cell_volumes,
							  const long int & cell_index,
							  std::vector<double> & relative_cell_position 
							  )
  const {
  int space_dimension=model.get_space_dimension();
  int component_index;
  std::vector<long int> origin_position(space_dimension);
  site_system.get_site_coordinates(
				   cell_origins[cell_index],
				   origin_position
				   );
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      relative_cell_position[component_index]
	=
	(double)cell_total_positions[space_dimension*cell_index+component_index]
	/(double)cell_volumes[cell_index]
	+(double)origin_position[component_index];
      if(relative_cell_position[component_index]<0.0)
	{
	  relative_cell_position[component_index]+=(double)system_dimensions[component_index];
	}
      else if(relative_cell_position[component_index]>system_dimensions[component_index])
	{
	  relative_cell_position[component_index]-=(double)system_dimensions[component_index];
	};
    };
};
//
void cell_system_class::show_cell_positions(
					    const model_parameters_cellular_potts_class & model,
					    const site_system_class & site_system,
					    const std::vector<long int> & configuration,
					    const std::vector<long long int> & cell_origins,
					    const std::vector<long long int> & cell_total_positions,
					    const std::vector<long long int> & cell_volumes
					    )
const   {
  long long int cell_index;
  long long int number_of_cells=model.get_number_of_cells();
  int space_dimension=model.get_space_dimension();
  std::vector<double> work_vector(space_dimension);
  io_cellular_potts io_method;
  std::string message;
  io_method.standard_output("==== Cell positions (CoM) ====");
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      work_vector=calculate_cell_position(
					  model,
					  site_system,
					  cell_origins,
					  cell_total_positions,
					  cell_volumes,
					  cell_index
					  );
      message ="Cell[";
      message+=io_method.longlongint_to_string(cell_index);
      message+="]=";
      message+=io_method.doublearray_to_string(work_vector);
      io_method.standard_output(message);
    };
};
//
void cell_system_class::check_cell_weight_polarity(
         const model_parameters_cellular_potts_class & model,
         site_system_class & site_system,
         const std::vector<long int> & configuration,
         const std::vector<long long int> & cell_origins,
         const std::vector<long long int> & cell_total_positions,
         const std::vector<long long int> & cell_volumes
         ) const
{
   io_cellular_potts io_method;
    /*
    // debug begin
    std::vector<long int> debug_coordinates(space_dimension,0);
    int counter=0;
    std::string message;
    // debug end
    */
  std::vector<long int> work_position(space_dimension);
  double work_norm;
  std::vector<double> work_polarity(space_dimension);
  std::vector<double> weight_polarity(space_dimension);
  for(int direction_index=0;direction_index<space_dimension;direction_index++)
  {
    weight_polarity[direction_index]=0.0;
  };
  for(long long int site_index=0;site_index<number_of_sites;site_index++)
    {
      if(configuration[site_index]!=buffer_cell)
      {
        site_system.origin_shift(
           model,
           cell_origins[configuration[site_index]],
           site_index,
           work_position
           );
        /*
        //debug begin
    message="debug totalcord["+ 
    io_method.int_to_string(counter)+ "]: ";
    counter++;
    site_system.get_site_coordinates(site_index,debug_coordinates);
    for(int direction_index=0;direction_index<space_dimension;direction_index++)
    {
      message+=io_method.longint_to_string(work_position[direction_index])+",";
    };
    message+="<-";
    for(int direction_index=0;direction_index<space_dimension;direction_index++)
    {
      message+=io_method.longint_to_string(debug_coordinates[direction_index])+",";
    };
    io_method.standard_output(message);
        //debug end
        */
          for(int direction_index=0;direction_index<space_dimension;direction_index++)
            {
              work_polarity[direction_index]
              =(double)work_position[direction_index]
              -(double)cell_total_positions[direction_index]
              /(double)cell_volumes[configuration[site_index]];
            };
            work_norm=tool.norm(work_polarity);
            if(work_norm!=0)
            {
              for(int direction_index=0;direction_index<space_dimension;direction_index++)
                {
                  weight_polarity[direction_index]
                  +=work_polarity[direction_index]
                  /work_norm;
                };
            };
            /*
        // debug begin
    message="debug direction[" + io_method.longlongint_to_string(counter)+"] " ;
    for(int direction_index=0;direction_index<space_dimension;direction_index++)
    {
      message+=io_method.double_to_string(
        (double)work_position[direction_index]-
        (double)cell_total_positions[direction_index]
        /(double)cell_volumes[configuration[site_index]]
          )+",";
    };
    message+=io_method.longint_to_string(cell_volumes[configuration[site_index]]);
    io_method.standard_output(message);
        // debug end
        */
      }; 
    };
    /*
    // debug begin
    message="debug center pos";
    for(int direction_index=0;direction_index<space_dimension;direction_index++)
    {
      message+=io_method.double_to_string((double)cell_total_positions[direction_index]
        /(double)cell_volumes[0])+",";
    };
    io_method.standard_output(message);
//    std::string message="debug: ";
*/
    std::string message="debug: ";
    for(int direction_index=0;direction_index<space_dimension;direction_index++)
    {
      message+=io_method.double_to_string(weight_polarity[direction_index])+",";
    };
  io_method.standard_output(message);
  // debug end
};
//
void cell_system_class::set_mobile_table()
{
  std::vector<long int> fixed_cells;
  fixed_cells.clear();
  std::string structure_item;
  io_method.get_input_longint_array(
				    "fixed_table",
				    "fixed_table.cells",
				    fixed_cells
				    );
  std::vector<long int>::iterator fixed_cell_index=fixed_cells.begin();
  io_method.standard_output("debug="+io_method.longint_to_string(fixed_cells.size()));
  while(fixed_cell_index!=fixed_cells.end())
    {
      mobile_table[(*fixed_cell_index)]=false;
      fixed_cell_index++;
    };
  io_method.standard_output("=== mobile cell list ===");
    for(
	long int cell_index=0;
	cell_index<number_of_cells;
	cell_index++
	)
      {
	if(mobile_table[cell_index])
	  {
	    io_method.standard_output(
				      "Cell ["
				      + io_method.longint_to_string(cell_index)
				      + "] : Mobile"
				      );
	  }
	else
	  {
	    io_method.standard_output(
				      "Cell ["
				      + io_method.longint_to_string(cell_index)
				      + "] : Unmobile"
				      );
	  };
      };
};
