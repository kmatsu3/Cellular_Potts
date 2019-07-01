#include "cellular_potts_motion_trace.hpp"
  /*======================
    Methods
   =======================*/
/*
void cell_displacement_class::update(
				     state_system_class & state
				     )
{
  //
  if(cell_tracking_flag=="on")
    {
      state.get_cell_displacements(
				   cell_displacements
				   );
      //
      std::vector<double>::iterator component_index=cell_displacements.begin();
      int counter = 0;
      // for debug should be muted!!!!
	if(cell_displacements.size()!=cell_total_displacements.size())
	{
	io_method.error_output(
	"cell_displacement_class",
			     "update",
			     "inconsistency between memories for displacements."
			     );
			     };
      // for debug should be muted!!!!
      while(component_index!=cell_displacements.end())
	{
//      std::string message2=io_method.int_to_string(counter);
//      io_method.standard_output("ok?"+message2);
	  cell_total_displacements[counter]
	= cell_total_displacements[counter]
	    + (*component_index);
	  component_index++;
	  counter++;
	  //debug
	  //      std::string message=io_method.double_to_string(cell_total_displacements[counter-1]);
	  //      io_method.standard_output(message);
	  //
	};
    };
};
*/
//
void cell_displacement_class::push(
				   std::vector<double> & input_cell_total_displacements
				   )
{
  for(
      long int cell_index=type_to_cell_pointer[0];
      cell_index<type_to_cell_pointer[1];
      cell_index++
      )
    {
      for(
	  int direction_index=0;
	  direction_index < space_dimension;
	  direction_index++
	  )
	{
	  cell_total_displacements[(cell_index-type_to_cell_pointer[0])*space_dimension+direction_index]
	    =input_cell_total_displacements[
					    cell_index*space_dimension
					    +direction_index
					    ];
	};
    };
};
//
void cell_displacement_class::get_average_mean_square_displacement(
								   double & cell_averaged_displacement
								   )
const{
  if(cell_tracking_flag=="on")
    {
      cell_averaged_displacement
	= std::inner_product(
			     cell_total_displacements.begin(),
			     cell_total_displacements.end(),
			     cell_total_displacements.begin(),
			     0.0
			     );
      /*
      double work_double=0.0;
      for(int cell_index=0;cell_index<number_of_cells;cell_index++)
	{
	  for(int direction_index=0;direction_index<space_dimension;direction_index++)
	    {
	      work_double
		+= cell_total_displacements[space_dimension*cell_index+direction_index]
		*  cell_total_displacements[space_dimension*cell_index+direction_index];
	    };
	  cell_averaged_displacement
	    += work_double;
	};
      */
      cell_averaged_displacement
	= cell_averaged_displacement/((double)number_of_cells);
    };
};
//
void cell_displacement_class::push_average_mean_square_displacement()
{
  if(cell_tracking_flag=="on")
    {
      double work_double=0.0;
      get_average_mean_square_displacement(work_double);
      average_mean_square_displacement.push_back(work_double);
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
        {
          first_cell_trajectory.push_back(
					  cell_total_displacements[
								   direction_index
								   ]
					  );
        };
    };
};
//
void cell_displacement_class::output_average_mean_square_displacement() const
{
  if(cell_tracking_flag=="on")
    {
      int counter=0;
      std::string message ="#";
      message+="/cell-averaged<";
      message+=io_method.int_to_string(average_mean_square_displacement.size());
      io_method.output_message(message,"cell_track.txt");
      std::vector<double>::const_iterator index=average_mean_square_displacement.begin();
      while(index!=average_mean_square_displacement.end())
	    { 
	     message ="";
	     message+=io_method.double_to_string(*index);
       for(int direction_index=0;direction_index<space_dimension;direction_index++)
          {
            message+=" ";
            message+=io_method.double_to_string(first_cell_trajectory[space_dimension*counter+direction_index]);
          }
	     io_method.output_message(message,"cell_track.txt");  
	     index++;
       counter++;
	   };
    };
};
//
void cell_displacement_class::output_displacement() const
{
  if(cell_tracking_flag=="on")
    {
      std::string message ="#";
      message+="/cell-averaged<";
      message+=io_method.int_to_string(cell_total_displacements.size());
      io_method.output_message(message,"cell_displacement.txt");
      std::vector<double>::const_iterator index=cell_total_displacements.begin();
      int counter=0;
      message = "";
      while(index!=average_mean_square_displacement.end())
      { 
        if(counter%space_dimension==space_dimension-1)
          {
            message+=io_method.double_to_string(*index);
            io_method.output_message(message,"cell_displacement.txt");  
            index++;
            message = "";
            counter = 0;
          }
          else
          {
            message += " ";
            message += io_method.double_to_string(*index);
            counter++;
            index++;
          }
      };
    };
};
//
long int cell_displacement_class::get_size_of_memory()
{
  return (long int)average_mean_square_displacement.size();
};
//
double cell_displacement_class::get_average_mean_square_displacement_component(
							      const long int & counter
							      )
{
  return average_mean_square_displacement[counter];
};
//
double cell_displacement_class:: get_first_cell_position(
							 const long int & counter,
							 const int & component_index
							 )
{
  return first_cell_trajectory[counter*space_dimension+component_index];
};
  /*======================
    Constructor
   =======================*/
cell_displacement_class::cell_displacement_class()
{
};
/*
cell_displacement_class::cell_displacement_class(
						 const model_parameters_cellular_potts_class & model,
						 const int & type_index,
						 const std::vector<long int> type_to_cell_pointer
						 )
{
  cell_tracking_flag=model.get_cell_tracking_flag();
  if(cell_tracking_flag=="on")
    {
      if(type_index==buffer_type)
	{
	  // Declarations for local values
	  int cell_index;
	  int direction_index;
	  // Loading of class members
	  space_dimension=model.get_space_dimension();
	  number_of_cells=model.get_number_of_cells();
	  //
	  for(cell_index=0;cell_index<number_of_cells;cell_index++)
	    {
	      for(direction_index=0;direction_index<space_dimension;direction_index++)
		{
		  cell_displacements.push_back(0.0);
		  cell_total_displacements.push_back(0.0);
		};
	    };
	};
      else if(
	      (type_index >=0)
	       &&
	      (type_index < type_to_cell_pointer.size())
	      )
	{
	  for(cell_index=type_to_cell_pointer(type_index);
	      cell_index<type_to_cell_pointer(type_index+1);
	      cell_index++)
	    {
	      for(direction_index=0;direction_index<space_dimension;direction_index++)
		{
		  cell_displacements.push_back(0.0);
		  cell_total_displacements.push_back(0.0);
		};
	    };
	}
      else
	{
	  io_method.error_output(
				 "cellular_potts_motion_trace",
				 "cell_displacement_class::cell_displacement_class"
				 "due to undefined cell type index is input!"
				 );
	};
    };
};
*/
  /*======================
    Initializer
   =======================*/
 //
void cell_displacement_class::initialize(
					 const model_parameters_cellular_potts_class & model,
					 const type_system_class & type_system,
					 const int & input_type_index,
					 const std::vector<long int> input_type_to_cell_pointer
					 )
{
  // Declarations for local values
  int cell_index;
  int direction_index;
  type_index=input_type_index;
  space_dimension=model.get_space_dimension();
  if(type_index==type_system.get_buffer_type())
    {
      number_of_cells=model.get_number_of_cells();
    }
  else
    {
      number_of_cells=type_system.get_number_of_cells(type_index);
      if(
	 (number_of_cells<1)
	 ||
	 (number_of_cells>model.get_number_of_cells())
	 )
	{
	  io_method.error_output(
				 "cellular_potts_motion_trace",
				 "cell_displacement_class::initialize",
				 "due to irrigal number of cells for type:"
				 + io_method.int_to_string(type_index)
				 + ", of value = "
				 + io_method.longint_to_string(number_of_cells)
				 + " !"
				 );
	};
    };
  type_to_cell_pointer.clear();
  // Loading of class members
  cell_tracking_flag=model.get_cell_tracking_flag();
  if(cell_tracking_flag=="on")
    {
      if(type_index==type_system.get_buffer_type())
	{
	  //
	  type_to_cell_pointer.push_back(0);
	  type_to_cell_pointer.push_back(number_of_cells);
	  for(cell_index=0;cell_index<number_of_cells;cell_index++)
	    {
	      for(direction_index=0;direction_index<space_dimension;direction_index++)
		{
		  cell_displacements.push_back(0.0);
		  cell_total_displacements.push_back(0.0);
		};
	    };
	}
      else if(
	      (type_index >=0)
	       &&
	      (type_index < (int)input_type_to_cell_pointer.size())
	      )
	{
	  type_to_cell_pointer.push_back(input_type_to_cell_pointer[type_index]);
	  type_to_cell_pointer.push_back(input_type_to_cell_pointer[type_index+1]);	  
	  for(
	      cell_index=type_to_cell_pointer[0];
	      cell_index<type_to_cell_pointer[1];
	      cell_index++
	      )
	    {
	      for(int direction_index=0;direction_index<space_dimension;direction_index++)
		{
		  cell_displacements.push_back(0.0);
		  cell_total_displacements.push_back(0.0);
		};
	    };
	}
      else
	{
	  io_method.error_output(
				 "cellular_potts_motion_trace",
				 "cell_displacement_class::cell_displacement_class",
				 "due to undefined cell type index is input"
				 + io_method.int_to_string(type_index)
				 + " !"
				 );
	};
    };
};
 //
cell_displacement_system::cell_displacement_system(
						   const model_parameters_cellular_potts_class & model,
						   const type_system_class & type_system
						   )
{
  buffer_type=type_system.get_buffer_type();
  number_of_cell_types=model.get_number_of_cell_types();
  number_of_cells=model.get_number_of_cells();
  space_dimension=model.get_space_dimension();
  cell_tracking_flag=model.get_cell_tracking_flag();
  for(
      long int cell_index=0;
      cell_index<number_of_cells;
      cell_index++
      )
    {
      for(int direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  cell_displacements.push_back(0.0);
	  cell_total_displacements.push_back(0.0);
	};
    };
  type_to_cell_pointer.clear();
  type_to_cell_pointer.push_back(0);
  std::string messages;
  for(
      int type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      int incriment = type_system.get_number_of_cells(type_index);
      if(incriment < 0 ) incriment =0;
      type_to_cell_pointer.push_back(0);
      type_to_cell_pointer[type_index+1]
	= type_to_cell_pointer[type_index]
	+ incriment;
    };
  //
  cell_displacement_for_types.clear();
  cell_displacement_class work_displacement;
  for(
      int type_index=0;
      type_index<number_of_cell_types;
      type_index++
      )
    {
      cell_displacement_for_types.push_back(work_displacement);
      cell_displacement_for_types[type_index].initialize(
							 model,
							 type_system,
							 type_index,
							 type_to_cell_pointer
							 );
    }
};
//
void cell_displacement_system::update(
				      state_system_class & state
				      )
{
  //
  if(cell_tracking_flag=="on")
    {
      state.get_cell_displacements(
				   cell_displacements
				   );
      //
      std::vector<double>::iterator component_index=cell_displacements.begin();
      int counter = 0;
      /*
	if(cell_displacements.size()!=cell_total_displacements.size())
	{
	io_method.error_output(
	"cell_displacement_class",
	"update",
	"inconsistency between memories for displacements."
	);
	};
      */
      while(component_index!=cell_displacements.end())
	{
	  //      std::string message2=io_method.int_to_string(counter);
	  //      io_method.standard_output("ok?"+message2);
	  cell_total_displacements[counter]
	    = cell_total_displacements[counter]
	    + (*component_index);
	  component_index++;
	  counter++;
	  //debug
	  //      std::string message=io_method.double_to_string(cell_total_displacements[counter-1]);
	  //      io_method.standard_output(message);
	  //
	};
    };
};
void cell_displacement_system::push()
{
  //
  if(cell_tracking_flag=="on")
    {
      for(
	  int type_index = 0;
	  type_index < number_of_cell_types;
	  type_index ++
	  )
	{
	  cell_displacement_for_types[type_index].push(cell_total_displacements);
	  cell_displacement_for_types[type_index].push_average_mean_square_displacement();
	};
    };
};
void cell_displacement_system::output()
{
  if(cell_tracking_flag=="on")
    {
      io_method.file_initialize("cell_track.txt");
      //     int counter=0;
      std::string message ="#";
      message+="/cell-averaged<";
      //      message+=io_method.int_to_string(average_mean_square_displacement.size());
      io_method.output_message(message,"cell_track.txt");
      /*
      while(index!=average_mean_square_displacement.end())
	    { 
	     message ="";
	     message+=io_method.double_to_string(*index);
       for(int direction_index=0;direction_index<space_dimension;direction_index++)
          {
            message+=" ";
            message+=io_method.double_to_string(first_cell_trajectory[space_dimension*counter+direction_index]);
          }
	     io_method.output_message(message,"cell_track.txt");  
	     index++;
       counter++;
       };
      */
      for(
	  long int memory_index=0;
	  memory_index<cell_displacement_for_types[0].get_size_of_memory();
	  memory_index++
	  )
	{
	  message = "";
	  for(
	      int type_index=0;
	      type_index<number_of_cell_types;
	      type_index++
	      )
	    {
	      message += " " + io_method.int_to_string(type_index);
	      message += " " + io_method.double_to_string(
							  cell_displacement_for_types[type_index].get_average_mean_square_displacement_component(
																		 memory_index
																		 )
							  );
	      for(
		  int component_index=0;
		  component_index<space_dimension;
		  component_index++
		  )
		{ 
		  message += " " + io_method.double_to_string(
							      cell_displacement_for_types[type_index].get_first_cell_position(
																 memory_index,
																 component_index
										      )
							      );
		};
	    };
	  io_method.output_message(message,"cell_track.txt");
	};
    };
};
