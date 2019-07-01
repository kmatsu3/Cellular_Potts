#include "cellular_potts_state.hpp"

void polarity_motion_class::initialize_polarities(
						  const model_parameters_cellular_potts_class & model,
						  const cell_system_class & cell_system,
						  const site_system_class & site_system,
						  const std::vector<double> & input_polarities,
						  const std::vector<long long int> & input_cell_total_positions,
						  const std::vector<long long int> & input_cell_volumes,
						  const std::vector<long long int> & input_cell_origins
						  )
{
  //
  long int cell_index;
  int direction_index;
  //
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      //
      work_vector_double_1=cell_system.calculate_cell_position(
							       model,
							       site_system,
							       input_cell_origins,
							       input_cell_total_positions,
							       input_cell_volumes,
							       cell_index
							       );
      //
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  memorized_positions[cell_index*space_dimension+direction_index]
	    =work_vector_double_1[direction_index];
	};
      //
    };
};
//
void polarity_motion_class::get_displacements(
					      const model_parameters_cellular_potts_class & model,
					      const cell_system_class & cell_system,
					      const site_system_class & site_system,
					      const std::vector<long long int> & input_cell_total_positions,
					      const std::vector<long long int> & input_cell_volumes,
					      const std::vector<long long int> & input_cell_origins,
					      std::vector<double> & cell_displacements
					      )
  {
  long int cell_index;
  int component_index;
  /* Calculate Cell positions with a origin set at a site in cell*/
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      //
      work_vector_double_1=cell_system.calculate_cell_position(
							       model,
							       site_system,
							       input_cell_origins,
							       input_cell_total_positions,
							       input_cell_volumes,
							       cell_index
							       );
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  cell_displacements[cell_index*space_dimension+component_index]
	    =work_vector_double_1[component_index];
	};
    };
  //
   for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      for(component_index=0;component_index<space_dimension;component_index++)
  	{
  	  cell_displacements[cell_index*space_dimension+component_index]
  	    =cell_displacements[cell_index*space_dimension+component_index]
	    -memorized_positions[cell_index*space_dimension+component_index];
  	};
    };
};
//
void polarity_motion_class::update_polarities(
					      const model_parameters_cellular_potts_class & model,
					      const cell_system_class & cell_system,
					      const site_system_class & site_system,
					      const std::vector<long long int> & input_cell_total_positions,
					      const std::vector<long long int> & input_cell_volumes,
					      const std::vector<long long int> & input_cell_origins,
					      const std::vector<double> & adhesion_field,
					      std::vector<double> & polarities
					      ) 
  {
    //
    std::string message="";
    //
  std::vector<double> work_vector_local_double(space_dimension);
  /* Calculate Cell positions with a origin set at a site in cell*/
  for(long int cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      //
      work_vector_local_double=cell_system.calculate_cell_position(
								   model,
								   site_system,
								   input_cell_origins,
								   input_cell_total_positions,
								   input_cell_volumes,
								   cell_index
								   );
      for(int component_index=0;component_index<space_dimension;component_index++)
	{
	  work_vector_double[cell_index*space_dimension+component_index]
	    =work_vector_local_double[component_index];
	};
    };
  //
   for(long int cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      // debug
      //if(cell_index==43) 
      //
      for(int component_index=0;component_index<space_dimension;component_index++)
  	{
  	  difference_polarities[cell_index*space_dimension+component_index]
  	    =(work_vector_double[cell_index*space_dimension+component_index]
  	      -memorized_positions[cell_index*space_dimension+component_index])
  	    /persistent_time[cell_types[cell_index]]
	    +adhesion_sensitivity[cell_types[cell_index]]
	    *adhesion_field[cell_index*space_dimension+component_index];
	  //	  if(cell_index==43)message+= io_method.double_to_string(difference_polarities[cell_index*space_dimension+component_index])+',';
	  //if(cell_index==43)message+= io_method.double_to_string(adhesion_field[cell_index*space_dimension+component_index])+',';
  	};
      // debug
      //if(cell_index==43)
      //{
      //  for(int direction_index=0;
      //      direction_index<space_dimension;
      //      direction_index++)
      //    {
      //      if(std::isnan(difference_polarities[cell_index*space_dimension+direction_index]))
      //	io_method.error_output("B0"+io_method.longint_to_string(direction_index),"",message);
      //	    }
      // debug
      //  for(int direction_index=0;
      //      direction_index<space_dimension;
      //      direction_index++)
      //    {
      //      if(std::isnan(adhesion_field[cell_index*space_dimension+direction_index]))
      //	/	    io_method.error_output("B1"+io_method.longint_to_string(direction_index),"",message);
      //    }
      //    for(int direction_index=0;
      //	  direction_index<space_dimension;
      //   direction_index++)
      //	if(std::isnan(difference_polarities[cell_index*space_dimension+direction_index]))
      //  {
      //    message = io_method.longint_to_string(cell_index)+',';
      //    message+= io_method.int_to_string(direction_index)+',';
      //    message+= io_method.double_to_string(difference_polarities[cell_index*space_dimension+direction_index]);
      //    message+= io_method.double_to_string(adhesion_field[cell_index*space_dimension+direction_index]);
      //    io_method.standard_output(message);
      //  };
      //};
      // debug
    };
  //
  for(long int cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      for(int component_index=0;component_index<space_dimension;component_index++)
	{
	  //
	  work_vector_double_1[component_index]=polarities[cell_index*space_dimension+component_index];
	  work_vector_double_2[component_index]=difference_polarities[cell_index*space_dimension+component_index];
	  //
	};
      //
      projection_to_parpendicular_direction
	=std::inner_product(
			    work_vector_double_1.begin(),
			    work_vector_double_1.end(),
			    work_vector_double_2.begin(),
			    0.0
			    )/
	std::inner_product(
			   work_vector_double_1.begin(),
			   work_vector_double_1.end(),
			   work_vector_double_1.begin(),
			   0.0
			   );
      //
      for(int component_index=0;component_index<space_dimension;component_index++)
	{
	  //
	  difference_polarities[cell_index*space_dimension+component_index]
	    =difference_polarities[cell_index*space_dimension+component_index]
	    -projection_to_parpendicular_direction
	    *polarities[cell_index*space_dimension+component_index];
	  //
	};
    };
  //
  //  fprintf(stderr,"bug is in hereafter");
  //
  if(space_dimension!=2)
    {
      for(long int cell_index=0;cell_index<number_of_cells;cell_index++)
	{
	  for(int component_index=0;component_index<space_dimension;component_index++)
	    {
	      //
	      work_vector_double_1[component_index]=polarities[cell_index*space_dimension+component_index]
		+difference_polarities[cell_index*space_dimension+component_index];
	      //
	    };
	  //	  fprintf(stderr,"(%f,%f,%f)::"
	  //	  ,polarities[cell_index*space_dimension+0]
	  //	  ,polarities[cell_index*space_dimension+1]
	  //	  ,polarities[cell_index*space_dimension+2]
	  //	  );
	  //	  fprintf(stderr,"(%f,%f,%f)::"
	  //	  ,work_vector_double_1[0]
	  //	  ,work_vector_double_1[1]
	  //	  ,work_vector_double_1[2]
	  //	  );
	  work_vector_double_2=tool.normalized_vector(work_vector_double_1);
	  // debug
	  //fprintf(stderr,"(%f,%f,%f)/(%f,%f)\n"
	  //		  ,work_vector_double_2[0]/work_vector_double_1[0]
	  //		  ,work_vector_double_2[1]/work_vector_double_1[1]
	  //		  ,work_vector_double_2[2]/work_vector_double_1[2]
	  //		  ,tool.norm(work_vector_double_1)
	  //		  ,tool.norm(work_vector_double_2)
	  //		  );
	  // debug
	  for(int component_index=0;component_index<space_dimension;component_index++)
	    {
	      //
	      polarities[cell_index*space_dimension+component_index]
		=work_vector_double_2[component_index];
	      //
	    };
	};
    }
  else
    {
      //  for(cell_index=0;cell_index<number_of_cells;cell_index++)
      //   {
      //   fprintf(
      //	      stderr,
      //	      "(%ld,%f,%f)\n",
      //	      cell_index,
      //	      polarities[cell_index*2+0]+difference_polarities[cell_index*2+0],
      //	      polarities[cell_index*2+1]+difference_polarities[cell_index*2+1]
      //	      );
      //    }
      for(long int cell_index=0;cell_index<number_of_cells;cell_index++)
	{
	  for(int component_index=0;component_index<space_dimension;component_index++)
	    {
	      //
	      work_vector_double_1[component_index]=polarities[cell_index*space_dimension+component_index];
	      work_vector_double_2[component_index]=difference_polarities[cell_index*space_dimension+component_index];
	      //
	    };
	  angle=tool.norm(work_vector_double_2)
	    /tool.norm(work_vector_double_1)
	    *external_product_sign(work_vector_double_1,work_vector_double_2);
	  work_vector_double_1[0]
	    =polarities[cell_index*space_dimension+0]*cos(angle)
	    -polarities[cell_index*space_dimension+1]*sin(angle);
	  work_vector_double_1[1]
	    =polarities[cell_index*space_dimension+0]*sin(angle)
	    +polarities[cell_index*space_dimension+1]*cos(angle);
	  polarities[cell_index*space_dimension+0]= work_vector_double_1[0];
	  polarities[cell_index*space_dimension+1]= work_vector_double_1[1];
	};
    };
  //
  //  for(long int cell_index=0;cell_index<number_of_cells;cell_index++)
  //  {
  //    fprintf(stderr,"(%ld,%f,%f)\n",cell_index,polarities[cell_index*2],polarities[cell_index*2+1]);
  //  }
};
//
/*
const std::vector<double> polarity_motion_class:: rotation_2d(
							     const std::vector<double> vectors,
							     const std::vector<double> difference_vectors,
							     const long int array_size
							     )
  const {
  std::vector<double> return_vector(vectors.size(),0.0);
  long int vector_index;
  std::vector<double> work;
  for (vector_index=0;vector_index<array_size;vector_index++)
    {
      rotation_angles[vector_index]
	=vectors[vector_index*2+1]*difference_vectors[vector_index*2+2]
	-vectors[vector_index*2+2]*difference_vectors[vector_index*2+1];
      rotation_angles[vector_index]
	=rotation_angles[vector_index]/fabs(rotation_angles[vector_index])
	*sqrt(difference_vectors[vector_index*2+1]*difference_vectors[vector_index*2+1]
	      +difference_vectors[vector_index*2+2]*difference_vectors[vector_index*2+2]);
      return_vector[vector_index*2+1]=cos(rotation_angles[vector_index])
	-sin(rotation_angles[vector_index]);
      return_vector[vector_index*2+2]=sin(rotation_angles[vector_index])
	+cos(rotation_angles[vector_index]);
    }
  return return_vector;
};
*/
//
inline double polarity_motion_class::external_product_sign(
							   const std::vector<double> & vector_1,
							   const std::vector<double> & vector_2
							   ) 
{
  //
  double external_product_12=0.0;
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
polarity_motion_class::polarity_motion_class(
					     const model_parameters_cellular_potts_class & model,
					     const cell_system_class & cell_system,
					     const type_system_class & cell_type_system
					     )
{
  long int cell_index;
  space_dimension=model.get_space_dimension();
  number_of_cells=model.get_number_of_cells();
  number_of_cell_types=model.get_number_of_cell_types();
  int direction_index;
  //
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      work_vector_double_1.push_back(0.0);
      work_vector_double_2.push_back(0.0);
    };
  //
  number_of_cells=model.get_number_of_cells();
  //
  cell_types.clear();
  for(cell_index=0;cell_index<number_of_cells;cell_index++)
    {
      //  rotaion_angles.push_back(0.0);
      cell_types.push_back(cell_system.get_type(cell_index));
      for(direction_index=0;direction_index<space_dimension;direction_index++)
	{
	  work_vector_double.push_back(0.0);
	  memorized_positions.push_back(0.0);
	  difference_polarities.push_back(0.0);
	  //	  rotation_axes.push_back(0.0);
	};
    };
  //
  number_of_cell_types=model.get_number_of_cell_types();
  //
  persistent_time.clear();
  adhesion_sensitivity.clear();
  int type_index;
  for(type_index=0;type_index<number_of_cell_types;type_index++)
    {
      persistent_time.push_back(0.0);
      adhesion_sensitivity.push_back(0.0);
    };
  cell_type_system.get_persistent_times(persistent_time);
  cell_type_system.get_adhesion_sensitivities(adhesion_sensitivity);
  //
};
