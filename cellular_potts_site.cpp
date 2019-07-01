#include "cellular_potts_site.hpp"
#include <cstdlib>
#include <cmath>
#include <list>
//

  /*======================
    Constructor
   =======================*/
  site_cellular_potts_class::site_cellular_potts_class()
  {
  };
/*======================
  Generating a configuration for initalization
 =======================*/
/*

void site_system_class::get_configuration_setting(
						  )
  {
    io_cellular_potts io_method;
    generation_method=io_method.get_input_string(
						 "site_setting_input" ,
						 "configuration_setting.method"
						 );
    generation_type=io_method.get_input_string(
					       "site_setting_input" ,
					       "configuration_setting.type"
					       );
  };
void site_system_class::generate_configuration(
					       const model_parameters_cellular_potts_class & model,
					       const type_system_class & cell_type_system
					       )
  {
    if(generation_method=="generation")
      {
	if(generation_type=="default") generation_type="random_checkerboard";
	if(generation_type=="random_checkerboard") random_checkerboard(
								       model
								       );
      };
    if(generation_method=="input")
      {
      };
  };
*/
  /*======================
   Site coordinate initializar
   =======================*/
void site_system_class::generate_sites(
				       const model_parameters_cellular_potts_class & model
				       )
{
  int work_int;
  long long int site_index;
  set_coordinate_number(model);
  //
  std::vector<long int> coordinates;
  for(work_int=0;work_int<model.space_dimension;work_int++)
    {
      coordinates.push_back(0);
    };
  //
  std::vector<long long int> work_nearest_list;
  for(work_int=0;work_int<coordinate_number;work_int++)
    {
      work_nearest_list.push_back(0);
    };
  //
  //
  for (site_index=0;site_index<model.number_of_sites;site_index++)
    {
      site_cellular_potts_class work_site;
      //
      coordinate_generater(
			   model,
			   site_index,
			   coordinates
			   );
      //
      if(site_index!=coordinate_to_site(model,coordinates)) 
	{
	  io_cellular_potts io_method;
	  std::string message;
	  message = "site coordinate mismatching occurs at site";
	  message+= io_method.longlongint_to_string(site_index);
	  io_method.standard_output(message);
	  std::abort();
	};
      //
      work_site.set_coordinate(coordinates);
      //
      nearest_neighbor_generator(
				 model,
				 coordinates,
				 work_nearest_list
				 );
      //
      work_site.set_nearest_list(work_nearest_list,coordinate_number);
      //
      site.push_back(work_site);
      //
    };
  //
  neighbor_coordinate_generator(
				model
				);
  //
  show_neighbor_coordinates(
			    model
			    );
  //
  neighbor_sites_generator(
			   model
			   );
  //
  //  show_neighbor_sites(
  //	      model
  //	      );
};
//
void site_system_class::show_site_list(
				       const model_parameters_cellular_potts_class & model
				       )
{
  io_cellular_potts io_method;
  std::string message;
  int component_index;
  int site_index;
  int neighbor_index;
  for (site_index=0;site_index<model.number_of_sites;site_index++)
    {
      message = "Site[";
      message+= io_method.longlongint_to_string(site_index);
      message+= "] : Coordinate(";
      for (component_index=0;component_index<(int)site[site_index].get_coordinate_dimension();component_index++)
	{
	  message+= io_method.longint_to_string(site[site_index].get_coordinate_component(component_index));
	  if(component_index<(int)site[site_index].get_coordinate_dimension()-1) message+=",";
	};
      message+=") Neighbor Sites :: ";
      for(neighbor_index=0;neighbor_index<(int)site[site_index].get_nearest_neighbor_size();neighbor_index++)
	{
	  if((int)site[site_index].get_nearest_neighbor_component(neighbor_index)!=no_site_flag){
	    message+="D[";
		message+=io_method.int_to_string(neighbor_index);
		message+="]=";
		message+=io_method.longlongint_to_string(site[site_index].get_nearest_neighbor_component(neighbor_index));
		if(neighbor_index<(int)site[site_index].get_nearest_neighbor_size()-1) message+=",";
	  };
	};
      //	  message+= "\n";
	  //	  fprintf(stderr,"check out ?%s",message.c_str());
	  io_method.standard_output(message);
    };
};
//
void site_system_class::nearest_neighbor_generator(
						   const model_parameters_cellular_potts_class & model,
						   const std::vector<long int> & coordinates,
						   std::vector<long long int> & nearest_list
						   )
{
  int neighbor_index;
  std::vector<long int> neighbor_coordinates;
  long int work_longint;
  for(neighbor_index=0;neighbor_index<model.space_dimension;neighbor_index++)
    {
      neighbor_coordinates.push_back(0);
    };
  if((int)nearest_list.size()!=coordinate_number)
    {
	  io_cellular_potts io_method;
	  std::string message;
	  message = "Irrigal nearest_list initialization!";
	  message+= io_method.longint_to_string(nearest_list.size());
	  io_method.standard_output(message);
	  std::abort();
    }
  for(neighbor_index=0;neighbor_index<coordinate_number;neighbor_index++)
    {
      get_nearest_site_coordinate(
				  model,
				  coordinates,
				  neighbor_index,
				  neighbor_coordinates
				  );
      //
      work_longint=coordinate_to_site(
				      model,
				      neighbor_coordinates
				      );
      //					
      nearest_list[neighbor_index]=work_longint;
    };
};
//
void site_system_class::get_nearest_site_coordinate(
						    const model_parameters_cellular_potts_class & model,
						    const std::vector<long int> & coordinates,
						    const int & shift_index,
						    std::vector<long int> & shifted_coordinates
						    )
{
  std::vector <long int> shift_vector;
  std::vector <long int> work_vector;
  int direction_index;
  for (direction_index=0;direction_index<space_dimension;direction_index++)
    {
      shift_vector.push_back(0);
      work_vector.push_back(0);
    };
  //
  if((int)shifted_coordinates.size()<=space_dimension)
    {
      for(direction_index=shift_vector.size()+1;direction_index<model.space_dimension;direction_index++)
	{
	  shifted_coordinates.push_back(0);
	};
    };
  //
  nearest_shift_definer(
			model,
			shift_index,
			shift_vector
			);
  // Boundary condition is applied here.
  for (direction_index=0;direction_index<model.space_dimension;direction_index++)
    {
      if(model.get_boundary_condition(direction_index)=="periodic")
	{
	  work_vector[direction_index]=(
					shift_vector.at(direction_index)
					+coordinates.at(direction_index)
					+model.get_system_dimension(direction_index)
					)
	    %model.get_system_dimension(direction_index);
	}else{
	  work_vector[direction_index]=shift_vector.at(direction_index)
	    +coordinates.at(direction_index);
      };
    };
  // 
  shifted_coordinates=work_vector;
};
//
void site_system_class::set_coordinate_number(
					      const  model_parameters_cellular_potts_class & model
					      )
{
  std::string error_code="error";
  if(model.get_lattice_structure()=="square") 
    { 
      coordinate_number=2*model.space_dimension;
      error_code="square";
    }
  if(error_code=="error")
    {
	  io_cellular_potts io_method;
	  std::string message;
	  message = "Input lattice structure is undefined";
	  message+= model.get_lattice_structure();
	  io_method.standard_output(message);
	  std::abort();
    };
};
//
int site_system_class::get_coordinate_number(
					     const  model_parameters_cellular_potts_class & model
					     )
  const {
  return coordinate_number;
};
//
void site_system_class::nearest_shift_definer(
					      const model_parameters_cellular_potts_class & model,
					      const int & shift_index,
					      std::vector<long int> & shift_vector
					      )
{
  std::string error_code="error";
  if(model.get_lattice_structure()=="square")
    {
      //
      int direction_index;
      int work_direction_int;
      int work_sign_int;
      error_code="square";
      //
      if(shift_index<0||shift_index>coordinate_number-1)
	{
	  io_cellular_potts io_method;
	  std::string message;
	  message = "Input shift index is out of range";
	  message+= io_method.int_to_string(shift_index);
	  io_method.standard_output(message);
	  std::abort();
	}
      //
      if((int)shift_vector.size()!=model.get_space_dimension())
	{
	  io_cellular_potts io_method;
	  std::string message;
	  message = "Irrigal input shift vector is found";
	  message+= io_method.int_to_string(shift_index);
	  io_method.standard_output(message);
	  std::abort();
	};
      //
      work_sign_int=(shift_index%2)*2-1;
      work_direction_int=(shift_index-work_sign_int)/2;
      for (direction_index=0;direction_index<model.space_dimension;direction_index++)
	{
	  if(direction_index==work_direction_int) 
	    { 
	      shift_vector[model.space_dimension-direction_index-1]=work_sign_int;
	    }else{
	    shift_vector[model.space_dimension-direction_index-1]=0;
	  };
	};
    };
  if(error_code=="error")
    {
	  io_cellular_potts io_method;
	  std::string message;
	  message = "Input lattice structure is undefined";
	  message+= model.get_lattice_structure();
	  io_method.standard_output(message);
	  std::abort();
    };
};
//
long long int site_system_class::get_nearest_neighbor_site(
							   const model_parameters_cellular_potts_class & model,
							   const long long int & site_index,
							   const int & shift_direction
							   )
  const {
  return site[site_index].get_nearest_neighbor_component(shift_direction);
};
//
unsigned int site_cellular_potts_class::get_coordinate_dimension()
  const {
  return coordinates.size();
};
//
long int site_cellular_potts_class::get_coordinate_component( const int & component_index)
  const {
  return coordinates.at(component_index);
};
//
std::vector<long int> site_cellular_potts_class::get_coordinate_components() 
  const {
  return coordinates;
};
//
long long int site_cellular_potts_class::get_neighbor_site( const int & neighbor_index)
const {
  return neighbor_sites.at(neighbor_index);
};
//
int site_cellular_potts_class::get_region()
  const {
  return region;
};
//
/*
const std::vector<long long int> site_cellular_potts_class::get_neighbor_sites() 
  const{
  return neighbor_sites;
    };
*/
void site_cellular_potts_class::get_neighbor_sites(
						   std::vector<long long int> & neighbors
						   ) 
  const {
  copy(neighbor_sites.begin(), neighbor_sites.end(), neighbors.begin());
};
//
unsigned int site_cellular_potts_class::get_nearest_neighbor_size()
  const {
  return nearest_neighbor_sites.size();
};
//
long int site_cellular_potts_class::get_nearest_neighbor_component(const int & component_index)
  const {
  return nearest_neighbor_sites.at(component_index);
};
//
void site_cellular_potts_class::set_coordinate(
					       const std::vector<long int> & coordinate_vector
					       )
{
  std::vector<long int>::const_iterator coordinate_index=coordinate_vector.begin(); 
  if(coordinates.size()!=0)
    {
	  io_cellular_potts io_method;
	  std::string message;
	  message = "Coordinates on a site is already set.(site_cellular_potts_class::set_coordinate)";
	  io_method.standard_output(message);
	  std::abort();
    }; 
  while(coordinate_index!=coordinate_vector.end())
    {
      coordinates.push_back(*coordinate_index);
      coordinate_index++;
    };
};
//
void site_cellular_potts_class::set_nearest_list(
						 const std::vector<long long int> & nearest_list_vector,
						 const int & coordinate_number
						 )
{
  if(nearest_neighbor_sites.size()!=0)
    {
	  io_cellular_potts io_method;
	  std::string message;
	  message = "Neighbors on a site is already set.(site_cellular_potts_class::set_nearest_list)";
	  io_method.standard_output(message);
	  std::abort();
    }; 
  std::vector<long long int>::const_iterator neighbor_index=nearest_list_vector.begin(); 
  while(neighbor_index!=nearest_list_vector.end())
    {
      nearest_neighbor_sites.push_back(*neighbor_index);
      neighbor_index++;
    };
};
//
void site_cellular_potts_class::set_neighbor_sites(
						   const std::vector<long long int> & input_neighbor_sites
						   )
{
  neighbor_sites.clear();
  std::vector<long long int>::const_iterator neighbor_index=input_neighbor_sites.begin();
  while(neighbor_index!=input_neighbor_sites.end())
    {
      neighbor_sites.push_back(*neighbor_index);
      neighbor_index++;
    };
};
//
void site_cellular_potts_class::set_region(
					   const int & input_region
					   )
{
  region=input_region;
};
//
long long int site_system_class::coordinate_to_site(
						    const model_parameters_cellular_potts_class & model,
						    const std::vector<long int> & coordinates
						    )
const {
  int component_index;
  long long int dimension_size=1;
  int work_int;
  long long int return_value=0;
  std::string out_of_system="off";
  for(component_index=0;component_index<model.space_dimension;component_index++)
    {
      if(
	 (coordinates[component_index]>model.system_dimensions[component_index]-1)
	 ||(coordinates[component_index]<0)
	 ) out_of_system="on";
    };
  for(component_index=0;component_index<model.space_dimension;component_index++)
    {
      work_int=model.space_dimension-component_index-1;
      return_value += (long long int)coordinates.at(work_int)*dimension_size;
      dimension_size = dimension_size * (long long int)model.system_dimensions.at(work_int);
    };
  if (out_of_system=="on")
    {
      return no_site_flag;
    }else{
    return return_value;
  };
};
//
void site_system_class::coordinate_generater(
					     const model_parameters_cellular_potts_class & model,
					     const long long int & site_index,
					     std::vector<long int> & coordinates
					     )
  {
    // Working memory is defined
    int component_index;
    long int dimension;
    int work_int;
    long long int work_longlongint;
    std::vector<long int> work_vector;
    work_longlongint=site_index;
    //
    for(component_index=0;component_index<model.get_space_dimension();component_index++)
      {
	work_vector.push_back(0);
      };
    //
    if((int)coordinates.size()!=model.get_space_dimension())
      {
	io_cellular_potts io_method;
	std::string message;
	message = "coordinate dimension is not correctly set (site_cellular_potts_class::coordinate_generator";
	io_method.standard_output(message);
	std::abort();
      };
    //
    for(component_index=0;component_index<model.space_dimension;component_index++)
      {
	work_int=model.space_dimension-component_index-1;
	dimension=model.system_dimensions.at(work_int);
	work_vector[component_index]=(long int)(work_longlongint%(long long int)dimension);
	work_longlongint=(work_longlongint-(long long int)work_vector.at(component_index))/(long long int)dimension;
      };
    //    coordinates=work_vector;
    for(component_index=0;component_index<model.space_dimension;component_index++)
	  {
	    work_int=model.space_dimension-component_index-1;
	    coordinates[work_int]=work_vector.at(component_index);
	  };
  };
//
void site_system_class::sub_coordinate_generater(
						 const std::vector<long int> & sub_system_dimensions,
						 const long long int & site_index,
						 std::vector<long int> & coordinates
						 )
   const {
    // Working memory is defined
    int component_index;
    long int dimension;
    int work_int;
    long long int work_longlongint;
    std::vector<long int> work_vector;
    work_longlongint=site_index;
    int sub_space_dimension=sub_system_dimensions.size();
    //
    for(component_index=0;component_index<sub_space_dimension;component_index++)
      {
	work_vector.push_back(0);
      };
    //
    if((int)coordinates.size()!=sub_space_dimension)
      {
	io_cellular_potts io_method;
	std::string message;
	message = "coordinate dimension is not correctly set (site_cellular_potts_class::coordinate_generator)";
	io_method.standard_output(message);
	std::abort();
      };
    //
    for(component_index=0;component_index<sub_space_dimension;component_index++)
      {
	work_int=sub_space_dimension-component_index-1;
	dimension=sub_system_dimensions.at(work_int);
	work_vector[component_index]=(long int)(work_longlongint%(long long int)dimension);
	work_longlongint=(work_longlongint-(long long int)work_vector.at(component_index))/(long long int)dimension;
      };
    //    coordinates=work_vector;
    for(component_index=0;component_index<sub_space_dimension;component_index++)
	  {
	    work_int=sub_space_dimension-component_index-1;
	    coordinates[work_int]=work_vector.at(component_index);
	  };
  };
//
int site_system_class::get_plane_dimension() 
  const {
  return plane_dimension;
};
//
void site_system_class::plane_to_sub_coordinates(
						 const long int & plane_index, 
						 const model_parameters_cellular_potts_class & model,
						 const std::vector<std::string> & plane_identifier,
						 std::vector<long int> & plane_coordinates
						 )
  const {
  int direction_index;
  int space_dimension=model.get_space_dimension();
  std::vector<long int> system_dimensions;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      system_dimensions.push_back(model.get_system_dimension(direction_index));
    };
  /*
  if(space_dimension<3){
    io_cellular_potts io_method;
    io_method.error_output(
			   "site class",
			   "plane_to_sub_coordinates",
			   "due to space dimension lower than 2."
			   );
  }
  */
  std::vector<long int> sub_system_dimensions;
  for(direction_index=0;direction_index<space_dimension;direction_index++)
    {
      if(plane_identifier[direction_index]=="off"){
	sub_system_dimensions.push_back(system_dimensions[direction_index]);
      };
    };
  sub_coordinate_generater(
			   sub_system_dimensions,
			   plane_index,
			   plane_coordinates
			   );
};
//
/*
  const void site_system_class::get_coordinates(
						const model_parameters_cellular_potts_class & model,
						const long long int & site_index,
						std::vector<long int> & coordinates
						)
  const {
    int component_index;
    int space_dimension=model.get_space_dimension();
    if((int)coordinates.size()==space_dimension)
      {
	for(component_index=0;component_index<space_dimension;component_index++)
	  {
	    coordinates[component_index]=site[site_index].get_coordinate_component(component_index);
	  };
      }else{
      io_cellular_potts io_method;
      io_method.error_output(
			     "site_system_class",
			     "get_coordinates",
			     "Dimension of coordinates is less than space dimension."
			     );
    };
  };
*/
std::vector<long int> site_system_class::get_coordinates(
							 const long long int & site_index
							 )
  const {
  return site[site_index].get_coordinate_components();
};
//
/*
std::vector<long int> site_system_class::get_site_coordinates(
							      const model_parameters_cellular_potts_class & model,
							      const long long int & site_index
							      )
  const {
  int component_index;
  int space_dimension=model.get_space_dimension();
  std::vector<long int> work_vector(space_dimension);
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      work_vector[component_index]=site[site_index].get_coordinate_component(component_index);
    };
  return work_vector;
};
*/
//
void site_system_class::get_site_coordinates(
					     const long long int & site_index,
					     std::vector<long int> & coordinates
					     )
  const {
  int component_index;
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      coordinates[component_index]=site[site_index].get_coordinate_component(component_index);
    };
};
//
int site_system_class::get_number_of_neighbor_sites() 
  const{
  if(number_of_neighbor_sites!=-1)
    {
      return number_of_neighbor_sites;
    } else {
    io_cellular_potts io_method;
    std::string message;
    message = "Uncalculated the number of neighbor sites. (site_system_class::get_number_of_neighbor_sites)";
    io_method.standard_output(message);
    std::abort();
  };
};
void site_system_class::neighbor_coordinate_generator(
						      const model_parameters_cellular_potts_class & model
						      )
{
  std::string lattice_structure=model.get_lattice_structure();
  std::string neighbor_definition=model.get_neighbor_definition();
  int space_dimension=model.get_space_dimension();
  int interaction_depth=model.get_interaction_depth();
  toolbox tool;
  neighbor_coordinates.clear();
  if((lattice_structure=="square")&&(neighbor_definition=="off_site"))
    {
      int neighbor_index;
      std::list<double> distance;
      int number_of_candidates=1;
      // caluclulation of cut off
      for(neighbor_index=0;neighbor_index<space_dimension;neighbor_index++)
	{
	  number_of_candidates=number_of_candidates*interaction_depth;
	};
      std::vector<long int> shift_vector;
      for(neighbor_index=0;neighbor_index<space_dimension;neighbor_index++)
	{
	  shift_vector.push_back(-1);
	};
      int counter=0;
      double candidate_distance;
      for(neighbor_index=0;neighbor_index<number_of_candidates;neighbor_index++)
	{
	  candidate_distance=neighbor_shift_distance(
						     model,
						     "positive",
						     neighbor_index,
						     shift_vector
						     );
	  //
	  if(tool.existence_list_double(distance,candidate_distance)==0)
	    {
	      distance.push_back(candidate_distance);
	      counter++;
	    };
	};
      distance.sort();
      double cut_off=tool.extend_little(tool.get_element_of_double_list(distance,interaction_depth));
      //
      number_of_candidates=1;
      for(neighbor_index=0;neighbor_index<space_dimension;neighbor_index++)
	{
	  number_of_candidates=number_of_candidates*(2*interaction_depth+1);
	};
      counter=0;
      for(neighbor_index=0;neighbor_index<number_of_candidates;neighbor_index++)
	{
	  candidate_distance=neighbor_shift_distance(
						     model,
						     "isotropic",
						     neighbor_index,
						     shift_vector
						     );
	  if((candidate_distance<cut_off)&&(tool.finite_abs_check(candidate_distance)==1))
	    {
	      neighbor_coordinates.push_back(shift_vector);
	      counter++;
	    };
	};
      number_of_neighbor_sites=counter;
      //
    } else {
    //
    io_cellular_potts io_method;
    std::string message;
    message = "Input lattice structure or neighbor_difinition is undefined (site_system_class::neighbor_coordinate_generator)";
    message+= lattice_structure;
    io_method.standard_output(message);
    std::abort();
  };
};
//
std::vector<long int> site_system_class::get_neighbor_coordinates(
								  const int & neighbor_index
								  )
  const {
  return neighbor_coordinates[neighbor_index];
};
//
void site_system_class::show_neighbor_coordinates(
						  const model_parameters_cellular_potts_class & model
						  )
  const{
    io_cellular_potts io_method;
    std::string message;
    int neighbor_index;
    int component_index;
    int space_dimension=model.get_space_dimension();
    //
    //    fprintf(stderr,"ok?%d\n",number_of_neighbor_sites);
    for(neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
      {
	//	fprintf(stderr,"%d\n",neighbor_index);
	message = "Neighbor site[";
	message+= io_method.int_to_string(neighbor_index);
	message+= "]=(";
	for(component_index=0;component_index<space_dimension;component_index++)
	  {
	    message+=io_method.double_to_string(neighbor_coordinates[neighbor_index][component_index]);
	    if(component_index<space_dimension-1) message+=",";
	  };
	message+=")";
	io_method.standard_output(message);
      };
};
//
void site_system_class::neighbor_sites_generator(
						 const model_parameters_cellular_potts_class & model
						 )
{
  long long int numbr_of_sites=model.get_number_of_sites();
  long long int site_index;
  int neighbor_index;
  std::vector<long long int> work_neighbor_sites(number_of_neighbor_sites,-1);
  int component_index;
  int space_dimension=model.get_space_dimension();
  std::vector<long int> coordinates(space_dimension,-1);
  std::vector<long int> shifted_coordinates(space_dimension,-1);
  std::vector<long int> work_coordinates(space_dimension,-1);
  for(site_index=0;site_index<numbr_of_sites;site_index++)
    {
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  coordinates[component_index]=site[site_index].get_coordinate_component(component_index);
	};
      //
      for(neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
	{
	  for(component_index=0;component_index<space_dimension;component_index++)
	    {
	      shifted_coordinates[component_index]=coordinates[component_index]
		+neighbor_coordinates[neighbor_index][component_index];
	    };
	  //
	  boundary_trimming(
			    model,
			    shifted_coordinates,
			    work_coordinates
			    );
	  //
	  work_neighbor_sites[neighbor_index]=coordinate_to_site(
								 model,
								 work_coordinates
								 );
	  //
	  if(work_neighbor_sites[neighbor_index]==no_site_flag)
	    {
	      work_neighbor_sites[neighbor_index]=site_index;
	    };
	};
      site[site_index].set_neighbor_sites(work_neighbor_sites);
    };
};
//
/*
const std::vector<long long int> site_system_class::get_neighbor_sites(
								       const long long int site_index
								       ) 
  const{
  return site[site_index].get_neighbor_sites();
    };
*/
void site_system_class::get_neighbor_sites(
					   const long long int & site_index,
					   std::vector<long long int> & neighbor_sites
					   ) 
  const{
  site[site_index].get_neighbor_sites(neighbor_sites);
};
//
long long int site_system_class::get_neighbor_site(
						   const long long int site_index,
						   const int neighbor_index
						   ) 
  const{
  return site[site_index].get_neighbor_site(neighbor_index);
    };
//
void site_system_class::boundary_trimming(
					  const model_parameters_cellular_potts_class & model,
					  const std::vector<long int> & input_coordinates,
					  std::vector<long int> & output_coordinates
					  ) 
  const {
  int direction_index;
  for (direction_index=0;direction_index<model.get_space_dimension();direction_index++)
    {
      if(model.get_boundary_condition(direction_index)=="periodic")
	{
	  output_coordinates[direction_index]=(
					       input_coordinates.at(direction_index)
					       +model.get_system_dimension(direction_index)
					       )
	    %model.get_system_dimension(direction_index);
	}else
	{
	  output_coordinates[direction_index]=input_coordinates.at(direction_index);
	};
    };
};
//
/* old shit function
const std::vector<long int> site_system_class::origin_shift(
							    const model_parameters_cellular_potts_class & model,
							    const long long int & origin_site_index,
							    const long long int & site_index
							    )
  {
  int component_index;
  //  origin_position=site[origin_site_index].get_coordinate_components();
  //  present_position=site[site_index].get_coordinate_components();
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      origin_position[component_index]=site[origin_site_index].get_coordinate_component(component_index);
      present_position[component_index]=site[site_index].get_coordinate_component(component_index);
    };
  //
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      if(boundary_conditions[component_index]=="periodic")
	{
	  double criteria_value;
	  criteria_value=(double)(present_position.at(component_index)-origin_position.at(component_index))/(double)system_dimensions.at(component_index)*2.0;
	    if(criteria_value<-1.0)
	      {
		return_vector[component_index]=present_position[component_index]-origin_position[component_index]+system_dimensions[component_index];
	      } else if(criteria_value>1.0)
	      {
		return_vector[component_index]=present_position[component_index]-origin_position[component_index]-system_dimensions[component_index];
	      } else 
	      {
		return_vector[component_index]=present_position[component_index]-origin_position[component_index];
	      };
	} 
      else
	{
	  return_vector[component_index]=present_position[component_index]-origin_position[component_index];
	};
    };
  // debug
  // fprintf(stderr,"debug:ok?");
  {
  io_cellular_potts io_method;
  std::string message;
  message ="site=";
  message+=io_method.longlongint_to_string(site_index);
  message+="/";
  message+=io_method.longintarray_to_string(position);
  message+="::";
  message+="origin=";
  message+=io_method.longlongint_to_string(origin_site_index);
  message+="/";
  message+=io_method.longintarray_to_string(origin_position);
  message+="::";
  message+="relative";
  message+="/";
  message+=io_method.longintarray_to_string(return_vector);
  io_method.standard_output(message);
  };
  // debug
  return return_vector;
};
*/
void site_system_class::origin_shift(
				     const model_parameters_cellular_potts_class & model,
				     const long long int & origin_site_index,
				     const long long int & site_index,
				     std::vector<long int> & return_vector
				     )
  {
  int component_index;
  //  origin_position=site[origin_site_index].get_coordinate_components();
  //  present_position=site[site_index].get_coordinate_components();
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      origin_position[component_index]=site[origin_site_index].get_coordinate_component(component_index);
      present_position[component_index]=site[site_index].get_coordinate_component(component_index);
    };
  //
  for(component_index=0;component_index<space_dimension;component_index++)
    {
      if(boundary_conditions[component_index]=="periodic")
	{
	  double criteria_value;
	  criteria_value=(double)(present_position[component_index]-origin_position[component_index])/(double)system_dimensions[component_index]*2.0;
	    if(criteria_value<-1.0)
	      {
		return_vector[component_index]=present_position[component_index]-origin_position[component_index]+system_dimensions[component_index];
	      } else if(criteria_value>1.0)
	      {
		return_vector[component_index]=present_position[component_index]-origin_position[component_index]-system_dimensions[component_index];
	      } else 
	      {
		return_vector[component_index]=present_position[component_index]-origin_position[component_index];
	      };
	} 
      else
	{
	  return_vector[component_index]=present_position[component_index]-origin_position[component_index];
	};
    };
};
//
void site_system_class::show_neighbor_sites(
					    const model_parameters_cellular_potts_class & model
					    ) 
  const{
  io_cellular_potts io_method;
  io_method.standard_output("==== Initial Neighbor site list ====");
  std::string message;
  long long int site_index;
  long long int number_of_sites=model.get_number_of_sites();
  int neighbor_index;
  for(site_index=0;site_index<number_of_sites;site_index++)
    {
      message ="Site[";
      message+=io_method.int_to_string(site_index);
      message+="]=(";
	for (neighbor_index=0;neighbor_index<number_of_neighbor_sites;neighbor_index++)
	  {
	    message+=io_method.longlongint_to_string(site[site_index].get_neighbor_site(neighbor_index));
	    if(neighbor_index<number_of_neighbor_sites-1) message+=",";
	  };
      message+=")";
      io_method.standard_output(message);
    }
};
//
double site_system_class::neighbor_shift_distance(
						  const model_parameters_cellular_potts_class & model,
						  const std::string & bound,
						  const int & neighbor_index,
						  std::vector<long int> & shift_vector
						  )
const {
  std::string lattice_structure=model.get_lattice_structure();
  //  int space_dimension=model.get_space_dimension();
  int interaction_depth=model.get_interaction_depth();
  if((int)shift_vector.size()!=space_dimension)
    {
      io_cellular_potts io_method;
      std::string message;
      message = "Irrigal dimension of shift_vector size (site_system_class::neighbor_shift_calculater)";
      message+= "Lattiec structure =" + io_method.int_to_string((int)shift_vector.size());
      io_method.standard_output(message);
      std::abort();
    };
  if((lattice_structure=="square")&&(bound=="positive"))
    {
      long int component_index;
      long int dimension=(long int)interaction_depth;
      long int work_int=-1;
      long int work_int_2=neighbor_index;
      double distance=0.0;
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  work_int=space_dimension-component_index-1;
	  shift_vector[work_int]=(long int)(work_int_2%dimension);
	  work_int_2=(work_int_2-shift_vector[work_int])/dimension;
	  distance=distance+(double)(shift_vector[work_int]*shift_vector[work_int]);
	};
      distance=pow(distance,0.5);
      return distance;
    }else if((lattice_structure=="square")&&(bound=="isotropic")){
      long int component_index;
      long int dimension=(long int)interaction_depth*2+1;
      long int work_int=-1;
      long int work_int_2=neighbor_index;
      double distance=0.0;
      for(component_index=0;component_index<space_dimension;component_index++)
	{
	  work_int=space_dimension-component_index-1;
	  shift_vector[work_int]=(long int)(work_int_2%dimension);
	  work_int_2=(work_int_2-shift_vector[work_int])/dimension;
	  shift_vector[work_int]=shift_vector[work_int]-interaction_depth;
	  distance=distance+(double)(shift_vector[work_int]*shift_vector[work_int]);
	};
      distance=pow(distance,0.5);
      return distance;
    }else{
    //
    io_cellular_potts io_method;
    std::string message;
    message = "Input lattice structure or bound is undefined (sit**e_system_class::neighbor_shift_calculater)";
    message+= "Lattiec structure =" + lattice_structure;
    message+= "/Bound condition =" + bound + ".";
    io_method.standard_output(message);
    std::abort();
  };
};
  /*======================
    Site system Constructor
   =======================*/
site_system_class::site_system_class(
				     const model_parameters_cellular_potts_class & model
				     )
{
    generation_method="generation";
    generation_type="default";
    generation_type_default="random_checkerboard";
    plane_dimension=2;
    number_of_neighbor_sites=-1;
    coordinate_number=default_coordinate_number;
    //
    space_dimension=model.get_space_dimension();
    int direction_index;
    for(direction_index=0;direction_index<space_dimension;direction_index++)
      {
	boundary_conditions.push_back("");
        system_dimensions.push_back(0);
        origin_position.push_back(0);
	present_position.push_back(0);
	//        return_vector.push_back(0);
      }
    boundary_conditions=model.get_boundary_conditions();
    system_dimensions=model.get_system_dimensions();
    //
};
void site_system_class::set_structure_site_item(
						int site_index,
						std::string child
						)
{
  io_cellular_potts io_method;
  structure_item=io_method.generate_structure(
					      "site",
					      site_index,
					      child
					      );
};

