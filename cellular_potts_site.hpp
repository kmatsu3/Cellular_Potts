#ifndef __SITE__
#define __SITE__
#include <vector>
#include "toolbox.hpp"
#include "cellular_potts_definition.hpp"
// don't include cellular_potts_cell.hpp for saving dependency
#include "cellular_potts_type.hpp"
#include "cellular_potts_io.hpp"
//
class site_cellular_potts_class
{
  /*======================
    Members
   =======================*/
private: std::vector<long int> coordinates;
private: std::vector<long long int> nearest_neighbor_sites;
private: std::vector<long long int> neighbor_sites;
private: int region;
  /*======================
    Methods
   =======================*/
  //
public: void set_coordinate(
			    const std::vector<long int> & coordinate_vector
			    );
  //
public: void set_nearest_list(
			      const std::vector<long long int> & nearest_list_vector,
			      const int & coordinate_number
			      );
public: void set_neighbor_sites(
				const std::vector<long long int> & input_neighbor_sites
				);
  //
public: void set_region(
			const int & input_region
			);
  //
public: unsigned int get_coordinate_dimension() const;
public: long int get_coordinate_component(const int & component_index) const;
public: std::vector<long int> get_coordinate_components() const;
public: unsigned int get_nearest_neighbor_size() const;
public: long int get_nearest_neighbor_component(const int & component_index) const;
public: long long int get_neighbor_site(const int & neighbor_index) const;
public: int get_region() const;
  // public: const std::vector<long long int> get_neighbor_sites() const;
public: void get_neighbor_sites(std::vector<long long int> & neighbors ) const;
  /*======================
    Constructor
   =======================*/
public: site_cellular_potts_class();
};
//
class site_system_class{
  /*======================
    Members
   =======================*/
private: std::vector<site_cellular_potts_class> site;
private: std::string generation_method;
private: std::string generation_type;
private: std::string generation_type_default;
private: int plane_dimension;
private: int number_of_neighbor_sites;
private: std::vector< std::vector<long int> > neighbor_coordinates;
private: int space_dimension;
  /*======================
    Temporary Members for methods
   =======================*/
private: std::string structure_item;
private: int coordinate_number;
public: static const long long int no_site_flag=-1; 
public: static const int default_coordinate_number=-1;
  /*======================
    Work Memories
   =======================*/
private:std::vector<std::string> boundary_conditions;
private:std::vector<long int> system_dimensions;
private:std::vector<long int> origin_position;
private:std::vector<long int> present_position;
private:std::vector<long int> return_vector;
  /*======================
    Methods
   =======================*/
public: site_system_class(
			  const model_parameters_cellular_potts_class & model
			  );
  //
public: void coordinate_generater(
				  const model_parameters_cellular_potts_class & model,
				  const long long int & site_index,
				  std::vector<long int> & coordinates
				  );
  //
  /*
public: const void get_coordinates(
				   const model_parameters_cellular_potts_class & model,
				   const long long int & site_index,
				   std::vector<long int> & coordinates
				   ) const;
  */
  //
public: std::vector<long int> get_coordinates(
					      const long long int & site_index
					      ) const;
  //
  /*
public: std::vector<long int> get_site_coordinates(
						   const model_parameters_cellular_potts_class &model,
						   const long long int & site_index
						   ) const;
  */
public: void get_site_coordinates(
				  const long long int & site_index,
				  std::vector<long int> & coordinates
				  ) const;
  //
public: int get_number_of_neighbor_sites() const;
  //
public: long long int coordinate_to_site(
					 const model_parameters_cellular_potts_class & model,
					 const std::vector<long int> & coordiates
					 ) const;
  //
public: void sub_coordinate_generater(
				      const std::vector<long int> & sub_system_dimensions,
				      const long long int & site_index,
				      std::vector<long int> & coordinates
				      ) const;
  //
public: int get_plane_dimension() const;
  //
public: void plane_to_sub_coordinates(
				      const long int & plane_index, 
				      const model_parameters_cellular_potts_class & model,
				      const std::vector<std::string> & plane_identifier,
				      std::vector<long int> & coordiates
				      ) const;			
  //
public: long long int get_nearest_neighbor_site(
						const model_parameters_cellular_potts_class & model,
						const long long int & site_index,
						const int & shift_direction
						) const;
  //
public: void nearest_neighbor_generator(
					const model_parameters_cellular_potts_class & model,
					const std::vector<long int> & coordinates,
					std::vector<long long int> & nearest_list
					);
  //
public: void neighbor_coordinate_generator(
					   const model_parameters_cellular_potts_class & model
					   );
  //
private: double neighbor_shift_distance(
					const model_parameters_cellular_potts_class & model,
					const std::string & bound,
					const int & neighbor_index,
					std::vector<long int> & shift_vector
					) const;
  //
public:std::vector<long int> get_neighbor_coordinates(
						      const int & neighbor_index
						      ) const;
  //
public:void show_neighbor_coordinates(
				      const model_parameters_cellular_potts_class & model
				      )const ;	
  //			   
public:void neighbor_sites_generator(
				     const model_parameters_cellular_potts_class & model
				     );
  //
  /*
public:const std::vector<long long int> get_neighbor_sites(
							   const long long int site_index
							   ) const;
  */
public:void get_neighbor_sites(
			       const long long int & site_index,
			       std::vector<long long int> & neighbor_sites
			       ) const;
  //
public:long long int get_neighbor_site(
				       const long long int site_index,
				       const int neighbor_index
				       ) const;
  //
public:void boundary_trimming(
			      const model_parameters_cellular_potts_class & model,
			      const std::vector<long int> & coordinates,
			      std::vector<long int> & shifted_coordinates
			      ) const;
  //
public:void origin_shift(
			 const model_parameters_cellular_potts_class & model,
			 const long long int & origin_site_index,
			 const long long int & site_index,
			 std::vector<long int> & return_vector
			 );
  //
public:void show_neighbor_sites(
				const model_parameters_cellular_potts_class & model
				) const;
  //
public: void generate_sites(
			    const model_parameters_cellular_potts_class & model
			    );
  //
public: void show_site_list(
			    const model_parameters_cellular_potts_class & model
			    );
  //
private: void set_structure_site_item(
				      int site_index,
				      std::string child
				      );
  //
private: void nearest_shift_definer(
				    const model_parameters_cellular_potts_class & model,
				    const int & shift_index,
				    std::vector<long int> & shift_vector
				    );
  //
private: void get_nearest_site_coordinate(
					  const  model_parameters_cellular_potts_class & model,
					  const std::vector<long int> & coordinates,
					  const int & shift_index,
					  std::vector<long int> & shifted_coordinates
					  );
  //
private: void set_coordinate_number(
				    const  model_parameters_cellular_potts_class & model
				    );
  //
public: int get_coordinate_number(
				  const  model_parameters_cellular_potts_class & model
				  ) const;
};
#endif //__SITE__
