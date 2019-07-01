#ifndef __MODEL_PARAMETERS_CELLULAR_POTTS__
#define __MODEL_PARAMETERS_CELLULAR_POTTS__
#include <string>
#include <vector>
class model_parameters_cellular_potts_class
{
  /*
    Definitions of model parameters
   */
public:long int number_of_cells;
public:long int number_of_cell_types;
public:int number_of_adhesion;
public:int number_of_adhesion_components;
public:long long int number_of_sites;
public:int space_dimension;
public:std::vector<long int> system_dimensions;
private:std::vector<std::string> boundary_conditions;
private:std::string lattice_structure;
private:int interaction_depth;
private:std::string neighbor_definition;
private:std::string cell_tracking_flag;
private:static const int space_dimension_upper_bound;
private:static const int space_dimension_lower_bound;
private:int number_of_regions;
// system_dimensions[space_dimension] should be corrected for any space dimension.
  /*
    Function declaration
   */
public: int show_model_parameters();
public: int input_model_parameters();
public: long long int evaluate_number_of_sites();
public: long long int get_number_of_sites () const;
public: long int get_number_of_cell_types () const;
public: long int get_number_of_cells () const;
public: int get_number_of_adhesion () const;
public: int get_number_of_adhesion_components () const;
public: int get_number_of_regions () const;
public: int get_space_dimension () const;
public: long int get_system_dimension (const int & direction_index) const;
public: std::vector<long int> get_system_dimensions () const;
public: std::string get_boundary_condition(const int direction_index) const;
public: std::vector<std::string> get_boundary_conditions() const;
public: std::string get_lattice_structure () const;
public: int get_interaction_depth () const;
public: std::string get_neighbor_definition () const;
public: std::string get_cell_tracking_flag () const;
public: std::string check_out_of_range(
					     const std::vector<long int> coordinates
					     )const;
private:
  void check_model_parameter(std::string model_parameter);
  /*
    Constructer declaration
  */
public:
  model_parameters_cellular_potts_class();
  /*
    Local variables declaration
  */
private:
};
#endif // __MODEL_PARAMETERS_CELLULAR_POTTS__
