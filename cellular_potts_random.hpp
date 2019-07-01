#ifndef __CPM_RANDOM__
#define __CPM_RANDOM__
#include <vector>
#include <boost/random.hpp>
#include <cmath>
#include "cellular_potts_definition.hpp"
#include "cellular_potts_io.hpp"
class random_system_seed_class{
  /*======================
    Members
   =======================*/
private: std::vector<unsigned long int> random_for_initializer_seeds;
private: std::vector<unsigned long int> random_for_trial_seeds;
private: std::vector<unsigned long int> random_for_site_choice_seeds;
private: std::vector<unsigned long int> random_for_neighbor_choice_seeds;
  /*======================
    Methods
   =======================*/
  //
private: void random_number_seed_input(
				       const model_parameters_cellular_potts_class & model
				       );
  //
public: const boost::random::seed_seq get_seeds(
						const std::string & flag
						) const;
  //
  /*======================
    Constructor
   =======================*/
public:random_system_seed_class(
				const model_parameters_cellular_potts_class & model
				);
};
//
class random_system_setting_class{
  /*======================
    Members
    =======================*/
private: const boost::random::mt19937 random_for_trial_engine;
private: const boost::random::mt19937 random_for_site_choice_engine;
private: const boost::random::mt19937 random_for_neighbor_choice_engine;
private: const boost::random::mt19937 random_for_initializer_engine;
private: const boost::random::uniform_01<> random_for_trial_distribution;
private: const boost::random::uniform_int_distribution<long long int> random_for_site_choice_distribution;
private: const boost::random::uniform_int_distribution<int> random_for_neighbor_choice_distribution;
private: const boost::random::uniform_01<> random_for_initializer_distribution;
  /*======================
    Methods
    =======================*/
  //
public: const boost::random::mt19937 get_engine(
						const std::string & flag
						) const;
  //
public: const boost::random::uniform_01<> get_distribution_double(
								  const std::string & flag
								  ) const;
  //
public: const boost::random::uniform_int_distribution<long long int> get_distribution_longlongint(
												  const std::string & flag
												  ) const;
  //
public: const boost::random::uniform_int_distribution<int> get_distribution_int(
										const std::string & flag
										) const;
  //
  /*======================
    Constructor
   =======================*/
public:random_system_setting_class(
				   const model_parameters_cellular_potts_class & model,
				   const boost::random::seed_seq & random_for_trial_seeds,
				   const boost::random::seed_seq & random_for_site_choice_seeds,
				   const boost::random::seed_seq & random_for_neighbor_choice_seeds,
				   const boost::random::seed_seq & random_for_initializer_seeds,
				   const long long int & number_of_sites,
				   const int & number_of_neighbors
				   );
};
#endif // __CPM_RANDOM___
