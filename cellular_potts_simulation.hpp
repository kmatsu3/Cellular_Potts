#ifndef __CPM_SIMULATION__
#define __CPM_SIMULATION__
#include <vector>
#include <boost/random.hpp>
#include <cmath>
#include "cellular_potts_definition.hpp"
#include "cellular_potts_io.hpp"
class simulation_system_class {
  /*======================
    Members
   =======================*/
  //simulation setting parameter
private: long long int number_of_monte_carlo_steps_for_observation;
private: long long int number_of_monte_carlo_steps_for_initialization;
private: long long int number_of_monte_carlo_steps_for_relaxation;
private: long long int number_of_observations;
private: long long int period_of_observation;
private: int number_of_steps_in_single_step;
private: long long int number_of_flips;
private: std::string observation_output_format;
private: int number_of_control_parameters;
private: long int number_of_sweep_steps;
  //temporal simulation setting parameter
private: long long int number_of_monte_carlo_steps_in_simulation;
  // random number generator
public: boost::random::variate_generator<boost::mt19937,boost::random::uniform_01<> > random_generator_for_trial;
public: boost::random::variate_generator<boost::mt19937,boost::random::uniform_int_distribution<long long int> > random_generator_for_site_choice;
public: boost::random::variate_generator<boost::mt19937,boost::random::uniform_int_distribution<int> > random_generator_for_neighbor_choice;
public: boost::random::variate_generator<boost::mt19937,boost::random::uniform_01<> > random_generator_for_initializer;
  /*======================
    Methods
   =======================*/
public: long long int get_number_of_monte_carlo_steps(std::string procedure) const;
public: int get_number_of_steps_in_single_step() const;
public: void set_simulation();
private: void output_simulation_schedule() const;
public: std::vector<double> get_trial_random_number(long long int number_of_random_numbers);
public: std::vector<long long int> get_site_choice_random_number(long long int number_of_random_numbers);
public: std::vector<int> get_neighbor_choice_random_number(long long int number_of_random_numbers);
public: std::vector<double> get_initializer_random_number(long long int number_of_random_numbers);
public: long long int get_number_of_flips() const;
  //
public: long long int get_period_of_observation() const;
  //
public: long long int get_number_of_observations() const;
  //
public: std::string get_observation_output_format() const;
  // Methods for the schedule management of Monte Carlo simulations
public: void monte_carlo_setting(
				 const std::string procedure
				 );
  //
public: int get_number_of_control_parameters() const;
public: long int get_number_of_sweep_steps() const;
  //	      
  /*======================
    Constructor
   =======================*/
public: simulation_system_class(
				const model_parameters_cellular_potts_class                  & model,
				const boost::random::mt19937                                 & random_engine_for_trial,
				const boost::random::uniform_01<>                            & random_distribution_for_trial,
				const boost::random::mt19937                                 & random_engine_for_site_choice,
				const boost::random::uniform_int_distribution<long long int> & random_distribution_for_site_choice,
				const boost::random::mt19937                                 & random_engine_for_neighbor_choice,
				const boost::random::uniform_int_distribution<int>           & random_distribution_for_neighbor_choice,
				const boost::random::mt19937                                 & random_engine_for_initializer,
				const boost::random::uniform_01<>                            & random_distribution_for_initializer
				);
};
#endif // __CPM_SIMULATION___
