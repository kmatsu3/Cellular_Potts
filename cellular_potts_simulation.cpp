#include "cellular_potts_simulation.hpp"
  /*======================
    Methods
   =======================*/
long long int simulation_system_class::get_number_of_monte_carlo_steps(std::string procedure) 
  const {
  long long int return_value=-1;
  if(procedure=="initialization")
    {
      return_value=number_of_monte_carlo_steps_for_initialization; 
    } 
  else if(procedure=="relaxation")
    {
      return_value=number_of_monte_carlo_steps_for_relaxation;
    } 
  else if(procedure=="observation")
    {
      return_value=number_of_monte_carlo_steps_for_observation;
    }
  else if(procedure=="temporal")
    {
      return_value=number_of_monte_carlo_steps_in_simulation;
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "simulation_system_class",
			     "get_number_of_monte_carlo_steps",
			     "due to undefined monte carlo procedure"
			     );      
    };
  return return_value;
};
int simulation_system_class::get_number_of_steps_in_single_step()
  const {
  return number_of_steps_in_single_step;
};
//
void simulation_system_class::set_simulation()
{
  io_cellular_potts io_method;
  number_of_monte_carlo_steps_for_observation=io_method.get_input_longlongint(
									      "model_input",
									      "simulation.number_of_mc_steps.observation"
									      );
  number_of_monte_carlo_steps_for_relaxation=io_method.get_input_longlongint(
									     "model_input",
									     "simulation.number_of_mc_steps.relaxation"
									     );
  number_of_monte_carlo_steps_for_initialization=io_method.get_input_longlongint(
										 "model_input",
										 "simulation.number_of_mc_steps.initialization"
										 );
  number_of_steps_in_single_step=io_method.get_input_int(
							 "model_input",
							 "simulation.number_of_steps_in_mc_step"
							 );
  number_of_observations=io_method.get_input_longlongint(
							 "model_input",
							 "simulation.number_of_observations"
							 );
  //
  number_of_control_parameters=io_method.get_input_int(
							"model_input",
							"simulation.number_of_control_parameters"
							);
  //
  number_of_sweep_steps=io_method.get_input_int(
						 "model_input",
						 "simulation.number_of_sweep_steps"
						 );
  //
  if(number_of_observations>0)
    {
      period_of_observation
	=number_of_monte_carlo_steps_for_observation
	/number_of_observations;
    }
  if(period_of_observation<1) 
    {
      period_of_observation=1;
    };
  observation_output_format=io_method.get_input_string(
						       "model_input",
						       "simulation.observation_output_format"
						       );
  output_simulation_schedule();
};
//
void simulation_system_class::output_simulation_schedule()
const {
  io_cellular_potts io_method;
  std::string message;
  io_method.standard_output("=== Simulation Schedule ===");
  message = "initialization mcs:" + io_method.longlongint_to_string(number_of_monte_carlo_steps_for_initialization);
  io_method.standard_output(message);
  message = "relaxation mcs:" + io_method.longlongint_to_string(number_of_monte_carlo_steps_for_relaxation);
  io_method.standard_output(message);
  message = "obsrvation mcs:" + io_method.longlongint_to_string(number_of_monte_carlo_steps_for_observation);
  io_method.standard_output(message);
  message = "times of obsrvation:" + io_method.longlongint_to_string(number_of_observations);
  io_method.standard_output(message);
}
//
std::vector<double> simulation_system_class::get_trial_random_number(long long int number_of_random_numbers) 
{
  std::vector<double> work_vector(number_of_random_numbers,0);
  long long int random_index;
  for(random_index=0;random_index<number_of_random_numbers;random_index++)
    {
      work_vector[random_index]=random_generator_for_trial();
    };
  return work_vector;
};
//
std::vector<double> simulation_system_class::get_initializer_random_number(long long int number_of_random_numbers)
 {
   std::vector<double> work_vector(number_of_random_numbers,0);
   long long int random_index;
   for(random_index=0;random_index<number_of_random_numbers;random_index++)
     {
       work_vector[random_index]=random_generator_for_initializer();
     };
   return work_vector;
};
//
std::vector<long long int> simulation_system_class::get_site_choice_random_number(long long int number_of_random_numbers) 
{
  std::vector<long long int> work_vector(number_of_random_numbers,0);
  long long int random_index;
  for(random_index=0;random_index<number_of_random_numbers;random_index++)
    {
      work_vector[random_index]=random_generator_for_site_choice();
    };
  return work_vector;
};
//
std::vector<int> simulation_system_class::get_neighbor_choice_random_number(long long int number_of_random_numbers) 
{
  std::vector<int> work_vector(number_of_random_numbers,0);
  long long int random_index;
  for(random_index=0;random_index<number_of_random_numbers;random_index++)
    {
      work_vector[random_index]=random_generator_for_neighbor_choice();
    };
  return work_vector;
};
//
long long int simulation_system_class::get_number_of_flips()
  const {
  return number_of_flips;
};
//
long long int simulation_system_class::get_number_of_observations()
  const {
  return number_of_observations;
};
//
long long int simulation_system_class::get_period_of_observation()
  const {
  return period_of_observation;
};
//
std::string simulation_system_class::get_observation_output_format() 
  const {
  return observation_output_format;
};
//
int simulation_system_class::get_number_of_control_parameters() 
  const {
  return number_of_control_parameters;
};
//
long int simulation_system_class::get_number_of_sweep_steps() 
  const {
  return number_of_sweep_steps;
};
//
void simulation_system_class::monte_carlo_setting(
						  const std::string procedure
						  )
{
  if(procedure=="initialization")
    {
      number_of_monte_carlo_steps_in_simulation=number_of_monte_carlo_steps_for_initialization; 
    } 
  else if(procedure=="relaxation")
    {
      number_of_monte_carlo_steps_in_simulation=number_of_monte_carlo_steps_for_relaxation;
    } 
  else if(procedure=="observation")
    {
      number_of_monte_carlo_steps_in_simulation=number_of_monte_carlo_steps_for_observation;
    }
  else
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "simulation_system_class",
			     "monte_carlo_setting",
			     "due to undefined monte carlo procedure"
			     );      
    };
};
//
  /*======================
    Constructor
   =======================*/
simulation_system_class::simulation_system_class(
						 const model_parameters_cellular_potts_class & model,
						 const boost::random::mt19937 & random_engine_for_trial,
						 const boost::random::uniform_01<> & random_distribution_for_trial,
						 const boost::random::mt19937 & random_engine_for_site_choice,
						 const boost::random::uniform_int_distribution<long long int> & random_distribution_for_site_choice,
						 const boost::random::mt19937 & random_engine_for_neighbor_choice,
						 const boost::random::uniform_int_distribution<int> & random_distribution_for_neighbor_choice,
						 const boost::random::mt19937 & random_engine_for_initializer,
						 const boost::random::uniform_01<> & random_distribution_for_initializer
						 ):  random_generator_for_trial(random_engine_for_trial,
										random_distribution_for_trial),
						     random_generator_for_site_choice(random_engine_for_site_choice,
										      random_distribution_for_site_choice),
						     random_generator_for_neighbor_choice(random_engine_for_neighbor_choice,
											  random_distribution_for_neighbor_choice),
						     random_generator_for_initializer(random_engine_for_initializer,
										      random_distribution_for_initializer)
{
  period_of_observation=0;
  set_simulation();
  long long int number_of_sites=model.get_number_of_sites();
  number_of_flips=number_of_sites*(long long int)number_of_steps_in_single_step;
};
