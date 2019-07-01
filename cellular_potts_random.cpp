#include "cellular_potts_random.hpp"
void random_system_seed_class::random_number_seed_input(
							const model_parameters_cellular_potts_class &model
							)
{
  io_cellular_potts io_method;
  {
    std::vector<unsigned long int> work_vector_unsignedlongint;
    // seed_for_trial
    io_method.get_input_unsignedlongint_array(
					      "model_input",
					      "random.seed_for_trial",
					      work_vector_unsignedlongint
					      );
    std::vector<unsigned long int>::iterator iterator_work=work_vector_unsignedlongint.begin();
    while(iterator_work!=work_vector_unsignedlongint.end())
      {
	random_for_trial_seeds.push_back(*iterator_work);
	iterator_work++;
      };
  };
  // seed for site choice
  {
    std::vector<unsigned long int> work_vector_unsignedlongint;
    // 
    io_method.get_input_unsignedlongint_array(
					      "model_input",
					      "random.seed_for_site_choice",
					      work_vector_unsignedlongint
					      );
    std::vector<unsigned long int>::iterator iterator_work=work_vector_unsignedlongint.begin();
    while(iterator_work!=work_vector_unsignedlongint.end())
      {
	random_for_site_choice_seeds.push_back(*iterator_work);
	iterator_work++;
      };
  };
  // seed for neighbor choice
  {
    std::vector<unsigned long int> work_vector_unsignedlongint;
    // 
    io_method.get_input_unsignedlongint_array(
					      "model_input",
					      "random.seed_for_neighbor_choice",
					      work_vector_unsignedlongint
					      );
    std::vector<unsigned long int>::iterator iterator_work=work_vector_unsignedlongint.begin();
    while(iterator_work!=work_vector_unsignedlongint.end())
      {
	random_for_neighbor_choice_seeds.push_back(*iterator_work);
	iterator_work++;
      };
  };
  // seed for initializer
  {
    std::vector<unsigned long int> work_vector_unsignedlongint;
    // 
    io_method.get_input_unsignedlongint_array(
					      "model_input",
					      "random.seed_for_initializer",
					      work_vector_unsignedlongint
					      );
    std::vector<unsigned long int>::iterator iterator_work=work_vector_unsignedlongint.begin();
    while(iterator_work!=work_vector_unsignedlongint.end())
      {
	random_for_initializer_seeds.push_back(*iterator_work);
	iterator_work++;
      };
  };
  //
};
//
const boost::random::seed_seq random_system_seed_class::get_seeds(
								  const std::string & flag
								  ) 
  const {
  boost::random::seed_seq return_value;
  if(flag=="trial")
    {
      boost::random::seed_seq work_seed_seq(random_for_trial_seeds.begin(), random_for_trial_seeds.end());
      return_value=work_seed_seq;
    } else if (flag=="site_choice")
    {
      boost::random::seed_seq work_seed_seq(random_for_site_choice_seeds.begin(), random_for_site_choice_seeds.end());
      return_value=work_seed_seq;
    } else if (flag=="neighbor_choice")
    {
      boost::random::seed_seq work_seed_seq(random_for_neighbor_choice_seeds.begin(), random_for_neighbor_choice_seeds.end());
      return_value=work_seed_seq;
    } else if (flag=="initializer")
    {
      boost::random::seed_seq work_seed_seq(random_for_initializer_seeds.begin(), random_for_initializer_seeds.end());
      return_value=work_seed_seq;
    } else 
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "random class",
			     "get seeds",
			     "due to irrigal flag."
			     );
    };
      return return_value;
};
		     
random_system_seed_class::random_system_seed_class(
						   const model_parameters_cellular_potts_class & model
						   )
{
  random_number_seed_input(model);
};
//
const boost::random::mt19937 random_system_setting_class::get_engine(
								     const std::string & flag
								     ) 
  const{
  boost::random::mt19937 return_value;
  if(flag=="trial")
    {
      return_value=random_for_trial_engine;
    } else if (flag=="site_choice")
    {
      return_value=random_for_site_choice_engine;
    } else if (flag=="neighbor_choice")
    {
      return_value=random_for_neighbor_choice_engine;
    } else if (flag=="initializer")
    {
      return_value=random_for_initializer_engine;
    } else 
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "random class",
			     "get engine",
			     "due to irrigal flag."
			     );
    };
      return return_value;
  
};
//
const boost::random::uniform_01<> random_system_setting_class::get_distribution_double(
										       const std::string & flag
										       ) 
  const{
      boost::random::uniform_01<> return_value;
  if(flag=="trial")
    {
      return_value=random_for_trial_distribution;
    } else if (flag=="initializer")
    {
      return_value=random_for_initializer_distribution;
    } else 
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "random class",
			     "get distribution double",
			     "due to irrigal flag."
			     );
    };
      return return_value;
  
};
//
const boost::random::uniform_int_distribution<long long int> random_system_setting_class::get_distribution_longlongint(
														       const std::string & flag
														       ) 
  const{
  boost::random::uniform_int_distribution<long long int> return_value;
  if(flag=="site_choice")
    {
      return_value=random_for_site_choice_distribution;
    } else 
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "random class",
			     "get distribution int",
			     "due to irrigal flag."
			     );
    };
      return return_value;
};
//
const boost::random::uniform_int_distribution<int> random_system_setting_class::get_distribution_int(
												     const std::string & flag
												     ) 
  const{
  boost::random::uniform_int_distribution<int> return_value;
  if(flag=="neighbor_choice")
    {
      return_value=random_for_neighbor_choice_distribution;
    } else 
    {
      io_cellular_potts io_method;
      io_method.error_output(
			     "random class",
			     "get distribution int",
			     "due to irrigal flag."
			     );
    };
      return return_value;
  
};
//
  /*======================
    Constructor
   =======================*/
random_system_setting_class::random_system_setting_class(
							 const model_parameters_cellular_potts_class & model,
							 const boost::random::seed_seq & random_for_trial_seeds,
							 const boost::random::seed_seq & random_for_site_choice_seeds,
							 const boost::random::seed_seq & random_for_neighbor_choice_seeds,
							 const boost::random::seed_seq & random_for_initializer_seeds,
							 const long long int & number_of_sites,
							 const int & number_of_neighbors
							 ): random_for_trial_engine(random_for_trial_seeds),
							    random_for_site_choice_engine(random_for_site_choice_seeds),
							    random_for_neighbor_choice_engine(random_for_neighbor_choice_seeds),
							    random_for_initializer_engine(random_for_initializer_seeds),
							    random_for_site_choice_distribution(0,number_of_sites-1),
							    random_for_neighbor_choice_distribution(0,number_of_neighbors-1)
{
};
