#ifndef H_NAMED_FUNC_UTILITIES
#define H_NAMED_FUNC_UTILITIES

#include <functional>
#include <vector>

#include "named_func.hpp"

namespace NamedFuncUtilities {

  //Returns a vector named func that is vector_named_func filtered with filter_named_func
  NamedFunc FilterNamedFunc(NamedFunc vector_named_func, NamedFunc filter_named_func);

  //Returns a vector named func that is vector_named_func filtered with filter_named_func. Warning: inputs cannot go out of scope
  NamedFunc FilterNamedFuncCached(const NamedFunc &vector_named_func, const NamedFunc &filter_named_func);

  //Returns a named func that is input_named_func with map_function applied (entrywise) to it
  NamedFunc MapNamedFunc(NamedFunc input_named_func, std::function<double(double)> map_function);

  //Returns a named func that is input_named_funcs with map_function applied 
  //to it
  NamedFunc MultiMapNamedFunc(std::vector<NamedFunc> input_named_funcs, 
      std::function<double(std::vector<double>)> map_function);

  //Returns a named func that is input_named_funcs with map_function applied 
  //to it. Warning: inputs cannot go out of scope
  NamedFunc MultiMapNamedFuncCached(
      std::vector<const NamedFunc*> input_named_funcs, 
      const std::function<double(std::vector<double>)> &map_function);

  //Returns a scalar named func that is vector_named_func with reduce_function applied to it
  NamedFunc ReduceNamedFunc(NamedFunc vector_named_func, 
      std::function<double(std::vector<double>)> reduce_function);

  //Returns a scalar named func that is vector_named_func with reduce_function applied to it. Warning: inputs cannot go out of scope
  NamedFunc ReduceNamedFuncCached(const NamedFunc &vector_named_func, 
      const std::function<double(std::vector<double>)> &reduce_function);

  //Returns a scalar named func that is the output of reduce_function applied to the vector created by vector_named_func
  NamedFunc MultiReduceNamedFunc(std::vector<NamedFunc> vector_named_func, 
      std::function<double(std::vector<std::vector<double>>)> reduce_function);

  //Returns a scalar named func that is the output of reduce_function applied to the vector created by vector_named_func. Warning: inputs cannot go out of scope
  NamedFunc MultiReduceNamedFuncCached(
      std::vector<const NamedFunc*> vector_named_func, 
      const std::function<double(std::vector<std::vector<double>>)> 
        &reduce_function);

  //Turns a vector<NamedFunc> into one usable for a cutflow table
  std::vector<NamedFunc> progressive_cuts(std::vector<NamedFunc> vector_NamedFunc);

  //This function adds all selections and reverses the selection at reverse
  NamedFunc Nreverse1(std::vector<NamedFunc> vector_NamedFunc, unsigned int reverse);

  //Returns a NamedFunc replacing one selection (marked by skip)
  NamedFunc Nreplace1(std::vector<NamedFunc> vector_NamedFunc, NamedFunc replace, unsigned int skip);

  //Returns a NamedFunc with all but one selection (marked by skip)
  NamedFunc Nminus1(std::vector<NamedFunc> vector_NamedFunc, unsigned int skip);

  //Returns a NamedFunc without all selections in the vector skip
  NamedFunc Nminusk(std::vector<NamedFunc> vector_NamedFunc, std::vector<unsigned int> skip);

  //get the sum of a vector
  double reduce_sum(std::vector<double> data);

  //get the maximum value of a vector
  double reduce_max(std::vector<double> data);

  //get the second highest value of a vector
  double reduce_sublead(std::vector<double> data);

  //return value of second array (data[1]) corresponding to maximum of first array (data[0])
  double reduce_maxfirst(std::vector<std::vector<double>> data);

  //return value of second array (data[1]) corresponding to second highest value of first array (data[0])
  double reduce_subleadfirst(std::vector<std::vector<double>> data);
}

#endif
