#ifndef H_NAMED_FUNC_UTILITIES
#define H_NAMED_FUNC_UTILITIES

#include <functional>
#include <vector>

#include "named_func.hpp"

namespace NamedFuncUtilities {

  //Returns a vector named func that is vector_named_func filtered with filter_named_func
  NamedFunc FilterNamedFunc(NamedFunc vector_named_func, NamedFunc filter_named_func);

  //Returns a named func that is input_named_func with map_function applied (entrywise) to it
  NamedFunc MapNamedFunc(NamedFunc input_named_func, std::function<double(double)> map_function);

  //Returns a scalar named func that is vector_named_func with reduce_function applied to it
  NamedFunc ReduceNamedFunc(NamedFunc vector_named_func, 
      std::function<double(std::vector<double>)> reduce_function);

  //Returns a scalar named func that is the output of reduce_function applied to the vector created by vector_named_func
  NamedFunc MultiReduceNamedFunc(std::vector<NamedFunc> vector_named_func, 
      std::function<double(std::vector<std::vector<double>>)> reduce_function);

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
