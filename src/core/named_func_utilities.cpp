/**
 * Utilities for a functional programming style using NamedFuncs
 */
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "core/baby.hpp"
#include "core/named_func.hpp"
#include "core/named_func_utilities.hpp"

namespace NamedFuncUtilities {

  //Returns a vector named func that is vector_named_func filtered with filter_named_func
  NamedFunc FilterNamedFunc(NamedFunc vector_named_func, NamedFunc filter_named_func) {
    return NamedFunc("FilterNamedFunc("+vector_named_func.Name()+","+filter_named_func.Name()+")",[vector_named_func,filter_named_func](const Baby &b) -> NamedFunc::VectorType{
      std::vector<double> mapped_named_func;
      std::vector<double> vector_named_func_ = vector_named_func.GetVector(b);
      std::vector<double> filter_named_func_ = filter_named_func.GetVector(b);
      for (unsigned i = 0; i < vector_named_func_.size(); i++) {
        if (filter_named_func_[i]) {
          mapped_named_func.push_back(vector_named_func_[i]);
        }
      }
      return mapped_named_func;
    });
  }
  
  //Returns a named func that is input_named_func with map_function applied (entrywise) to it
  NamedFunc MapNamedFunc(NamedFunc input_named_func, std::function<double(double)> map_function) {
    if (input_named_func.IsScalar()) {
      return NamedFunc("MapNamedFunc("+input_named_func.Name()+")",[input_named_func,map_function](const Baby &b) -> NamedFunc::ScalarType{
        return map_function(input_named_func.GetScalar(b));
      });
    }
    return NamedFunc("MapNamedFunc("+input_named_func.Name()+")",[input_named_func,map_function](const Baby &b) -> NamedFunc::VectorType{
      std::vector<double> mapped_named_func;
      std::vector<double> input_named_func_ = input_named_func.GetVector(b);
      for (unsigned i = 0; i < input_named_func_.size(); i++) {
        mapped_named_func.push_back(map_function(input_named_func_[i]));
      }
      return mapped_named_func;
    });
  }
  
  //Returns a scalar named func that is vector_named_func with reduce_function applied to it
  NamedFunc ReduceNamedFunc(NamedFunc vector_named_func, std::function<double(std::vector<double>)> reduce_function) {
    return NamedFunc("ReduceNamedFunc("+vector_named_func.Name()+")",[vector_named_func,reduce_function](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<double> vector_named_func_ = vector_named_func.GetVector(b);
      return reduce_function(vector_named_func_);
    });
  }
  
  //Returns a scalar named func that is the output of reduce_function applied to the vector created by vector_named_func
  NamedFunc MultiReduceNamedFunc(std::vector<NamedFunc> vector_named_func, std::function<double(std::vector<std::vector<double>>)> reduce_function) {
    std::string named_func_name = "MultiReduceNamedFunc(";
    bool first = true;
    for (NamedFunc named_func : vector_named_func) {
      if (!first) named_func_name += ",";
      named_func_name += named_func.Name();
      first = false;
    }
    named_func_name += ")";
    return NamedFunc(named_func_name,[vector_named_func,reduce_function](const Baby &b) -> NamedFunc::ScalarType{
      std::vector<std::vector<double>> vector_named_func_;
      for (NamedFunc named_func : vector_named_func) {
        vector_named_func_.push_back(named_func.GetVector(b));
      }
      return reduce_function(vector_named_func_);
    });
  }

  //get the sum of a vector
  double reduce_sum(std::vector<double> data) {
    double sum = 0;
    for (double item : data)
      sum += item;
    return sum;
  }

  //get the maximum value of a vector
  double reduce_max(std::vector<double> data) {
    if (data.size()==0) return 0;
    double max = data[0];
    for (double item : data)
      if (item > max)
        max = item;
    return max;
  }

  //get the second highest value of a vector
  double reduce_sublead(std::vector<double> data) {
    if (data.size()<2) return 0;
    double max = 0, next = 0;
    if (data[0] > data[1]) {
      max = data[0];
      next = data[1];
    }
    else {
      max = data[1];
      next = data[0];
    }
    for (unsigned i = 2; i < data.size(); i++) {
      if (data[i] > max) {
        next = max;
        max = data[i];
      }
      else if (data[i] > next) {
        next = data[i];
      }
    }
    return next;
  }

  //return value of second array (data[1]) corresponding to maximum of first array (data[0])
  double reduce_maxfirst(std::vector<std::vector<double>> data) {
    if (data.size()<2) return 0;
    if (data[0].size() != data[1].size()) return 0;
    if (data[0].size()==0) return 0;
    double max = data[0][0];
    double ret = data[1][0];
    for (unsigned i = 1; i < data[0].size(); i++) {
      if (data[0][i] > max) {
        max = data[0][i];
        ret = data[1][i];
      }
    }
    return ret;
  }
  
  //return value of second array (data[1]) corresponding to second highest value of first array (data[0])
  double reduce_subleadfirst(std::vector<std::vector<double>> data) {
    if (data.size()<2) return 0;
    if (data[0].size() != data[1].size()) return 0;
    if (data[0].size()<2) return 0;
    double max = 0, next = 0, max_ret = 0, next_ret = 0;
    if (data[0][0] > data[0][1]) {
      max = data[0][0];
      next = data[0][1];
      max_ret = data[1][0];
      next_ret = data[1][1];
    }
    else {
      max = data[0][1];
      next = data[0][0];
      max_ret = data[1][1];
      next_ret = data[1][0];
    }
    for (unsigned i = 2; i < data[0].size(); i++) {
      if (data[0][i] > max) {
        next = max;
        max = data[0][i];
        next_ret = max_ret;
        max_ret = data[1][i];
      }
      else if (data[0][i] > next) {
        next = data[0][i];
        next_ret = data[1][i];
      }
    }
    return next_ret;
  }

}
