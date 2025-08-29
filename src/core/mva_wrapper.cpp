/*! \class MVAWrapper
 * \brief A class to wrap TMVA::Reader for draw_pico usage
*/

#include "core/mva_wrapper.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/named_func.hpp"

#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"
#include "Rtypes.h"

/*!\brief Standard constructor
*/
MVAWrapper::MVAWrapper(std::string name) :
  name_(name),
  variables_(std::vector<NamedFunc>()),
  variable_names_(std::vector<std::string>()),
  alt_variables_(std::unordered_map<std::string,std::vector<NamedFunc>>()),
  variable_values_(std::vector<float>()),
  spectators_(std::vector<NamedFunc>()),
  spectator_names_(std::vector<std::string>()),
  spectator_values_(std::vector<float>()),
  mva_reader_(TMVA::Reader()),
  booked_(false)
{}

/*!\brief Assigns a variable for the MVA and associates a NamedFunc
 * Currently, this class uses pointers to elements of a vector, which is kludge-y but will
 * work as long as the vector length is not modified
 * \param[in] name - variable name in MVA
 * \param[in] varaible - NamedFunc that will be evaluated to create this variable
*/
MVAWrapper & MVAWrapper::SetVariable(std::string name, NamedFunc variable) {
  if (!booked_) {
    variables_.push_back(variable);
    variable_names_.push_back(name);
    variable_values_.push_back(0);
  }
  else {
    throw std::runtime_error("Cannot add variables after booking.");
  }
  return *this;
}

/*!\brief Assigns a variant version of a variable for the MVA and associates a 
 * NamedFunc
 * \param[in] variation - variation name
 * \param[in] varaible - NamedFunc that will be evaluated to create variable
*/
MVAWrapper & MVAWrapper::SetAltVariable(std::string variation, 
                                        NamedFunc variable) {
  if (!booked_) {
    if (alt_variables_.count(variation)==0)
      alt_variables_[variation] = std::vector<NamedFunc>();
    alt_variables_[variation].push_back(variable);
  }
  else {
    throw std::runtime_error("Cannot add variables after booking.");
  }
  return *this;
}

MVAWrapper & MVAWrapper::SetSpectator(std::string name, NamedFunc variable) {
  if (!booked_) {
    spectators_.push_back(variable);
    spectator_names_.push_back(name);
    spectator_values_.push_back(0);
  }
  else {
    throw std::runtime_error("Cannot add spectator after booking.");
  }
  return *this;
}

/*!\brief Books an MVA and makes it ready for evaluation
 * \param[in] weights_filename - filename for the MVA weights files
*/
MVAWrapper & MVAWrapper::BookMVA(std::string weights_filename) {
  std::cout << "Creating an MVA reader. For the reader to work, ensure "
               "PlotMaker is set to single-threaded." << std::endl;
  for (auto& alt_var : alt_variables_) {
    if (alt_var.second.size() != variables_.size())
      throw std::runtime_error("Alternate variable set size differs.");
  }
  for (unsigned i = 0; i < variables_.size(); i++) {
    mva_reader_.AddVariable(variable_names_[i], &(variable_values_[i]));
  }
  for (unsigned i = 0; i < spectators_.size(); i++) {
    mva_reader_.AddSpectator(spectator_names_[i], &(spectator_values_[i]));
  }

  //MVA type currently hardcoded to BDT
  mva_reader_.BookMVA("BDT",weights_filename.c_str());
  booked_ = true;
  return *this;
}

/*!\brief Returns a NamedFunc that yields discriminant for an event
*/
NamedFunc MVAWrapper::GetDiscriminant(std::string variation) {
  return NamedFunc((name_+"Score").c_str(),[&](const Baby &b) 
      -> NamedFunc::ScalarType{ 
    for (unsigned i = 0; i < variables_.size(); i++) {
      if (variation=="") {
        variable_values_[i] = variables_[i].GetScalar(b);
      }
      else {
        variable_values_[i] = alt_variables_[variation][i].GetScalar(b);
      }
    }
    return mva_reader_.EvaluateMVA("BDT");
  }).EnableCaching(true);
}
