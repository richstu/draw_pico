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
  booked_(false),
  use_references_(false),
  variables_ref_(std::vector<const NamedFunc*>()),
  alt_variables_ref_(
      std::unordered_map<std::string,std::vector<const NamedFunc*>>())
{}

/*!\brief All-in-one constructor for thread_local use
*/
MVAWrapper::MVAWrapper(std::string name, std::vector<std::string> input_names, 
                       std::vector<NamedFunc> inputs, 
                       std::vector<std::string> variation_names, 
                       std::vector<std::vector<NamedFunc>> variation_inputs,
                       std::string weights_filename) :
    name_(name),
    variables_(std::vector<NamedFunc>()),
    variable_names_(std::vector<std::string>()),
    alt_variables_(std::unordered_map<std::string,std::vector<NamedFunc>>()),
    variable_values_(std::vector<float>()),
    spectators_(std::vector<NamedFunc>()),
    spectator_names_(std::vector<std::string>()),
    spectator_values_(std::vector<float>()),
    mva_reader_(TMVA::Reader()),
    booked_(false),
    use_references_(false),
    variables_ref_(std::vector<const NamedFunc*>()),
    alt_variables_ref_(
        std::unordered_map<std::string,std::vector<const NamedFunc*>>()) {
  if (input_names.size() != inputs.size())
    throw std::runtime_error("MVA: inputs and name sizes don't match.");
  if (variation_names.size() != variation_inputs.size())
    throw std::runtime_error("MVA: inputs and name sizes don't match.");
  for (unsigned ivar = 0; ivar < inputs.size(); ivar++) {
    this->SetVariable(input_names[ivar], inputs[ivar]);
  }
  for (unsigned ivar = 0; ivar < variation_names.size(); ivar++) {
    for (unsigned iin = 0; iin < variation_inputs[ivar].size(); iin++) {
      this->SetAltVariable(variation_names[ivar], variation_inputs[ivar][iin]);
    }
  }
  this->BookMVA(weights_filename);
}

/*!\brief Assigns a variable for the MVA and associates a NamedFunc
 * Currently, this class uses pointers to elements of a vector, which is kludge-y but will
 * work as long as the vector length is not modified
 * \param[in] name - variable name in MVA
 * \param[in] varaible - NamedFunc that will be evaluated to create this variable
*/
MVAWrapper & MVAWrapper::SetVariable(std::string name, 
                                     const NamedFunc &variable) {
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
                                        const NamedFunc &variable) {
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

/*!\brief Assigns a variable for the MVA and associates a NamedFunc pointer
 * \param[in] name - variable name in MVA
 * \param[in] varaible - NamedFunc that will be evaluated to create this variable
*/
MVAWrapper & MVAWrapper::SetVariableRef(std::string name, 
                                        const NamedFunc &variable) {
  if (!booked_) {
    use_references_ = true;
    variables_ref_.push_back(&variable);
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
MVAWrapper & MVAWrapper::SetAltVariableRef(std::string variation, 
                                           const NamedFunc &variable) {
  if (!booked_) {
    use_references_ = true;
    if (alt_variables_ref_.count(variation)==0)
      alt_variables_ref_[variation] = std::vector<const NamedFunc*>();
    alt_variables_ref_[variation].push_back(&variable);
  }
  else {
    throw std::runtime_error("Cannot add variables after booking.");
  }
  return *this;
}

MVAWrapper & MVAWrapper::SetSpectator(std::string name, 
                                      const NamedFunc &variable) {
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
    if (alt_var.second.size() != variables_.size()) {
      std::cout << "Variation " << alt_var.first << std::endl;
      throw std::runtime_error("Alternate variable set size differs.");
    }
  }
  for (auto& alt_var : alt_variables_ref_) {
    if (alt_var.second.size() != variables_ref_.size()) {
      std::cout << "Variation " << alt_var.first << std::endl;
      throw std::runtime_error("Alternate variable set size differs.");
    }
  }
  for (unsigned i = 0; i < variable_names_.size(); i++) {
    mva_reader_.AddVariable(variable_names_[i], &(variable_values_[i]));
  }
  for (unsigned i = 0; i < spectator_names_.size(); i++) {
    mva_reader_.AddSpectator(spectator_names_[i], &(spectator_values_[i]));
  }

  //MVA type currently hardcoded to BDT
  mva_reader_.BookMVA("BDT",weights_filename.c_str());
  booked_ = true;
  return *this;
}

/*!\brief Gets discriminant for an event
*/
double MVAWrapper::GetDiscriminantScore(const Baby &b, std::string variation) {
  for (unsigned i = 0; i < variable_values_.size(); i++) {
    if (variation=="") {
      if (use_references_) {
        variable_values_[i] = variables_ref_[i]->GetScalar(b);
      }
      else {
        variable_values_[i] = variables_[i].GetScalar(b);
      }
    }
    else {
      if (use_references_) {
        variable_values_[i] = alt_variables_ref_[variation][i]->GetScalar(b);
      }
      else {
        variable_values_[i] = alt_variables_[variation][i].GetScalar(b);
      }
    }
  }
  return mva_reader_.EvaluateMVA("BDT");
}

/*!\brief Returns a NamedFunc that yields discriminant for an event
*/
NamedFunc MVAWrapper::GetDiscriminant(std::string variation) {
  return NamedFunc((name_+"Score").c_str(),[variation, this](const Baby &b) 
      -> NamedFunc::ScalarType{ 
    return GetDiscriminantScore(b, variation);
  }).EnableCaching(true);
}
