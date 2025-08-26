#ifndef H_MVA_WRAPPER
#define H_MVA_WRAPPER

#include <string>
#include <vector>
#include <unordered_map>

#include "core/baby.hpp"
#include "core/named_func.hpp"

#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"
#include "Rtypes.h"

class MVAWrapper{
public:
  MVAWrapper(std::string name);
  MVAWrapper(std::string name, std::vector<std::string> input_names, 
             std::vector<NamedFunc> inputs, 
             std::vector<std::string> variation_names, 
             std::vector<std::vector<NamedFunc>> variation_inputs,
             std::string weights_filename);
  MVAWrapper(const MVAWrapper &) = default;
  MVAWrapper & operator=(const MVAWrapper &) = default;
  MVAWrapper(MVAWrapper &&) = default;
  MVAWrapper & operator=(MVAWrapper &&) = default;
  ~MVAWrapper() = default;

  //MVAWrapper SetWeightsFile(std::string filename);
  MVAWrapper & SetVariable(std::string name, const NamedFunc &variable);
  MVAWrapper & SetAltVariable(std::string variation, 
                              const NamedFunc &variable);
  MVAWrapper & SetSpectator(std::string name, const NamedFunc &variable);
  MVAWrapper & BookMVA(std::string weights_filename);
  MVAWrapper & SetVariableRef(std::string name, const NamedFunc &variable);
  MVAWrapper & SetAltVariableRef(std::string variation, 
                                 const NamedFunc &variable);

  double GetDiscriminantScore(const Baby &b, std::string variation="");
  NamedFunc GetDiscriminant(std::string variation="");

private:
  MVAWrapper() = delete;

  std::string name_;
  std::vector<NamedFunc> variables_;
  std::vector<std::string> variable_names_;
  std::unordered_map<std::string, std::vector<NamedFunc>> alt_variables_;
  std::vector<Float_t> variable_values_;
  std::vector<NamedFunc> spectators_;
  std::vector<std::string>  spectator_names_;
  std::vector<float> spectator_values_;
  TMVA::Reader mva_reader_;
  bool booked_;
  bool use_references_;
  std::vector<const NamedFunc*> variables_ref_;
  std::unordered_map<std::string, std::vector<const NamedFunc*>> 
      alt_variables_ref_;

};

#endif
