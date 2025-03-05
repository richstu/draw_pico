#ifndef H_MVA_WRAPPER
#define H_MVA_WRAPPER

#include <string>
#include <vector>

#include "core/named_func.hpp"

#include "TMVA/Reader.h"
#include "TMVA/Configurable.h"
#include "Rtypes.h"

class MVAWrapper{
public:
  MVAWrapper(std::string name);
  MVAWrapper(const MVAWrapper &) = default;
  MVAWrapper & operator=(const MVAWrapper &) = default;
  MVAWrapper(MVAWrapper &&) = default;
  MVAWrapper & operator=(MVAWrapper &&) = default;
  ~MVAWrapper() = default;

  //MVAWrapper SetWeightsFile(std::string filename);
  MVAWrapper & SetVariable(std::string name, NamedFunc variable);
  MVAWrapper & SetSpectator(std::string name, NamedFunc variable);
  MVAWrapper & BookMVA(std::string weights_filename);

  NamedFunc GetDiscriminant();

private:
  MVAWrapper() = delete;

  std::string name_;
  std::vector<NamedFunc> variables_;
  std::vector<std::string> variable_names_;
  std::vector<Float_t> variable_values_;
  std::vector<NamedFunc> spectators_;
  std::vector<std::string>  spectator_names_;
  std::vector<float> spectator_values_;
  TMVA::Reader mva_reader_;
  bool booked_;
  Long64_t previous_event_;
  float cached_value_;

};

#endif
