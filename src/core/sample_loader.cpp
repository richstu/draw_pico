/*! \class SampleLoader
  \brief A class to load sample (Process) information from a text file

*/

/*! \class SampleLoader::DataSample
  \brief Small class containing properties to be assigned to a Process

*/

#include "core/sample_loader.hpp"

#include <cmath>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "core/named_func.hpp"
#include "core/palette.hpp"
#include "core/process.hpp"
#include "core/utilities.hpp"

using namespace std;

/*!\brief Standard constructor

    \param[in] name Name for this DataSample/Process
  */
SampleLoader::DataSample::DataSample(const string &name):
  name_(name),
  type_(Process::Type::background),
  color_(1),
  files_(set<string>()),
  selection_("1"),
  ntuple_type_("Baby_pico"){
}

/*!\brief Standard constructor
  */
SampleLoader::SampleLoader():
  macros_(map<string, set<string>>()),
  overwrite_macro_keys_(set<string>()),
  samples_(vector<DataSample>()),
  palette_(Palette()),
  loaded_namedfuncs_(map<string,shared_ptr<NamedFunc>>()),
  verbose_(false){
}

SampleLoader SampleLoader::operator()() const{
  return *this;
}

/*!\brief Allows caller to define macros used in sample paths

    \param[in] macro_key Name of the macro
    \param[in] values What the macro expands to when used 
  */
SampleLoader & SampleLoader::SetMacro(const string &macro_key, 
                                      const set<string> &values) {
  overwrite_macro_keys_.insert(macro_key);
  if (macros_.count(macro_key)==0) {
    macros_.insert(make_pair(macro_key, values));
  }
  else {
    macros_[macro_key] = values;
  }
  return *this;
}

/*!\brief Load color Palette from file

    \param[in] file Name of file to load from
    \param[in] palette Name of palette within file
  */
SampleLoader & SampleLoader::LoadPalette(const string &file, 
                                         const string &palette) {
  palette_ = Palette(file, palette);
  return *this;
}

/*!\brief Set verbosity

    \param[in] verbose Verbose enable
  */
SampleLoader & SampleLoader::Verbose(const bool verbose) {
  verbose_ = verbose;
  return *this;
}

/*!\brief Makes a NamedFunc accessible for sample selections

    \param[in] name Name of the NamedFunc
    \param[in] named_func NamedFunc
  */
SampleLoader & SampleLoader::LoadNamedFunc(const string &name, 
                                           const NamedFunc &named_func) {
  if (loaded_namedfuncs_.count(name)==0) {
    loaded_namedfuncs_.insert(make_pair(name, make_shared<NamedFunc>(named_func)));
  }
  else {
    loaded_namedfuncs_[name] = make_shared<NamedFunc>(named_func);
  }
  return *this;
}

/*!\brief Loads samples from file and returns vector of Processes

    \param[in] file_name Name of file to load from
    \param[in] config_name Name of configuration in file
  */
vector<shared_ptr<Process>> SampleLoader::LoadSamples(const string &file_name, 
    const string &config_name){
  ParseFile(file_name, config_name);
  return GetSamples();
}

/*!\brief Retreives sample information from file

    \param[in] file_name Name of file to load from
    \param[in] config_name Name of configuration in file
  */
SampleLoader & SampleLoader::ParseFile(const string &file_name, 
    const string &config_name){
  ifstream file(file_name);
  string line;
  string current_config = "";
  string current_sample = "";
  int line_num = 0;
  while(getline(file, line)){
    ++line_num;
    if (line.substr(0,4) == "#add") { //macro line
      size_t space_one = line.find(' ');
      size_t space_two = line.find(' ', space_one+1);
      if (space_one==string::npos || space_two==string::npos) {
        ERROR("Could not find spaces in macro from line "+to_string(line_num));
      }
      string macro_key = line.substr(space_one+1, space_two-space_one-1);
      string set_value = line.substr(space_two+1, line.size()-space_two-1);
      if (overwrite_macro_keys_.count(macro_key)==0) {
        //macro not overwritten by SetMacro method; add value
        if (verbose_) 
          cout << "Adding macro key " << macro_key << ", value " << set_value << endl;
        if (macros_.count(macro_key)==0) {
          macros_.insert(make_pair(macro_key, set<string>({set_value})));
        }
        else {
          macros_[macro_key].insert(set_value);
        }
      }
    }
    ReplaceAll(line, " ", "");
    ReplaceAll(line, "\t", "");
    size_t start  = line.find("[[");
    size_t end = line.find("]]");
    if(start<end && start != string::npos && end != string::npos){
      current_config = line.substr(start+2, end-(start+1)-1);
    }else if(current_config == config_name
             && line.size()
             && line.at(0)!='#'){
      start  = line.find("[");
      end = line.find("]");
      if(start<end && start != string::npos && end != string::npos){
        current_sample = line.substr(start+1, end-start-1);
        if (!SampleExists(current_sample)) {
          samples_.push_back(DataSample(current_sample));
        }
      } else if (current_sample != "") {
        size_t pos = line.find("=");
        if(pos == string::npos) continue;
        string prop_name = line.substr(0,pos);
        string value = line.substr(pos+1);
        SetProperty(current_sample, prop_name, value);
      }
    }
  }
  return *this;
}

/*!\brief Returns vector of Processes for currently loaded samples

  */
vector<shared_ptr<Process>> SampleLoader::GetSamples() {
  vector<shared_ptr<Process>> processes;
  for (DataSample this_sample : samples_) {
    //currently only picos supported; easy enough to add more Baby children, though
    if (this_sample.ntuple_type_ == "Baby_pico") {
      if (verbose_) {
        std::cout << this_sample.name_ << "\n";
      }
      processes.push_back(Process::MakeShared<Baby_pico>(this_sample.name_,
          this_sample.type_, this_sample.color_, this_sample.files_, 
          this_sample.selection_));
      processes.back()->SetLineColor(this_sample.color_);
    }
  }
  return processes;
}

/*!\brief Expands string by replacing macros with all possible values

    \param[in] str_to_expand String for which to replace macros
  */
set<string> SampleLoader::ExpandMacros(const string & str_to_expand) {
  vector<set<string>> expanded_string = {{str_to_expand}};
  for (pair<string, set<string>> macro : macros_) {
    string macro_expand = "$"+macro.first;
    unsigned i = 0;
    while (i < expanded_string.size()) {
      if (expanded_string[i].size() == 1) {
        string unexpanded_substring;
        for (string expanded_substring_set_string : expanded_string[i]) {
          unexpanded_substring = expanded_substring_set_string;
        }
        size_t macro_location = unexpanded_substring.find(macro_expand);
        if (macro_location != string::npos) {
          string prefix = unexpanded_substring.substr(0,macro_location);
          string suffix = unexpanded_substring.substr(
              macro_location+macro_expand.size());
          expanded_string.insert(expanded_string.begin()+i, {prefix});
          expanded_string.insert(expanded_string.begin()+i+1, macro.second);
          expanded_string[i+2] = {suffix};
          i += 3;
        }
        else i++;
      }
      else i++; //no support for nested macros
    } //loop over expanded_string
  }
  set<string> string_combinations = {""};
  for (set<string> combination_set : expanded_string) {
    set<string> next_string_combinations;
    for (string prefix : string_combinations) {
      for (string suffix : combination_set) {
        next_string_combinations.insert(prefix+suffix);
      }
    }
    string_combinations = next_string_combinations;
  }
  //debug printout
  if (verbose_) {
    for (string comb : string_combinations) {
      cout << comb << endl;
    }
  }
  return string_combinations;
}

/*!\brief Sets property of a given sample to a given value

    \param[in] sample_name Name of sample to modify
    \param[in] property Name of property to modify
    \param[in] value What to set said property to
  */
SampleLoader & SampleLoader::SetProperty(const string &sample_name,
                               const string &property,
                               const string &value){
  if (!SampleExists(sample_name)) {
    samples_.push_back(DataSample(sample_name));
    cout << "WARNING: implicitly creating new sample to set property" << endl;
  }
  unsigned int sample_index = SampleIndex(sample_name);
  if(property == "ProcessType"){
    samples_[sample_index].type_ = static_cast<Process::Type>(stoi(value));
  }else if(property == "PaletteColor"){
    samples_[sample_index].color_ = palette_(value);
  }else if(property == "Color"){
    samples_[sample_index].color_ = stoi(value);
  }else if(property == "Files"){
    samples_[sample_index].files_ = ExpandMacros(value);
  }else if(property == "Selection"){
    if (loaded_namedfuncs_.count(value) != 0) {
      samples_[sample_index].selection_ = *loaded_namedfuncs_[value];
    } else {
      samples_[sample_index].selection_ = value;
    }
  }else{
    DBG("Did not understand property name "<<property);
  }
  return *this;
}

/*!\brief Returns true if there already exists a sample with sample_name

    \param[in] sample_name name of sample to check
  */
bool SampleLoader::SampleExists(const string sample_name) {
  for (DataSample sample : samples_) {
    if (sample.name_ == sample_name)
      return true;
  }
  return false;
}


/*!\brief Returns index of sample if there already exists a sample with
          sample_name and 999 otherwise

    \param[in] sample_name name of sample to check
  */
unsigned int SampleLoader::SampleIndex(const string sample_name) {
  for (unsigned int isample = 0; isample<samples_.size(); isample++) {
    if (samples_[isample].name_ == sample_name) {
      return isample;
    }
  }
  return 999;
}
