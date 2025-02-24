#ifndef H_SAMPLE_LOADER
#define H_SAMPLE_LOADER

#include <cstddef>

#include <set>
#include <string>
#include <map>
#include <memory>
#include <vector>

#include "core/named_func.hpp"
#include "core/palette.hpp"
#include "core/process.hpp"

class SampleLoader{
public:
  SampleLoader();
  SampleLoader(const SampleLoader &) = default;
  SampleLoader& operator=(const SampleLoader &) = default;
  SampleLoader(SampleLoader &&) = default;
  SampleLoader& operator=(SampleLoader &&) = default;
  ~SampleLoader() = default;

  SampleLoader operator()() const;

  SampleLoader & SetMacro(const std::string &macro_key, 
                          const std::set<std::string> &values);
  SampleLoader & LoadPalette(const std::string &file, 
                             const std::string &palette);
  SampleLoader & LoadNamedFunc(const std::string &name, 
                               const NamedFunc &named_func);
  SampleLoader & ParseFile(const std::string &file_name, 
                           const std::string &config_name);
  SampleLoader & Verbose(const bool verbose);
  std::vector<std::shared_ptr<Process>> GetSamples(
      const std::set<Process::Type> types={Process::Type::data, 
      Process::Type::background, Process::Type::signal});

  std::vector<std::shared_ptr<Process>> LoadSamples(
                           const std::string &file_name, 
                           const std::string &config_name,
                           const std::set<Process::Type> types
                           ={Process::Type::data, Process::Type::background, 
                           Process::Type::signal});

  std::vector<std::shared_ptr<Process>> LoadMCSamples(
                           const std::string &file_name, 
                           const std::string &config_name);

  SampleLoader & SetProperty(const std::string &sample_name, 
                             const std::string &property,
                             const std::string &value);

private:
  class DataSample {
  public:
    DataSample(const std::string &name="");

    std::string name_;
    Process::Type type_;
    int color_;
    std::set<std::string> files_;
    NamedFunc selection_;
    std::string ntuple_type_;
  };

  std::map<std::string, std::set<std::string>> macros_;
  std::set<std::string> overwrite_macro_keys_;
  //std::map<std::string, DataSample> samples_;
  std::vector<DataSample> samples_;
  Palette palette_;
  std::map<std::string, std::shared_ptr<NamedFunc>> loaded_namedfuncs_;
  bool verbose_;

  std::set<std::string> ExpandMacros(const std::string &str_to_expand);
  bool SampleExists(const std::string sample_name);
  unsigned int SampleIndex(const std::string sample_name);


};

#endif
