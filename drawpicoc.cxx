//interface for ctypes python bindings
#include <iostream>
#include <memory>
#include <new>
#include <set>
#include <vector>

#include "TError.h"

#include "core/axis.hpp"
#include "core/hist1d.hpp"
#include "core/named_func.hpp"
#include "core/plot_maker.hpp"
#include "core/plot_opt.hpp"
#include "core/process.hpp"
#include "core/sample_loader.hpp"
#include "higgsino/hig_functions.hpp"

extern "C"
{
  //general
  void SuppressRootWarnings() {
    gErrorIgnoreLevel = 6000;
  }

  //Axis
  void* NewAxis(int nbins, double xmin, double xmax, void* var, const char* title, 
      double* cut_vals, int cut_vals_len, double* hard_cut_vals, 
      int hard_cut_vals_len) {
    std::set<double> cut_vals_set;
    std::set<double> hard_cut_vals_set;
    for (int i = 0; i < cut_vals_len; i++) {
      cut_vals_set.insert(cut_vals[i]);
    }
    for (int i = 0; i < hard_cut_vals_len; i++) {
      hard_cut_vals_set.insert(hard_cut_vals[i]);
    }
    return static_cast<void*>(new(std::nothrow) Axis(nbins, xmin, xmax, 
        *static_cast<NamedFunc*>(var), title, cut_vals_set, hard_cut_vals_set));
  }

  void DeleteAxis(void* axis) {
    delete static_cast<Axis*>(axis);
  }

  //Hist1D
  void Hist1DWeight(void* hist1d, void* weight) {
    static_cast<Hist1D*>(hist1d)->Weight(*static_cast<NamedFunc*>(weight));
  }

  void Hist1DTag(void* hist1d, const char* tag) {
    static_cast<Hist1D*>(hist1d)->Tag(tag);
  }

  void Hist1DLuminosityTag(void* hist1d, const char* tag) {
    static_cast<Hist1D*>(hist1d)->LuminosityTag(tag);
  }

  //NamedFunc
  void* NewNamedFunc(const char* function) {
    void* ret = static_cast<void*>(new(std::nothrow) NamedFunc(function));
    return ret;
  }

  void DeleteNamedFunc(void* named_func) {
    delete static_cast<NamedFunc*>(named_func);
  }

  void* CopyNamedFunc(void* named_func) {
    void* ret = static_cast<void*>(new(std::nothrow) NamedFunc(*static_cast<NamedFunc*>(named_func)));
    return ret;
  }

  void* NamedFuncAnd(void* f, void* g) {
    return static_cast<void*>(new(std::nothrow) NamedFunc(
        *(static_cast<NamedFunc*>(f)) && *(static_cast<NamedFunc*>(g))));
  }

  //PlotMaker
  void* NewPlotMaker() {
    return static_cast<void*>(new(std::nothrow) PlotMaker());
  }
  
  void DeletePlotMaker(void* plot_maker) {
    delete static_cast<PlotMaker*>(plot_maker);
  }

  void* PlotMakerPushHist1D(void* plot_maker, void* axis, void* cut, 
      void* processes, void** plot_options, int plot_options_len) {
    std::vector<PlotOpt> plot_options_vector;
    for (int i = 0; i < plot_options_len; i++) {
      plot_options_vector.push_back(*static_cast<PlotOpt*>(plot_options[i]));
    }
    return static_cast<void*>(&(static_cast<PlotMaker*>(plot_maker)->Push<Hist1D>(
        *static_cast<Axis*>(axis), *static_cast<NamedFunc*>(cut),
        *static_cast<std::vector<std::shared_ptr<Process>>*>(processes),
        plot_options_vector)));
  }

  void PlotMakerMakePlots(void* plot_maker, double luminosity, 
      const char* subdir) {
    static_cast<PlotMaker*>(plot_maker)->MakePlots(luminosity, subdir);
  }

  //PlotOpt
  void* NewPlotOpt(const char* file_name, const char* config_name) {
    return static_cast<void*>(new(std::nothrow) PlotOpt(file_name, config_name));
  }

  void DeletePlotOpt(void* plot_opt) {
    delete static_cast<PlotOpt*>(plot_opt);
  }

  //SampleLoader
  void* NewSampleLoader() {
    return static_cast<void*>(new(std::nothrow) SampleLoader());
  }

  void DeleteSampleLoader(void* sample_loader) {
    delete static_cast<SampleLoader*>(sample_loader);
  }

  void SampleLoaderLoadNamedFunc(void* sample_loader, const char* name, 
      void* named_func) {
    static_cast<SampleLoader*>(sample_loader)->LoadNamedFunc(name, 
        *static_cast<NamedFunc*>(named_func));
  }

  void* SampleLoaderLoadSamples(void* sample_loader, const char* file_name,
      const char* config_name) {
    std::vector<std::shared_ptr<Process>>* processes = new(std::nothrow) 
        std::vector<std::shared_ptr<Process>>(static_cast<SampleLoader*>(
        sample_loader)->LoadSamples(file_name, config_name));
    return static_cast<void*>(processes);
  }

  //ProcessList
  void DeleteProcessList(void* process_list) {
    delete static_cast<std::vector<std::shared_ptr<Process>>*>(process_list);
  }

  //Higfuncs
  void* HigFuncsMetTrigger() {
    return static_cast<void*>(new(std::nothrow) NamedFunc(Higfuncs::met_trigger));
  }

  void* HigFuncsFinalPassFilters() {
    return static_cast<void*>(new(std::nothrow) NamedFunc(Higfuncs::final_pass_filters));
  }

  void* HigFuncsFinalWeight() {
    return static_cast<void*>(new(std::nothrow) NamedFunc(Higfuncs::final_weight));
  }

}
