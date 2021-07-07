#ifndef H_PLOT_MAKER
#define H_PLOT_MAKER

#include <functional>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "ROOT/RResultPtr.hxx"

#include "core/figure.hpp"
#include "core/plot_opt.hpp"

using ProcessedDataFrame = ROOT::RDF::RInterface<ROOT::Detail::RDF::RJittedFilter, void>;

class Process;

class PlotMaker{
public:
  PlotMaker();
  PlotMaker(const PlotMaker &) = default;
  PlotMaker& operator=(const PlotMaker &) = default;
  PlotMaker(PlotMaker &&) = default;
  PlotMaker& operator=(PlotMaker &&) = default;
  ~PlotMaker() = default;

  template<typename FigureType, typename... Args>
  FigureType & Push(Args&&... args){
    figures_.emplace_back(static_cast<Figure*>(new FigureType(args...)));
    return *static_cast<FigureType*>(figures_.back().get());
  }

  void MakePlots(double luminosity,
                 const std::string &subdir = "");
  void MakePlotsRdf(double luminosity,
                    const std::string &subdir = "");

  const std::unique_ptr<Figure> & GetFigure(std::string tag) const;
  const std::vector<std::unique_ptr<Figure> > & Figures() const;
  template<typename FigureType>
  FigureType * GetLast(){
    FigureType *out = nullptr;
    for(auto f = figures_.crbegin(); f != figures_.crend(); ++f){
      if((out = dynamic_cast<FigureType*>(f->get()))) return out;
    }
    return nullptr;
  }
  void Clear();

  void SetEventVetoData(void * eventVetoData);
  void DefineRdfColumns(
      std::function<ProcessedDataFrame(ProcessedDataFrame&)> rdf_column_func);
  void DefineRdfColumnsData(
      std::function<ProcessedDataFrame(ProcessedDataFrame&)> rdf_column_func);
  void DefineRdfColumnsMC(
      std::function<ProcessedDataFrame(ProcessedDataFrame&)> rdf_column_func);

  std::string tree_name_;
  bool multithreaded_;
  bool min_print_;
  bool print_2d_figures_;
  long max_entries_;
  void * event_veto_data_;

private:
  std::vector<std::unique_ptr<Figure> > figures_;//!<Figures to be produced
  std::vector<std::function<ProcessedDataFrame(ProcessedDataFrame&)>> rdf_column_funcs_mc_;//!<Functions to be applied to define columns in RDataFrame for MC
  std::vector<std::function<ProcessedDataFrame(ProcessedDataFrame&)>> rdf_column_funcs_data_;//!<Functions to be applied to define columns in RDataFrame for data

  void GetYields();
  long GetYield(Baby *baby_ptr);

  std::set<Baby*> GetBabies() const;
  std::set<const Process *> GetProcesses() const;
  std::set<Figure::FigureComponent*> GetComponents(const Process *process) const;
};

#endif
