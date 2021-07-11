#ifndef H_EFFICIENCY_PLOT
#define H_EFFICIENCY_PLOT

#include <vector>
#include <utility>
#include <memory>
#include <set>
#include <limits>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TH1D.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"

#include "core/figure.hpp"
#include "core/process.hpp"
#include "core/axis.hpp"
#include "core/plot_opt.hpp"

class EfficiencyPlot final: public Figure{
public:
  class SingleEfficiencyPlot final: public Figure::FigureComponent{
  public:
    SingleEfficiencyPlot(const EfficiencyPlot &stack,
                 const std::shared_ptr<Process> &process,
                 const TH1D &denominator_hist, 
                 const TH1D &numerator_hist);
    ~SingleEfficiencyPlot() = default;

    TH1D raw_denominator_hist_;//!<Histogram storing distribution before stacking and luminosity weighting
    TH1D raw_numerator_hist_;//!<Histogram storing distribution before stacking and luminosity weighting

    void RecordEvent(const Baby &baby) final;
    void BookResult(ROOT::RDF::RNode data_frame, int &rdf_plot_idx) final;

    void GetResult() final;

  private:
    SingleEfficiencyPlot() = delete;
    SingleEfficiencyPlot(const SingleEfficiencyPlot &) = delete;
    SingleEfficiencyPlot& operator=(const SingleEfficiencyPlot &) = delete;
    SingleEfficiencyPlot(SingleEfficiencyPlot &&) = delete;
    SingleEfficiencyPlot& operator=(SingleEfficiencyPlot &&) = delete;

    NamedFunc proc_and_hist_cut_, numerator_cut_;
    NamedFunc::VectorType cut_vector_, wgt_vector_, val_vector_, numerator_cut_vector_;
    ROOT::RDF::RResultPtr<TH1D> booked_raw_numerator_hist_ptr_;
    ROOT::RDF::RResultPtr<TH1D> booked_raw_denominator_hist_ptr_;
  };

  EfficiencyPlot(const Axis &xaxis, const NamedFunc &denominator_cut, const NamedFunc &numerator_cut,
                 const std::vector<std::shared_ptr<Process> > &processes,
                 const bool draw_histograms=true, const std::vector<PlotOpt> &plot_options = {PlotOpt()});
  EfficiencyPlot(EfficiencyPlot &&) = default;
  EfficiencyPlot& operator=(EfficiencyPlot &&) = default;
  ~EfficiencyPlot() = default;

  //functions overwriting the virtual functions from Figure
  void Print(double luminosity,
             const std::string &subdir) final;
  std::string GetTag() const final;
  std::set<const Process*> GetProcesses() const final;
  FigureComponent * GetComponent(const Process *process) final;

  std::string Name() const;
  std::string Title() const;
  std::string YAxisText() const; 

  //functions for specifying additional paramters
  EfficiencyPlot & Weight(const NamedFunc &weight);
  EfficiencyPlot & Tag(const std::string &tag);
  EfficiencyPlot & LuminosityTag(const std::string &tag);
  EfficiencyPlot & FixTitle(const std::string &title);
  EfficiencyPlot & YTitle(const std::string &ytitle);

  Axis xaxis_;//!<Specification of content: plotted variable, binning, etc.
  NamedFunc cut_;//!<Event selection
  NamedFunc numerator_cut_;//!<Event selection
  std::vector<PlotOpt> plot_options_;//!<Styles with which to draw plot
  bool draw_histograms_;//!<flag for drawing histograms in addition to efficiency plot

  NamedFunc weight_;//!<Event weight
  std::string tag_;//!<Filename tag to identify plot
  std::string ytitle_;//!<Title used on y-axis of plot
  std::string luminosity_tag_;//!<Filename tag to identify plot
  std::string fixed_title_;//!<title fixed using FixTitle method, empty string indicates default

  std::unique_ptr<TGraphAsymmErrors> background_ratio_plot_;
  std::vector<std::unique_ptr<TGraphAsymmErrors>> signal_ratio_plots_;
  std::vector<std::unique_ptr<TGraphAsymmErrors>> data_ratio_plots_;

private:
  std::vector<std::unique_ptr<SingleEfficiencyPlot> > backgrounds_;//!<Background components of the figure
  std::vector<std::unique_ptr<SingleEfficiencyPlot> > signals_;//!<Signal components of the figure
  std::vector<std::unique_ptr<SingleEfficiencyPlot> > datas_;//!<Data components of the figure

  Color_t background_color_;//!<color used for background plots
  std::vector<std::string> signal_names_;//!<names of the signal samples, for legends
  std::vector<std::string> data_names_;//!<names of the data samples, for legends
  mutable double luminosity_;//!<Luminosity currently being drawn

  EfficiencyPlot(const EfficiencyPlot &) = delete;
  EfficiencyPlot& operator=(const EfficiencyPlot &) = delete;
  EfficiencyPlot() = delete;

  void SetRatioPlotDrawOptions(std::unique_ptr<TGraphAsymmErrors> &ratio_plot, 
                               double y_min, double y_max);
  void SetLinearPlotDrawOptions(std::unique_ptr<TH1D> &linear_plot, 
                                double y_max, bool is_denominator);
  //std::vector<TLine> GetCutLines(double y_min, double y_max, bool adjust_bottom) const;
  //TODO: implement

  const std::vector<std::unique_ptr<SingleEfficiencyPlot> >& GetComponentList(const Process *process);

};

#endif //H_EFFICIENCY_PLOT
