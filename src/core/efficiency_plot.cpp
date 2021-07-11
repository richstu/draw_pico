/*! \class EfficiencyPlot 
 
  \brief A set of (overlayed) efficiency plots

  EfficiencyPlot contains the information necessary to produce an efficiency
  plot, possibly featuring the denominator and numerator histograms or multiple
  processes. The
  content and style of the plot are maintained (mostly) independently, so that
  once the histograms have been filled with data, redrawing with multiple styles
  has minimal overhead.

  To produce a plot, the component histograms in EfficiencyPlot::background_denominators_,
  EfficiencyPlot::background_numerators_, EfficiencyPlot::signal_denominators_,
  EfficiencyPlot::signal_numerators_, EfficiencyPlot::data_denominators_, and
  EfficiencyPlot::data_numerators_ must be filled (e.g., by
  PlotMaker). Once the data is ready, a call to EfficiencyPlot::Print() will
  generate the formatted plot for each style contained in
  EfficiencyPlot::plot_options_.

*/

/*! \class EfficiencyPlot::SingleEfficiencyPlot

  \brief Container for TH1Ds and TGraphAsymmErrors associated with a single Process

  EfficiencyPlot::SingleEfficiencyPlot is mostly a "dumb" container used by EfficiencyPlot for
  convenience. It contains a pointer to a single Process and some TH1Ds/TGraphAsymmErrors.
*/

#include "core/efficiency_plot.hpp"

#include <cmath>

#include <algorithm>
#include <sstream>

#include <sys/stat.h>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TH1.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TBox.h"
#include "TLegendEntry.h"
#include "TFile.h"

#include "core/plot_opt.hpp"
#include "core/utilities.hpp"

using namespace PlotOptTypes;

/*!\brief Standard constructor

  \param[in] process Process used to fill histogram

  \param[in] hist A fully styled, and typically unfilled histogram to start from
*/
EfficiencyPlot::SingleEfficiencyPlot::SingleEfficiencyPlot(
                                   const EfficiencyPlot &figure,
                                   const std::shared_ptr<Process> &process,
                                   const TH1D &denominator_hist,
                                   const TH1D &numerator_hist):
  FigureComponent(figure, process),
  raw_denominator_hist_(denominator_hist),
  raw_numerator_hist_(numerator_hist),
  proc_and_hist_cut_(figure.cut_ && process->cut_),
  numerator_cut_(figure.numerator_cut_),
  cut_vector_(),
  wgt_vector_(),
  val_vector_(),
  numerator_cut_vector_(){
  raw_numerator_hist_.Sumw2();
  raw_denominator_hist_.Sumw2();
  raw_numerator_hist_.SetBinErrorOption(TH1::kPoisson);
  raw_denominator_hist_.SetBinErrorOption(TH1::kPoisson);
}

void EfficiencyPlot::SingleEfficiencyPlot::RecordEvent(const Baby &baby){
  const EfficiencyPlot& stack = static_cast<const EfficiencyPlot&>(figure_);
  size_t min_vec_size;
  bool have_vec = false;

  const NamedFunc &cut = proc_and_hist_cut_;
  if(cut.IsScalar()){
    if(!cut.GetScalar(baby)) return;
  }else{
    cut_vector_ = cut.GetVector(baby);
    if(!HavePass(cut_vector_)) return;
    have_vec = true;
    min_vec_size = cut_vector_.size();
  }

  if (numerator_cut_.IsVector()) {
    numerator_cut_vector_ = numerator_cut_.GetVector(baby);
    have_vec = true;
    min_vec_size = numerator_cut_vector_.size();
  }

  const NamedFunc &wgt = stack.weight_;
  NamedFunc::ScalarType wgt_scalar = 0.;
  if(wgt.IsScalar()){
    wgt_scalar = wgt.GetScalar(baby);
  }else{
    wgt_vector_ = wgt.GetVector(baby);
    if(!have_vec || wgt_vector_.size() < min_vec_size){
      have_vec = true;
      min_vec_size = wgt_vector_.size();
    }
  }

  const NamedFunc &val = stack.xaxis_.var_;
  NamedFunc::ScalarType val_scalar = 0.;
  if(val.IsScalar()){
    val_scalar = val.GetScalar(baby);
  }else{
    val_vector_ = val.GetVector(baby);
    if(!have_vec || val_vector_.size() < min_vec_size){
      have_vec = true;
      min_vec_size = val_vector_.size();
    }
  }

  if(!have_vec){
    //avoid negative weight events failing numerator since these might make numerator > denominator
    if (wgt_scalar >= 0 || numerator_cut_.GetScalar(baby)) {
      raw_denominator_hist_.Fill(val_scalar, wgt_scalar);
      if (numerator_cut_.GetScalar(baby))
        raw_numerator_hist_.Fill(val_scalar, wgt_scalar);
    }
  }else{
    for(size_t i = 0; i < min_vec_size; ++i){
      if(cut.IsVector() && !cut_vector_.at(i)) continue;
      //avoid negative weight events failing numerator since these might make numerator > denominator
      //TODO: is there some alternative to deal with negative weights?
      if(wgt.IsScalar()) {
        if (numerator_cut_.IsScalar()) {
          if (wgt_scalar < 0 && !numerator_cut_.GetScalar(baby)) continue;
        }
        else {
          if (wgt_scalar < 0 && !numerator_cut_vector_.at(i)) continue;
        }
      }
      else {
        if (numerator_cut_.IsScalar()) {
          if (wgt_vector_.at(i) < 0 && !numerator_cut_.GetScalar(baby)) continue;
        }
        else {
          if (wgt_vector_.at(i) < 0 && !numerator_cut_vector_.at(i)) continue;
        }
      }
      //fill denominator and maybe numerator histograms
      raw_denominator_hist_.Fill(val.IsScalar() ? val_scalar : val_vector_.at(i),
                                 wgt.IsScalar() ? wgt_scalar : wgt_vector_.at(i));
      if (numerator_cut_.IsScalar()) {
        if (!numerator_cut_.GetScalar(baby)) continue;
      }
      else {
        if (!numerator_cut_vector_.at(i)) continue;
      }
      raw_numerator_hist_.Fill(val.IsScalar() ? val_scalar : val_vector_.at(i),
                               wgt.IsScalar() ? wgt_scalar : wgt_vector_.at(i));
    } //vector index
  }
}

void EfficiencyPlot::SingleEfficiencyPlot::BookResult(
    ROOT::RDF::RNode data_frame, int &rdf_plot_idx) {

  const EfficiencyPlot& stack = static_cast<const EfficiencyPlot&>(figure_);
  //no support for vector columns for now
  const NamedFunc &cut = proc_and_hist_cut_;
  const NamedFunc &num = numerator_cut_;
  const NamedFunc &wgt = stack.weight_;
  const NamedFunc &val = stack.xaxis_.var_;
  std::vector<double> bins = stack.xaxis_.Bins();

  std::string den_hist_name = "rdf_hist" + std::to_string(rdf_plot_idx);
  std::string num_hist_name = "rdf_hist" + std::to_string(rdf_plot_idx+1);
  rdf_plot_idx += 2;
  auto figure_filtered_frame = data_frame.Filter(cut.Name());
  booked_raw_denominator_hist_ptr_ = figure_filtered_frame.Histo1D(
      {den_hist_name.c_str(),stack.xaxis_.title_.c_str(),
      static_cast<int>(stack.xaxis_.Nbins()),&bins[0]},val.Name(),wgt.Name());
  figure_filtered_frame = figure_filtered_frame.Filter(num.Name());
  booked_raw_numerator_hist_ptr_ = figure_filtered_frame.Histo1D(
      {num_hist_name.c_str(),stack.xaxis_.title_.c_str(),
      static_cast<int>(stack.xaxis_.Nbins()),&bins[0]},val.Name(),wgt.Name());
}

void EfficiencyPlot::SingleEfficiencyPlot::GetResult() {
  TH1D* booked_raw_numerator_hist_ = (static_cast<TH1D*>(booked_raw_numerator_hist_ptr_->Clone()));
  TH1D* booked_raw_denominator_hist_ = (static_cast<TH1D*>(booked_raw_denominator_hist_ptr_->Clone()));
  raw_numerator_hist_.Add(booked_raw_numerator_hist_);
  raw_denominator_hist_.Add(booked_raw_denominator_hist_);
}

/*! \brief Standard constructor

  \param[in] denominator_cut cut applied to both numerator and denominator of efficiency plot

  \param[in] numerator_cut cut applied to numerator histogram only

  \param[in] processes List of process for the component histograms

  \param[in] definition Specification of contents (plotted variable, binning,
  etc.)

  \param[in] plot_options Styles with which to draw plot
*/
EfficiencyPlot::EfficiencyPlot(const Axis &xaxis, const NamedFunc &denominator_cut,
                               const NamedFunc &numerator_cut,
                               const std::vector<std::shared_ptr<Process> > &processes,
                               const bool draw_histograms,
                               const std::vector<PlotOpt> &plot_options):
  Figure(),
  xaxis_(xaxis),
  cut_(denominator_cut),
  numerator_cut_(numerator_cut),
  plot_options_(plot_options),
  draw_histograms_(draw_histograms),
  weight_("weight"),
  tag_(""),
  ytitle_(""),
  luminosity_tag_(""),
  fixed_title_(""),
  backgrounds_(),
  signals_(),
  datas_(),
  luminosity_(){

  std::string x_title = xaxis_.title_;
  if(xaxis_.units_ != "") x_title += " ["+xaxis_.units_+"]";

  TH1D empty("", (";"+x_title+";").c_str(), xaxis_.Nbins(), &xaxis_.Bins().at(0));
  empty.SetStats(false);
  empty.Sumw2(true);
  bool first_background = true;
  for(const auto &process: processes){
    std::unique_ptr<SingleEfficiencyPlot> eff_plot(new SingleEfficiencyPlot(*this, process, empty, empty));
    eff_plot->raw_denominator_hist_.SetFillColor(process->GetFillColor());
    eff_plot->raw_denominator_hist_.SetFillStyle(process->GetFillStyle());
    eff_plot->raw_denominator_hist_.SetLineColor(process->GetLineColor());
    eff_plot->raw_denominator_hist_.SetLineStyle(process->GetLineStyle());
    eff_plot->raw_denominator_hist_.SetLineWidth(process->GetLineWidth());
    eff_plot->raw_denominator_hist_.SetMarkerColor(process->GetMarkerColor());
    eff_plot->raw_denominator_hist_.SetMarkerStyle(process->GetMarkerStyle());
    eff_plot->raw_denominator_hist_.SetMarkerSize(process->GetMarkerSize());
    eff_plot->raw_numerator_hist_.SetFillColor(process->GetFillColor());
    eff_plot->raw_numerator_hist_.SetFillStyle(process->GetFillStyle());
    eff_plot->raw_numerator_hist_.SetLineColor(process->GetLineColor());
    eff_plot->raw_numerator_hist_.SetLineStyle(process->GetLineStyle());
    eff_plot->raw_numerator_hist_.SetLineWidth(process->GetLineWidth());
    eff_plot->raw_numerator_hist_.SetMarkerColor(process->GetMarkerColor());
    eff_plot->raw_numerator_hist_.SetMarkerStyle(process->GetMarkerStyle());
    eff_plot->raw_numerator_hist_.SetMarkerSize(process->GetMarkerSize());

    switch(process->type_){
    case Process::Type::data:
      datas_.push_back(move(eff_plot));
      data_names_.push_back(process->name_);
      break;
    case Process::Type::background:
      backgrounds_.push_back(move(eff_plot));
      if (first_background) {
        background_color_ = process->GetFillColor();
        first_background = false;
      }
      break;
    case Process::Type::signal:
      signals_.push_back(move(eff_plot));
      signal_names_.push_back(process->name_);
      break;
    default:
      break;
    }
  }
}

/*! \brief Produce and save formatted plots at given luminosity

  \param[in] luminosity The integrated luminosity with which to draw the plot
*/
void EfficiencyPlot::Print(double luminosity,
                   const std::string &subdir){
  //TODO: make EfficiencyPlot::Print as versatile and nice as Hist1D
  luminosity_ = luminosity;

  //gROOT->ForceStyle();
  //RefreshScaledHistos();
  //SetRanges();
  //ApplyStyles();
  //AdjustFillStyles();
  
  //setup canvas/pad
	gStyle->SetOptStat(0);
  std::unique_ptr<TCanvas> full(new TCanvas(("canvas_"+Name()).c_str(),"full",1024,1024));
	full->cd();
  std::unique_ptr<TPad> pad(new TPad("main_pad","pad",0.,0.,1.0,1.0));
	pad->Draw();
	pad->cd();
	pad->SetGrid();
  pad->SetLogy(false);
	gPad->SetMargin(0.15,0.15,0.15,0.15);

  //draw ratio plots and histograms
  double ratio_y_max = 0.0;
  double linear_y_max = 0.0;
  std::unique_ptr<TH1D> total_background_denominator;
  std::unique_ptr<TH1D> total_background_numerator;
  std::vector<std::unique_ptr<TH1D>> signal_numerator_plots;
  std::vector<std::unique_ptr<TH1D>> signal_denominator_plots;
  std::vector<Color_t> signal_colors;
  std::vector<std::unique_ptr<TH1D>> data_numerator_plots;
  std::vector<std::unique_ptr<TH1D>> data_denominator_plots;
  std::vector<Color_t> data_colors;

  //first loop generates ratio plots and finds maximum for scaling
  if (backgrounds_.size() > 0) {
    total_background_denominator = std::unique_ptr<TH1D>(
        static_cast<TH1D*>(backgrounds_.at(0)->raw_denominator_hist_.Clone()));
    total_background_numerator = std::unique_ptr<TH1D>(
        static_cast<TH1D*>(backgrounds_.at(0)->raw_numerator_hist_.Clone()));
    for (unsigned int bkgd_idx = 0; bkgd_idx < backgrounds_.size(); bkgd_idx++) {
      total_background_denominator->Add(&(backgrounds_.at(bkgd_idx)->raw_denominator_hist_));
      total_background_numerator->Add(&(backgrounds_.at(bkgd_idx)->raw_numerator_hist_));
    }
    total_background_numerator->Scale(luminosity_);
    total_background_denominator->Scale(luminosity_);
    //double this_linear_y_max = total_background_denominator->GetMaximum();
    double this_linear_y_max = total_background_denominator->GetBinContent(total_background_denominator->GetMaximumBin());
    if (this_linear_y_max > linear_y_max) linear_y_max = this_linear_y_max;
    TFile* out_file = TFile::Open("ntuples/effplotdebugging.root","UPDATE");
    total_background_numerator->Write((Name()+"_num").c_str());
    total_background_denominator->Write((Name()+"_den").c_str());
    out_file->Close();
    background_ratio_plot_ = std::unique_ptr<TGraphAsymmErrors>(
        new TGraphAsymmErrors(total_background_numerator.get(),total_background_denominator.get(),"cp"));
    if (background_ratio_plot_->GetN() > 0) {
      double this_ratio_y_max = TMath::MaxElement(background_ratio_plot_->GetN(), background_ratio_plot_->GetY());
      if (this_ratio_y_max > ratio_y_max) ratio_y_max = this_ratio_y_max;
    }
  }
  for (unsigned int signal_idx = 0; signal_idx < signals_.size(); signal_idx++) {
    TH1D* numerator_hist = &signals_.at(signal_idx)->raw_numerator_hist_;
    TH1D* denominator_hist = &signals_.at(signal_idx)->raw_denominator_hist_;
    denominator_hist->Scale(luminosity_);
    numerator_hist->Scale(luminosity_);
    //double this_linear_y_max = denominator_hist->GetMaximum();
    double this_linear_y_max = denominator_hist->GetBinContent(denominator_hist->GetMaximumBin());
    if (this_linear_y_max > linear_y_max) linear_y_max = this_linear_y_max;
    signal_colors.push_back(numerator_hist->GetLineColor());
    signal_numerator_plots.push_back(
        std::unique_ptr<TH1D>(static_cast<TH1D*>(numerator_hist->Clone())));
    signal_denominator_plots.push_back(
        std::unique_ptr<TH1D>(static_cast<TH1D*>(denominator_hist->Clone())));
    signal_ratio_plots_.push_back(
        std::unique_ptr<TGraphAsymmErrors>(new TGraphAsymmErrors(numerator_hist,denominator_hist,"cp")));
    if (signal_ratio_plots_.back()->GetN() > 0) {
      double this_ratio_y_max = TMath::MaxElement(signal_ratio_plots_.back()->GetN(), signal_ratio_plots_.back()->GetY());
      if (this_ratio_y_max > ratio_y_max) ratio_y_max = this_ratio_y_max;
    }
  }
  for (unsigned int data_idx = 0; data_idx < datas_.size(); data_idx++) {
    TH1D* numerator_hist = &datas_.at(data_idx)->raw_numerator_hist_;
    TH1D* denominator_hist = &datas_.at(data_idx)->raw_denominator_hist_;
    //double this_linear_y_max = denominator_hist->GetMaximum();
    double this_linear_y_max = denominator_hist->GetBinContent(denominator_hist->GetMaximumBin());
    if (this_linear_y_max > linear_y_max) linear_y_max = this_linear_y_max;
    data_colors.push_back(numerator_hist->GetLineColor());
    data_numerator_plots.push_back(
        std::unique_ptr<TH1D>(static_cast<TH1D*>(numerator_hist->Clone())));
    data_denominator_plots.push_back(
        std::unique_ptr<TH1D>(static_cast<TH1D*>(denominator_hist->Clone())));
    data_ratio_plots_.push_back(
        std::unique_ptr<TGraphAsymmErrors>(new TGraphAsymmErrors(numerator_hist,denominator_hist,"cp")));
    if (data_ratio_plots_.back()->GetN() > 0) {
      double this_ratio_y_max = TMath::MaxElement(data_ratio_plots_.back()->GetN(), data_ratio_plots_.back()->GetY());
      if (this_ratio_y_max > ratio_y_max) ratio_y_max = this_ratio_y_max;
    }
  }

  //second loop: draw
  std::string draw_options = "AP"; //change after axes are drawn once
  if (backgrounds_.size() > 0) {
    SetRatioPlotDrawOptions(background_ratio_plot_, 0, ratio_y_max);
    background_ratio_plot_->SetMarkerColor(background_color_);
    background_ratio_plot_->SetLineColor(background_color_);
    background_ratio_plot_->Draw(draw_options.c_str());
    draw_options = "P";
    if (draw_histograms_) {
      SetLinearPlotDrawOptions(total_background_denominator, linear_y_max/ratio_y_max, true);
      SetLinearPlotDrawOptions(total_background_numerator, linear_y_max/ratio_y_max, false);
      total_background_denominator->SetLineColor(background_color_);
      total_background_numerator->SetLineColor(background_color_);
      total_background_denominator->Draw("same hist");
      total_background_numerator->Draw("same hist");
    }
  }
  for (unsigned int signal_idx = 0; signal_idx < signals_.size(); signal_idx++) {
    SetRatioPlotDrawOptions(signal_ratio_plots_.at(signal_idx), 0, ratio_y_max);
    signal_ratio_plots_.at(signal_idx)->SetMarkerColor(signal_colors.at(signal_idx));
    signal_ratio_plots_.at(signal_idx)->SetLineColor(signal_colors.at(signal_idx));
    signal_ratio_plots_.at(signal_idx)->Draw(draw_options.c_str());
    draw_options = "P";
    if (draw_histograms_) {
      SetLinearPlotDrawOptions(signal_denominator_plots.at(signal_idx), linear_y_max/ratio_y_max, true);
      SetLinearPlotDrawOptions(signal_numerator_plots.at(signal_idx), linear_y_max/ratio_y_max, false);
      signal_denominator_plots.at(signal_idx)->SetLineColor(signal_colors.at(signal_idx));
      signal_numerator_plots.at(signal_idx)->SetLineColor(signal_colors.at(signal_idx));
      signal_denominator_plots.at(signal_idx)->Draw("same hist");
      signal_numerator_plots.at(signal_idx)->Draw("same hist");
    }
  }
  for (unsigned int data_idx = 0; data_idx < datas_.size(); data_idx++) {
    SetRatioPlotDrawOptions(data_ratio_plots_.at(data_idx), 0, ratio_y_max);
    data_ratio_plots_.at(data_idx)->SetMarkerColor(data_colors.at(data_idx));
    data_ratio_plots_.at(data_idx)->SetLineColor(data_colors.at(data_idx));
    data_ratio_plots_.at(data_idx)->Draw(draw_options.c_str());
    draw_options = "P";
    if (draw_histograms_) {
      SetLinearPlotDrawOptions(data_denominator_plots.at(data_idx), linear_y_max/ratio_y_max, true);
      SetLinearPlotDrawOptions(data_numerator_plots.at(data_idx), linear_y_max/ratio_y_max, false);
      data_denominator_plots.at(data_idx)->SetLineColor(data_colors.at(data_idx));
      data_numerator_plots.at(data_idx)->SetLineColor(data_colors.at(data_idx));
      data_denominator_plots.at(data_idx)->Draw("same hist");
      data_numerator_plots.at(data_idx)->Draw("same hist");
    }
  }
  //draw ratio plots once more to get them on top, TODO: fix this
  if (backgrounds_.size() > 0) {
    background_ratio_plot_->Draw(draw_options.c_str());
  }
  for (unsigned int signal_idx = 0; signal_idx < signals_.size(); signal_idx++) {
    signal_ratio_plots_.at(signal_idx)->Draw(draw_options.c_str());
  }
  for (unsigned int data_idx = 0; data_idx < datas_.size(); data_idx++) {
    data_ratio_plots_.at(data_idx)->Draw(draw_options.c_str());
  }

  //draw overlay (title) text
	TLatex t;
	t.SetTextColor(kBlack);
	t.SetTextSize(0.04);
  PlotOpt this_opt_ = plot_options_[0];
  //TODO: fix this to allow multiple styles in one go
  if (this_opt_.Title() == TitleType::preliminary)
    t.DrawLatexNDC(0.155,0.87,"#font[62]{CMS} #scale[0.8]{#font[52]{Preliminary}}");
  else if (this_opt_.Title() == TitleType::simulation_preliminary)
    t.DrawLatexNDC(0.155,0.87,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Preliminiary}}");
  else if (this_opt_.Title() == TitleType::simulation)
    t.DrawLatexNDC(0.155,0.87,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation}}");
  else if (this_opt_.Title() == TitleType::supplementary)
    t.DrawLatexNDC(0.155,0.87,"#font[62]{CMS} #scale[0.8]{#font[52]{Supplementary}}");
  else if (this_opt_.Title() == TitleType::simulation_supplementary)
    t.DrawLatexNDC(0.155,0.87,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Supplementary}}");
  else if (this_opt_.Title() == TitleType::data)
    t.DrawLatexNDC(0.155,0.87,"#font[62]{CMS}");
  else if (this_opt_.Title() == TitleType::info) {
    if (datas_.size() == 0) {
      t.DrawLatexNDC(0.155,0.87,"#font[62]{CMS} #scale[0.8]{#font[52]{Simulation Preliminiary}}");
    }
    else {
      t.DrawLatexNDC(0.155,0.87,"#font[62]{CMS} #scale[0.8]{#font[52]{Preliminary}}");
    }
    //TODO: add info title
  }
	t.SetTextAlign(31);
  TString lumi_string = RoundNumber(luminosity_,1) + " fb^{-1}";
  if (luminosity_tag_ != "")
    lumi_string = luminosity_tag_ + " fb^{-1}";
	t.DrawLatexNDC(0.845,0.87,("#font[42]{"+lumi_string+"}").Data());
	t.SetTextAlign(33);
	t.SetTextSize(0.025);
  if (Title().size() >= 66) t.SetTextSize(0.015);
  if (Title().size() >= 120) t.SetTextSize(0.015*120.0/static_cast<float>(Title().size()));
	t.DrawLatexNDC(0.825,0.83,("#font[42]{"+Title()+"}").c_str());
	//if (Title().size() < 66) {
	//	t.DrawLatexNDC(0.825,0.83,("#font[42]{"+Title()+"}").c_str());
	//}
	//else {
	//	//title too long, draw on separate lines
	//	unsigned int right_split_comma_pos = 999;
	//	for (unsigned int string_pos = 65; string_pos > 0; string_pos --) {
	//		if (title[string_pos]==',') {
	//			//split title at this comma
	//			right_split_comma_pos = string_pos;
	//			break;
	//		}
	//	}
	//	if (right_split_comma_pos == 999) {
	//		//unable to find splitting comma, just let it draw off histogram
	//		t.DrawLatexNDC(0.825,0.83,("#font[42]{"+title+"}").c_str());
	//	}
	//	else {
	//		t.DrawLatexNDC(0.825,0.83,("#font[42]{"+title.substr(0,right_split_comma_pos+1)+"}").c_str());
	//		t.DrawLatexNDC(0.825,0.78,("#font[42]{"+title.substr(right_split_comma_pos+1,title.size()-right_split_comma_pos-1)+"}").c_str());
	//	}
	//}

	//draw right axis
	pad->Update();
	TGaxis *norm_axis = new TGaxis(gPad->GetUxmax(), gPad->GetUymin(), gPad->GetUxmax(), 
                                 gPad->GetUymax(), 0, gPad->GetUymax()*2.0/ratio_y_max*linear_y_max, 505, "+L");
  norm_axis->SetTickLength(0.3);
  norm_axis->SetLabelSize(0.03);
  norm_axis->SetTitle("Events/bin");
  //norm_axis->SetTitleColor(kBlue);
  norm_axis->SetTitleFont(42);
  norm_axis->SetTitleOffset(1.6);
  //norm_axis->SetLineColor(kBlue);
  //norm_axis->SetLabelColor(kBlue);
  norm_axis->Draw();
    
  //draw legend
  unsigned int number_legend_entries = signal_names_.size()+data_names_.size();
  if (backgrounds_.size() > 0) number_legend_entries++;
  TLegend * legend = new TLegend(0.35, 0.65, 0.8, 0.8);
  legend->SetNColumns(3);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  if (backgrounds_.size() > 0) {
    legend->AddEntry(background_ratio_plot_.get(),"Background","l");
  }
  for (unsigned int signal_idx = 0; signal_idx < signals_.size(); signal_idx++) {
    legend->AddEntry(signal_ratio_plots_.at(signal_idx).get(),signal_names_.at(signal_idx).c_str(),"l");
  }
  for (unsigned int data_idx = 0; data_idx < datas_.size(); data_idx++) {
    legend->AddEntry(data_ratio_plots_.at(data_idx).get(),data_names_.at(data_idx).c_str(),"l");
  }
  if (number_legend_entries < 12) {
    for (unsigned int extra_entries_idx = number_legend_entries; extra_entries_idx < 12; extra_entries_idx++) {
      legend->AddEntry(static_cast<TObject*>(0), "", "");
    }
  }
  legend->Draw();

  //TODO: implement GetCutLines
  //std::vector<TLine> cut_vals = GetCutLines(GetMinDraw(), GetMaxDraw());
  //for(auto &cut: cut_vals) cut.Draw();

  if(subdir != "") mkdir(("plots/"+subdir).c_str(), 0777);
  std::string base_name = subdir != ""
    ? "plots/"+subdir+"/"+Name()
    : "plots/"+Name();
  //TODO: support other file extensions
  std::string full_name = base_name+".pdf";
  if (Contains(tag_,"FixName:")){
    std::string tagName=tag_;
    ReplaceAll(tagName, "FixName:", "");
    base_name = subdir != ""
      ? "plots/"+subdir+"/"+tagName
      : "plots/"+tagName;
    full_name = base_name+".pdf";
  }
  full->Print(full_name.c_str());
  std::cout << "open " << full_name << std::endl;
}

std::string EfficiencyPlot::GetTag() const{
  return tag_;
}

std::set<const Process*> EfficiencyPlot::GetProcesses() const{
  std::set<const Process*> processes;
  for(const auto &proc: backgrounds_){
    processes.insert(proc->process_.get());
  }
  for(const auto &proc: signals_){
    processes.insert(proc->process_.get());
  }
  for(const auto &proc: datas_){
    processes.insert(proc->process_.get());
  }
  return processes;
}

Figure::FigureComponent * EfficiencyPlot::GetComponent(const Process *process){
  const auto &component_list = GetComponentList(process);
  for(const auto &component: component_list){
    if(component->process_.get() == process){
      return component.get();
    }
  }
  DBG("Could not find histogram for process "+process->name_+".");
  return nullptr;
}

std::string EfficiencyPlot::Name() const{
  std::string cut = "";
  if(cut_.Name() != "1") cut = "__"+cut_.Name();

  std::string numerator_cut = "";
  if(numerator_cut_.Name() != "1") numerator_cut = "__"+numerator_cut_.Name();

  std::string weight = "";
  if(weight_.Name() != "weight") weight = "__"+weight_.Name();
  
  if(tag_ == ""){
    return CodeToPlainText(xaxis_.var_.Name()+cut+numerator_cut+weight);
  }else if (Contains(tag_,"ShortName:")){
    std::string tagName=tag_;
    ReplaceAll(tagName, "ShortName:", "");
    return CodeToPlainText(tagName+"__"+xaxis_.var_.Name()+weight);
  }else{
    return CodeToPlainText(tag_+"__"+xaxis_.var_.Name()+cut+numerator_cut+weight);
  }
}

std::string EfficiencyPlot::Title() const{
  if (fixed_title_ != "") return fixed_title_;
  bool cut = (cut_.Name() != "" && cut_.Name() != "1");
  bool weight = weight_.Name() != "weight";
  if(cut && weight){
    return CodeToRootTex(cut_.Name())+" (weight="+CodeToRootTex(weight_.Name())+")";
  }else if(cut){
    return CodeToRootTex(cut_.Name());
  }else if(weight){
    return CodeToRootTex("weight="+weight_.Name());
  }else{
    return "";
  }
}

std::string EfficiencyPlot::YAxisText() const{
  std::string y_axis_title = "Efficiency ";
  if (ytitle_ != "")
    y_axis_title += ytitle_;
  else
    y_axis_title += numerator_cut_.Name();
  //y_axis_title += CodeToRootTex(numerator_cut_.Name());
  return y_axis_title;
}

EfficiencyPlot & EfficiencyPlot::Weight(const NamedFunc &weight){
  weight_ = weight;
  return *this;
}

EfficiencyPlot & EfficiencyPlot::Tag(const std::string &tag){
  tag_ = tag;
  return *this;
}

EfficiencyPlot & EfficiencyPlot::YTitle(const std::string &ytitle){
  ytitle_ = ytitle;
  return *this;
}

EfficiencyPlot & EfficiencyPlot::LuminosityTag(const std::string &tag){
  luminosity_tag_ = tag;
  return *this;
}

EfficiencyPlot & EfficiencyPlot::FixTitle(const std::string &title){
  fixed_title_ = title;
  return *this;
}

void EfficiencyPlot::SetRatioPlotDrawOptions(std::unique_ptr<TGraphAsymmErrors> &ratio_plot, 
                                             double ymin, double ymax) {
  double xmin = xaxis_.Bins().at(0);
  double xmax = xaxis_.Bins().at(xaxis_.Bins().size()-1);
  std::string yaxis_text = YAxisText();
  ratio_plot->SetMarkerStyle(kFullCircle);
  ratio_plot->SetMarkerSize(1);
  ratio_plot->SetLineWidth(2);
  //ratio_plot->GetXaxis()->SetRangeUser(xmin,xmax);
  ratio_plot->GetXaxis()->SetLimits(xmin,xmax);
  ratio_plot->GetXaxis()->SetTitleSize(0.04);
  ratio_plot->GetXaxis()->SetTitleOffset(1.0);
  ratio_plot->GetYaxis()->SetRangeUser(ymin, 1.4*ymax);
  if (yaxis_text.size() > 40) {
    ratio_plot->GetYaxis()->SetTitleSize(0.025);
    ratio_plot->GetYaxis()->SetTitleOffset(1.6);
  }
  else {
    ratio_plot->GetYaxis()->SetTitleSize(0.04);
    ratio_plot->GetYaxis()->SetTitleOffset(1.4);
  }
  ratio_plot->SetTitle((";"+xaxis_.Title()+";"+yaxis_text).c_str());
}

void EfficiencyPlot::SetLinearPlotDrawOptions(std::unique_ptr<TH1D> &linear_plot, 
                                double ymax, bool is_denominator) {
  linear_plot->Scale(0.5/ymax);
  if (is_denominator) linear_plot->SetLineStyle(2);
  else linear_plot->SetLineStyle(1);
  linear_plot->SetLineWidth(2);
  linear_plot->SetFillStyle(0);
}



///*!\brief Get vertical lines at cut values
//
//  \param[in] y_min Lower bound of y-axis
//
//  \param[in] y_max Upper bound of y-axis
//
//  \return Lines at x-coordinate of cut value and y-coordinates running from
//  bottom of plot to bottom of legend
//*/
//vector<TLine> EfficiencyPlot::GetCutLines(double y_min, double y_max, bool adjust_bottom) const{
//  double bottom = y_min;
//  if(adjust_bottom){
//    switch(this_opt_.YAxis()){
//    default:
//      DBG("Bad YAxis type " << static_cast<int>(this_opt_.YAxis()));
//      /* FALLTHRU */
//    case YAxisType::linear: bottom = y_min >= 0. ? 0. : y_min; break;
//    case YAxisType::log:    bottom = y_min > this_opt_.LogMinimum() ? y_min : this_opt_.LogMinimum(); break;
//    }
//  }
//  vector<TLine> out(xaxis_.cut_vals_.size());
//  for(double cut: xaxis_.cut_vals_){
//    out.emplace_back(cut, bottom, cut, y_max);
//    out.back().SetNDC(false);
//    out.back().SetLineStyle(2);
//    out.back().SetLineColor(kBlack);
//    out.back().SetLineWidth(3);
//  }
//
//  return out;
//}

const std::vector<std::unique_ptr<EfficiencyPlot::SingleEfficiencyPlot> >& EfficiencyPlot::GetComponentList(const Process *process){
  switch(process->type_){
  case Process::Type::data:
    return datas_;
  case Process::Type::background:
    return backgrounds_;
  case Process::Type::signal:
    return signals_;
  default:
    ERROR("Did not understand process type "+std::to_string(static_cast<long>(process->type_))+".");
    return backgrounds_;
  }
}
