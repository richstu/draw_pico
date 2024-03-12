/*! \class Datacard
 
  \brief A datacard for use with Combine

  Datacard contains machinery to create a Higgs combine datacard

  TODO: add more description here as this class is developed
  TODO: add usage instructions for users
  TODO: add description for datacardprocess subclass

*/

/*! \class Datacard::DatacardProcess

  \brief Container for TH1Ds(->RooDataHist?) associated with a single Process
  
  TODO: may be modified if we need RooDataSet instead
*/

/*! \class Datacard::Systematic

  \brief Container for systematics (nuisance parameters) other than PDF parameters
  
  These come in two flavours: weight systematics where the nominal weight is 
  replaced by an alternative, and selection systematics where one selection
  is replaced by an alternative
*/

/*! \class Datacard::SelectionList

  \brief Container for selections to be applied to a particular category
*/

#include "core/datacard.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/stat.h>

#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooHistPdf.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooCategory.h"
#include "TFile.h"
#include "TH1.h"

#include "core/axis.hpp"
#include "core/baby.hpp"
#include "core/figure.hpp"
#include "core/named_func.hpp"
#include "core/process.hpp"
#include "core/RooMultiPdf.hpp"

//----------------------------------------------------------------------------
//SelectionList
//----------------------------------------------------------------------------

/*!\brief Standard constructor
*/
Datacard::SelectionList::SelectionList(const std::string &channel_name) :
    channel_name_(channel_name) {}

/*!\brief Method to add a selection
 \param[in] name       Selection name, needed for systematics referencing
 \param[in] selection  Selection
*/
Datacard::SelectionList & Datacard::SelectionList::AddSelection(
    const std::string &name, const NamedFunc &selection) {
  if (!selection.IsScalar()) {
    throw std::invalid_argument(("Selection NamedFunc "+selection.Name()
                                 +" is not a scalar.").c_str());
  }
  name_.push_back(name);
  selection_.push_back(selection);
  return *this;
}

//----------------------------------------------------------------------------
//Systematic
//----------------------------------------------------------------------------

/*!\brief Constructor for reweighting-based systematics
 \param[in] name              Name of systematic
 \param[in] alternate_weight  Alternate weight
*/
Datacard::Systematic::Systematic(const std::string &name, 
                                 const NamedFunc &alternate_weight) :
  is_weight_systematic_(true),
  name_(name),
  selection_name_(""),
  content_(alternate_weight) {
  if (!alternate_weight.IsScalar()) {
    throw std::invalid_argument(("Systematic NamedFunc "+alternate_weight.Name()
                                 +" is not a scalar.").c_str());
  }
}

/*!\brief Constructor for reweighting-based systematics
 \param[in] name                  Name of systematic
 \param[in] selection_name        Name of selection to replace
 \param[in] alternate_selection   Alternate selection
*/
Datacard::Systematic::Systematic(const std::string &name,
                                 const std::string &selection_name, 
                                 const NamedFunc &alternate_selection) :
  is_weight_systematic_(false),
  name_(name),
  selection_name_(selection_name),
  content_(alternate_selection) {
  if (!alternate_selection.IsScalar()) {
    throw std::invalid_argument(("Systematic NamedFunc "
                                 +alternate_selection.Name()
                                 +" is not a scalar.").c_str());
  }
}

//----------------------------------------------------------------------------
//DatacardProcess
//----------------------------------------------------------------------------

/*!\brief Standard constructor
  \param[in] figure   Parent figure
  \param[in] process  Process used to fill histogram
*/
Datacard::DatacardProcess::DatacardProcess(const Figure &figure,
    const std::shared_ptr<Process> &process,
    const Axis &axis) :
    FigureComponent(figure, process),
    weight_(RooRealVar("weight","",-50.0,50.0)) {
  if (process_->type_ == Process::Type::data) {
    is_data_ = true;
    process_->name_ = "data_obs";
  }
  else {
    is_data_ = false;
  }
  replace_with_param_ = false;
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  float axis_min = axis.Bins().at(0);
  float axis_max = axis.Bins().at(axis.Bins().size()-1);
  for (unsigned ichan = 0; ichan < datacard->n_channels_; ichan++) {
    var_.push_back(RooRealVar(
        (axis.var_.Name()+"_"+datacard->channel_name_[ichan]).c_str(),"",
        axis.Bins().at(0),axis.Bins().at(axis.Bins().size()-1)));
    raw_histogram_nom_.push_back(TH1D("","",axis.Nbins(),&axis.Bins().at(0)));
    raw_histogram_sys_.push_back(std::vector<TH1D>());
    if (is_data_)
      raw_dataset_nom_.push_back(RooDataSet("","",RooArgSet(var_[ichan],weight_),RooFit::WeightVar(weight_)));

    //for now, only do normalization systematics
    for (unsigned isyst = 0; isyst < datacard->n_systematics_; isyst++) {
      raw_histogram_sys_.back().push_back(TH1D("","",1,axis_min,axis_max));
    }
  }
}

void Datacard::DatacardProcess::RecordEvent(const Baby &baby) {
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  float nominal_weight = datacard->nominal_weight_.GetScalar(baby);
  for (unsigned ichan = 0; ichan < datacard->n_channels_; ichan++) {
    if (datacard->channel_selection_[ichan].GetScalar(baby)) {
      float nominal_value = datacard->variable_.GetScalar(baby);
      raw_histogram_nom_[ichan].Fill(nominal_value, nominal_weight);
      if (is_data_) {
        var_[ichan].setVal(nominal_value);
        raw_dataset_nom_[ichan].add(RooArgSet(var_[ichan]),nominal_weight);
      }
    }
    for (unsigned isyst = 0; isyst < datacard->n_systematics_; isyst++) {
      if (datacard->systematic_selection_[ichan][isyst].GetScalar(baby)) {
        raw_histogram_sys_[ichan][isyst].Fill(
            datacard->variable_.GetScalar(baby),
            datacard->systematic_weight_[ichan][isyst].GetScalar(baby));
      }
    } //loop over systematics
  } //loop over channels
}

/*! \brief Returns the name of the RooWorkspace saved for a given channel
  \param[in] channel    index of channel to write
*/
std::string Datacard::DatacardProcessNonparametric::WSName(unsigned int channel) {
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  if (is_data_) {
    return "WS_data_obs_"+datacard->channel_name_[channel];
  }
  //is MC
  return "WS_"+name_+"_"+datacard->channel_name_[channel];
}

/*! \brief Returns the name of the PDF/RooDataHist/RooDataSet saved for a given channel
  \param[in] channel    index of channel to write
*/
std::string Datacard::DatacardProcessNonparametric::PDFName(unsigned int channel) {
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  if (is_data_) {
    return "data_obs_"+datacard->channel_name_[channel];
  }
  else if (replace_with_param_) { //MC replaced by parametric PDF
    return param_pdf_[channel]->GetName();
  }
  //is MC RooHistPdf
  return "pdf_"+name_+"_"+datacard->channel_name_[channel];
}

/*! \brief Returns yield (rate) of events in a given channel
  \param[in] channel    index of channel
*/
float Datacard::DatacardProcessNonparametric::Yield(unsigned int channel, 
                                                    unsigned int systematic) {
  if (systematic==999)
    return raw_histogram_nom_[channel].Integral();
  float ratio = raw_histogram_nom_[channel].Integral()/
                raw_histogram_sys_[channel][systematic].Integral();
  if (ratio<1.0) ratio = 1.0/ratio;
  return ratio;
}

/*! \brief Performs fit to data and freezes PDF parameters
  \param[in] channel   index of channel to fit
*/
void Datacard::DatacardProcessNonparametric::FitAndFreeze(unsigned int channel) {
  if (is_data_ || !replace_with_param_)
    throw std::runtime_error("FitAndFreeze called for function with no parametric form.");
  RooDataHist fit_hist("fit_hist","",RooArgList(var_[channel]), &raw_histogram_nom_[channel]);
  param_pdf_[channel]->fitTo(fit_hist);
  RooArgSet* params = param_pdf_[channel]->getParameters(fit_hist);
  //can't use normal range for loop because ancient ROOT version...
  RooFIter param_iterator = params->fwdIterator();
  bool param_loop = true;
  while (param_loop) {
    RooAbsArg* param = param_iterator.next();
    if (param != nullptr) {
      static_cast<RooRealVar*>(param)->setConstant(true);
    }
    else {
      param_loop = false;
    }
  }
}

/*! \brief Writes contained PDFs to a workspace and saves to ROOT file. 
    Assumes ROOT file open
  \param[in] channel    index of channel to write
*/
void Datacard::DatacardProcessNonparametric::WriteWorkspace(unsigned int channel) {
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  RooWorkspace ws(WSName(channel).c_str());
  if (is_data_) {
    if (datacard->save_data_as_hist_) {
      RooDataHist data_obs(PDFName(channel).c_str(),"",RooArgList(var_[channel]),
          &raw_histogram_nom_[channel]);
      ws.import(data_obs);
    }
    else { //save data as RooDataSet
      raw_dataset_nom_[channel].SetName(PDFName(channel).c_str());
      ws.import(raw_dataset_nom_[channel]);
    }
  }
  else { //is MC
    if (replace_with_param_) { //replace histogram with analytic PDF
      FitAndFreeze(channel);
      ws.import(*param_pdf_[channel]);
    }
    else { //don't replace with param, just make RooHistPdf
      RooDataHist hist((PDFName(channel)+"_hist").c_str(),"",RooArgList(var_[channel]),
          &raw_histogram_nom_[channel]);
      RooHistPdf pdf(PDFName(channel).c_str(),"",RooArgList(var_[channel]),hist);
      //add a separate norm RooRealVar for non-signal processes (signal controlled by r)
      if (!is_signal_) {
        float nom_yield = raw_histogram_nom_[channel].Integral();
        RooRealVar norm((PDFName(channel)+"_norm").c_str(),"",nom_yield,
            0,10.0*nom_yield);
        ws.import(norm);
      }
      ws.import(pdf);
    }
  }
  ws.Write();
}

//----------------------------------------------------------------------------
//DatacardProcessParametric
//----------------------------------------------------------------------------

/*! \brief Add parametric process to datacard (ex. background constrained by fit)
  \param[in] name   Name of process
  \param[in] pdf    PDF for process by channel
*/
Datacard::DatacardProcessParametric::DatacardProcessParametric(
    const std::string &name, 
    std::vector<std::shared_ptr<RooAbsPdf>> &pdf, const Figure &figure) :
    DatacardProcess(figure) {

  std::vector<std::vector<std::shared_ptr<RooAbsPdf>>> double_vec;
  for (std::shared_ptr<RooAbsPdf> channel_pdf : pdf) {
    double_vec.push_back(std::vector<std::shared_ptr<RooAbsPdf>>());
    double_vec.back().push_back(channel_pdf);
    is_profiled_.push_back(false);
  }

  //initialize properties
  name_ = name;
  pdf_ = double_vec;
  is_data_ = false;
  is_signal_ = false; //currently, don't allow signal
}

/*! \brief Add parametric process to datacard (ex. background constrained by fit)
  \param[in] name   Name of process
  \param[in] pdf    PDF for process by channel and profile
*/
Datacard::DatacardProcessParametric::DatacardProcessParametric(
    const std::string &name, 
    std::vector<std::vector<std::shared_ptr<RooAbsPdf>>> &pdf, 
    const Figure &figure) : 
    DatacardProcess(figure) {

  //check for errors and which channels are profiled
  for (std::vector<std::shared_ptr<RooAbsPdf>>& channel_pdfs : pdf) {
    if (channel_pdfs.size()==0)
      throw std::invalid_argument(("Process "+name+" has a channel with 0 PDFs.").c_str());
    else if (channel_pdfs.size()==1)
      is_profiled_.push_back(false);
    else
      is_profiled_.push_back(true);
  }

  //initialize properties
  name_ = name;
  pdf_ = pdf;
  is_data_ = false;
  is_signal_ = false; //currently, don't allow signal
}


/*! \brief Helper function that wraps a vector of PDFs as a vector of single-entry vectors of PDFs
  \param[in] pdf  vector of (pointers to) PDFs
*/
std::vector<std::vector<std::shared_ptr<RooAbsPdf>>> 
    Datacard::DatacardProcessParametric::MakeDoubleVector(
    std::vector<std::shared_ptr<RooAbsPdf>> pdfs) {
  std::vector<std::vector<std::shared_ptr<RooAbsPdf>>> double_vec;
  for (std::shared_ptr<RooAbsPdf> pdf : pdfs) {
    double_vec.push_back(std::vector<std::shared_ptr<RooAbsPdf>>());
    double_vec.back().push_back(pdf);
  }
  return double_vec;
}

/*! \brief Dummy function needed for interface
*/
void Datacard::DatacardProcessParametric::RecordEvent(const Baby &baby) {
  baby.SampleType();
  return;
}

/*! \brief Returns the name of the RooAbsPdf saved for a given channel
  \param[in] channel    index of channel
*/
std::string Datacard::DatacardProcessParametric::PDFName(unsigned int channel) {
  if (!is_profiled_[channel]) {
    return pdf_[channel][0]->GetName();
  }
  //is profiled
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  return "pdf_"+name_+"_"+datacard->channel_name_[channel];
}

/*! \brief Returns the name of the RooWorkspace saved for a given channel
  \param[in] channel    index of channel
*/
std::string Datacard::DatacardProcessParametric::WSName(unsigned int channel) {
  //if (figure_==nullptr)
  //  throw std::runtime_error(("WSName called for "+name_+" before datacard assignment.").c_str());
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  //is MC
  return "WS_"+name_+"_"+datacard->channel_name_[channel];
}

/*! \brief Dummy function that returns a yield (rate) of 1.0 for parametric 
    function as per combine convention
*/
float Datacard::DatacardProcessParametric::Yield(unsigned int channel, unsigned int systematic) {
  channel += 0;
  systematic += 0;
  return 1.0;
}

/*! \brief Writes contained PDFs to a workspace and saves to ROOT file. 
    Assumes ROOT file open
  \param[in] channel    index of channel to write
*/
void Datacard::DatacardProcessParametric::WriteWorkspace(unsigned int channel) {
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  RooWorkspace ws(WSName(channel).c_str());
  if (!is_profiled_[channel]) {
    RooRealVar norm((PDFName(channel)+"_norm").c_str(),"",data_norm_[channel],
        0,3.0*data_norm_[channel]);
    ws.import(*pdf_[channel][0]);
    ws.import(norm);
  }
  else {
    std::string index_name = "pdfindex_" + name_ +"_" + datacard->channel_name_[channel];
    std::string model_name = "pdfmodels_" + name_ + "_"+ datacard->channel_name_[channel];
    RooRealVar norm((PDFName(channel)+"_norm").c_str(),"",data_norm_[channel],
        0,3.0*data_norm_[channel]);
    RooCategory pdfindex(index_name.c_str(), index_name.c_str());
    RooArgList models = RooArgList(model_name.c_str());
    for (std::shared_ptr<RooAbsPdf> pdf: pdf_[channel]){
      models.add(*pdf);
    }
    RooMultiPdf profile(PDFName(channel).c_str(), PDFName(channel).c_str(), pdfindex, models);
    //ws.import(*(param_profile_ind_process_[iproc_eff][channel]));
    ws.import(pdfindex);
    ws.import(profile);
    ws.import(norm);
  }
  ws.Write();
}

//----------------------------------------------------------------------------
//Datacard
//----------------------------------------------------------------------------

/*! \brief Constructor for shape-based datacard. Cut-and-count constructor 
           still to be implemented
  \param[in] name         datacard name
  \param[in] channels     selection list for each channel
  \param[in] systematics  list of systematics
  \param[in] processes    standard processes including data, signal, and 
                          MC-based backgrounds
  \param[in] weight       nominal weight
  \param[in] axis         associated histogram axis
*/
Datacard::Datacard(const std::string &name,
                   const std::vector<SelectionList> &channels, 
                   const std::vector<Systematic> &systematics,
                   const std::vector<std::shared_ptr<Process>> &processes,
                   const NamedFunc &weight,
                   const Axis &axis) :
    name_(name),
    n_channels_(channels.size()),
    n_processes_(processes.size()),
    n_systematics_(systematics.size()),
    nominal_weight_(weight),
    variable_(axis.var_),
    channel_selection_(),
    systematic_name_(),
    systematic_selection_(),
    systematic_weight_(),
    datacard_process_(),
    save_data_as_hist_(false) {

  if (!weight.IsScalar()) {
    throw std::invalid_argument(("Weight NamedFunc "+weight.Name()
                                 +" is not a scalar.").c_str());
  }
  if (!variable_.IsScalar()) {
    throw std::invalid_argument(("Signal extraction NamedFunc "+variable_.Name()
                                 +" is not a scalar.").c_str());
  }

  //initialize NamedFunc look-up vectors
  for (const SelectionList &channel_map : channels) {
    //nominal selection
    NamedFunc channel_selection(1);
    for (unsigned isel = 0; isel < channel_map.selection_.size(); isel++) {
      channel_selection = channel_selection && channel_map.selection_[isel];
    }
    channel_name_.push_back(channel_map.channel_name_);
    channel_selection_.push_back(channel_selection);
    systematic_selection_.push_back(std::vector<NamedFunc>());
    systematic_weight_.push_back(std::vector<NamedFunc>());
    for (const Systematic &systematic : systematics) {
      systematic_name_.push_back(systematic.name_);
      if (systematic.is_weight_systematic_) {
        systematic_selection_.back().push_back(channel_selection);
        systematic_weight_.back().push_back(systematic.content_);
      }
      else {
        NamedFunc systematic_selection(1);
        for (unsigned isel = 0; isel < channel_map.selection_.size(); isel++) {
          if (channel_map.name_[isel] == systematic.selection_name_) {
            systematic_selection = systematic_selection && systematic.content_;
          }
          else {
            systematic_selection = systematic_selection && channel_map.selection_[isel];
          }
        }
        systematic_selection_.back().push_back(systematic_selection);
        systematic_weight_.back().push_back(nominal_weight_);
      }
    }
  }

  //initialize processes
  for (const std::shared_ptr<Process> & process : processes) {
    datacard_process_.push_back(std::unique_ptr<DatacardProcess>(new DatacardProcess(*this, process, axis)));
  }
  
}

/*! \brief Produce and save datacard txt file and associate root files with 
           roofit workspace
  \param[in] luminosity   currently unused since lumi implemented via weight
  \param[in] subdir       subdirectory to save files to
*/
void Datacard::Print(double luminosity, const std::string &subdir) {

  //create subdirectory if needed
  std::string subdir_mod = subdir;
  if (subdir != "") {
    mkdir(("datacards/"+subdir).c_str(), 0777);
    subdir_mod = subdir + "/";
  }

  //dummy logic to avoid unused variable needed for interface
  luminosity += 0;

  //save RooFit workspaces in root file
  TFile root_file(("datacards/"+subdir_mod+name_+".root").c_str(),"RECREATE"); 
  std::vector<float> data_norm;
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    float lower_edge = datacard_process_[0]->raw_histogram_nom_[ichan]
        .GetXaxis()->GetBinLowEdge(1);
    float upper_edge = datacard_process_[0]->raw_histogram_nom_[ichan]
        .GetXaxis()->GetBinUpEdge(datacard_process_[0]->raw_histogram_nom_[ichan]
        .GetXaxis()->GetNbins());
    for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
      if (iproc < datacard_process_.size()) {
        std::string proc_name = datacard_process_[iproc]->process_->name_;
        std::string pdf_name = "pdf_"+proc_name+"_"+channel_name_[ichan];
        RooWorkspace ws(("WS_"+proc_name+"_"+channel_name_[ichan]).c_str());
        if (datacard_process_[iproc]->process_->type_ == Process::Type::data)  {
          proc_name = "data_obs";
          pdf_name = "pdf_"+proc_name+"_"+channel_name_[ichan];
          data_norm.push_back(
              datacard_process_[iproc]->raw_histogram_nom_[ichan].Integral());
          if (save_data_as_hist_) {
            RooRealVar variable((variable_.Name()+channel_name_[ichan]).c_str(),
                "",lower_edge,upper_edge);
            RooDataHist data_obs(pdf_name.c_str(),"",RooArgList(variable),
                &datacard_process_[iproc]->raw_histogram_nom_[ichan]);
            std::cout << "importing data_hist" << std::endl;
            ws.import(data_obs);
          }
          else {
            datacard_process_[iproc]->raw_dataset_nom_[ichan].SetName(
                pdf_name.c_str());
            std::cout << "importing dataset" << std::endl;
            //ws.import(datacard_process_[iproc]->weight_);
            ws.import(datacard_process_[iproc]->raw_dataset_nom_[ichan]);
          }
        }
        else {
          if (datacard_process_[iproc]->replace_with_param_) {
            std::cout << "importing replaced pdf" << std::endl;
            ws.import(*(datacard_process_[iproc]->param_pdf_[ichan]));
          }
          else {
            RooRealVar variable((variable_.Name()+channel_name_[ichan]).c_str(),
                "",lower_edge,upper_edge);
            RooDataHist hist((pdf_name+"_hist").c_str(),"",RooArgList(variable),
                &datacard_process_[iproc]->raw_histogram_nom_[ichan]);
            RooHistPdf pdf(pdf_name.c_str(),"",RooArgList(variable),hist);
            //float nom_yield = datacard_process_[iproc]->raw_histogram_nom_[ichan].Integral();
            //RooRealVar norm((pdf_name+"_norm").c_str(),"",nom_yield,
            //    0,100.0*nom_yield);
            //ws.import(norm);
            std::cout << "importing hist pdf" << std::endl;
            ws.import(pdf);
          }
        }
        ws.Write();
      }
      else {
        unsigned iproc_eff = iproc-datacard_process_.size();
        std::string proc_name = param_process_name_[iproc_eff];
        if (!param_process_profile_dec[iproc_eff]){
	  //	  std::cout << "Channel = " << ichan << ", Process = " <<proc_name << ", not profiled"<<std::endl;
	  std::string pdf_name = "pdf_"+proc_name+"_"+channel_name_[ichan] + "_" + param_func_name_[iproc_eff][ichan];
          RooRealVar norm((pdf_name+"_norm").c_str(),"",data_norm[ichan],
              0,3.0*data_norm[ichan]);
          RooWorkspace ws(("WS_"+proc_name+"_"+channel_name_[ichan]).c_str());
          ws.import(*param_process_[iproc_eff][ichan]);
          ws.import(norm);
          ws.Write();
        }
        else{
	  //	  std::cout << "Channel = " << ichan << ", Process = " << proc_name << ", profiles" << std::endl;
	  std::string profile_name = "profile_"+proc_name+"_"+channel_name_[ichan];
          RooRealVar norm((profile_name+"_norm").c_str(),"",data_norm[ichan],
			  0,3.0*data_norm[ichan]);
          RooWorkspace ws(("WS_"+proc_name+"_"+channel_name_[ichan]).c_str());
          RooCategory pdfindex(("pdfindex_" + proc_name +"_" + channel_name_[ichan]).c_str(), ("pdfindex_" + proc_name + "_"+channel_name_[ichan]).c_str());
	  auto models = RooArgList(("pdfmodels_" + proc_name +"_"+ channel_name_[ichan]).c_str());
	  for (auto pdf: param_profile_process_[iproc_eff][ichan]){
	    models.add(*pdf);
	  }
	  RooMultiPdf profile(("profile_" + proc_name + "_"+channel_name_[ichan]).c_str(), ("profile_" + proc_name + "_"+channel_name_[ichan]).c_str(), pdfindex, models);
	//	  ws.import(*(param_profile_ind_process_[iproc_eff][ichan]));
	  ws.import(pdfindex);
	  ws.import(profile);
          ws.import(norm);
          ws.Write();
        }
      }
    }
  }
  root_file.Close();
  std::cout << "open datacards/"+subdir_mod+name_+".root" << std::endl;

  //save datacard txt file
  std::ofstream datacard_file;
  datacard_file.open(("datacards/"+subdir_mod+name_+".txt").c_str(),std::ios::out);
  //header
  datacard_file << "max  " << n_channels_ << " number of categories\n";
  datacard_file << "jmax " << n_processes_-2 << " number of samples minus one\n";
  datacard_file << "kmax " << n_systematics_ <<" number of nuisance parameters\n";
  datacard_file << "-----------------------------------------------------------------"
                << "-----------------------------------------------------------------\n";
  //shape locations
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
      if (iproc < datacard_process_.size()) {
        std::string proc_name = datacard_process_[iproc]->process_->name_;
        if (datacard_process_[iproc]->process_->type_ == Process::Type::data)
          proc_name = "data_obs";
        datacard_file << "shapes " << std::left << std::setw(19)
            << proc_name << std::left << std::setw(19) << channel_name_[ichan] 
            << name_+".root " << "WS_"+proc_name+"_"+channel_name_[ichan] << ":"
            << "pdf_"+proc_name+"_"+channel_name_[ichan] << "\n";
      }
      else {
        unsigned iproc_eff = iproc-datacard_process_.size();
        if (!param_process_profile_dec[iproc_eff]){
            datacard_file << "shapes " << std::left << std::setw(19)
            << param_process_name_[iproc_eff] << std::left << std::setw(19) 
            << channel_name_[ichan] << name_+".root " 
            << "WS_"+param_process_name_[iproc_eff]+"_"+channel_name_[ichan] 
            << ":pdf_"+param_process_name_[iproc_eff]+"_"+channel_name_[ichan]+"_"+param_func_name_[iproc_eff][ichan]
            << "\n";
        }
        else {
          datacard_file << "shapes " << std::left << std::setw(19)
          << param_process_name_[iproc_eff] << std::left << std::setw(19) 
          << channel_name_[ichan] << name_+".root " 
          << "WS_"+param_process_name_[iproc_eff]+"_"+channel_name_[ichan] 
          << ":profile_"+param_process_name_[iproc_eff]+"_"+channel_name_[ichan] 
          << "\n";
        }
      }
    }
  }
  datacard_file << "-----------------------------------------------------------------"
                << "-----------------------------------------------------------------\n";
  //dummy observations
  datacard_file << std::left << std::setw(13) << "bin";
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    datacard_file << std::left << std::setw(15) << channel_name_[ichan];
  }
  datacard_file << "\n";
  datacard_file << std::left << std::setw(13) << "observation";
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    datacard_file << std::left << std::setw(15) << "-1";
  }
  datacard_file << "\n";
  datacard_file << "-----------------------------------------------------------------"
                << "-----------------------------------------------------------------\n";
  //process rates
  datacard_file << std::left << std::setw(33) << "bin";
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
      if (iproc < datacard_process_.size()) {
        if (datacard_process_[iproc]->process_->type_ == Process::Type::data)
          continue;
        datacard_file << std::left << std::setw(19) << channel_name_[ichan];
      }
      else {
        datacard_file << std::left << std::setw(19) << channel_name_[ichan];
      }
    }
  }
  datacard_file << "\n";
  datacard_file << std::left << std::setw(33) << "process";
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
      if (iproc < datacard_process_.size()) {
        if (datacard_process_[iproc]->process_->type_ == Process::Type::data)
          continue;
        std::string proc_name = datacard_process_[iproc]->process_->name_;
        datacard_file << std::left << std::setw(19) << proc_name;
      }
      else {
        unsigned iproc_eff = iproc-datacard_process_.size();
        datacard_file << std::left << std::setw(19) << param_process_name_[iproc_eff];
      }
    }
  }
  datacard_file << "\n";
  datacard_file << std::left << std::setw(33) << "process";
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    int signal_number = -1;
    int background_number = 1;
    for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
      if (iproc < datacard_process_.size()) {
        if (datacard_process_[iproc]->process_->type_ == Process::Type::data)
          continue;
        else if (datacard_process_[iproc]->process_->type_ == Process::Type::signal) {
          datacard_file << std::left << std::setw(19) << signal_number;
          signal_number -= 1;
        }
        else { //background
          datacard_file << std::left << std::setw(19) << background_number;
          background_number += 1;
        }
      }
      else { 
        //unsigned iproc_eff = iproc-datacard_process_.size();
        //if (param_process_type_[iproc_eff]==Process::Type::signal) {
        //  datacard_file << std::left << std::setw(19) << signal_number;
        //  signal_number -= 1;
        //}
        //else {
        datacard_file << std::left << std::setw(19) << background_number;
        background_number += 1;
        //}
      }
    }
  }
  datacard_file << "\n";
  datacard_file << std::left << std::setw(33) << "rate";
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
      if (iproc < datacard_process_.size()) {
        if (datacard_process_[iproc]->process_->type_ == Process::Type::data)
          continue;
        datacard_file << std::left << std::setw(19) << 
            datacard_process_[iproc]->raw_histogram_nom_[ichan].Integral();
      }
      else {
        datacard_file << std::left << std::setw(19) << "1";
      }
    }
  }
  datacard_file << "\n";
  datacard_file << "-----------------------------------------------------------------"
                << "-----------------------------------------------------------------\n";
  //systematics
  for (unsigned isyst = 0; isyst < n_systematics_; isyst++) {
    datacard_file << std::left << std::setw(25) << systematic_name_[isyst];
    datacard_file << std::left << std::setw(8) << "lnN";
    for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
      for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
        if (iproc < datacard_process_.size()) {
          if (datacard_process_[iproc]->process_->type_ == Process::Type::data)
            continue;
          float ratio = datacard_process_[iproc]->raw_histogram_nom_[ichan].Integral()/
              datacard_process_[iproc]->raw_histogram_sys_[ichan][isyst].Integral();
          if (ratio < 1) ratio = 1.0/ratio;

	  if (ratio == 1.) datacard_file << std::left << std::setw(19) << "-";//TODO make a way to make 1.0 systematics just be -
          else datacard_file << std::left << std::setw(19) << ratio;
        }
      }
    }
    datacard_file << "\n";
  }
  //no params yets

  //Discrete profile indices
  datacard_file << "\n";
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
      if (iproc >= datacard_process_.size()) {
        unsigned iproc_eff = iproc-datacard_process_.size();
        if (param_process_profile_dec[iproc_eff]) {
          datacard_file << std::left << std::setw(40) << "pdfindex_"+param_process_name_[iproc_eff]+"_"+channel_name_[ichan]; 
	  datacard_file << std::left << std::left << std::setw(19) << "discrete" << "\n";
        }
      }
    }
  }

  datacard_file.close();
  std::cout << "open datacards/"+subdir_mod+name_+".txt" << std::endl;
}

/*! \brief Add parametric process to datacard (ex. background constrained by fit)
  \param[in] name   Name of process
  \param[in] pdf    PDF for process
  If more than one PDF are found in the same channel, the discrete profile is turned on for the process. 
  And a RooMultiPdf object will be saved in the workspace.
*/
Datacard& Datacard::MakeProcessParametric(
    const std::string &name, 
    std::vector<std::shared_ptr<RooAbsPdf>> &pdf) {
  if (pdf.size() != n_channels_) 
    throw std::invalid_argument(("Wrong PDFs for parametric process "+name).c_str());

  bool found_proc = false;
  for (std::unique_ptr<DatacardProcess> &datacard_process : datacard_process_) {
    if (datacard_process->process_->name_ == name) {
      datacard_process->replace_with_param_ = true;
      datacard_process->param_pdf_ = pdf;
      found_proc = true;
    }
  }
  if (!found_proc)
    throw std::invalid_argument(("No process with name "+name).c_str());

  return *this;
}

/*! \brief Adds parametric processes
  \param[in] process  parametric process to add
*/
Datacard& Datacard::AddParametricProcess(
    const std::string &name, 
    std::vector<std::shared_ptr<RooAbsPdf>> &pdf) {
  //check for errors 
  if (pdf.size() != n_channels_)
    throw std::invalid_argument(
        ("Process "+name+" does not have enough channels.").c_str());
  //add
  datacard_process_parametric_.push_back(
        std::unique_ptr<DatacardProcessParametric>(
        new DatacardProcessParametric(name, pdf, *this)));
  datacard_process_.push_back(datacard_process_parametric_.back().get());
  n_processes_++;
  return *this;
}

/*! \brief Adds parametric processes
  \param[in] process  parametric process to add
*/
Datacard& Datacard::AddParametricProcess(
    const std::string &name, 
    std::vector<std::vector<std::shared_ptr<RooAbsPdf>>> &pdf) {
  //check for errors 
  if (pdf.size() != n_channels_)
    throw std::invalid_argument(
        ("Process "+name+" does not have enough channels.").c_str());
  //add
  datacard_process_parametric_.push_back(
        std::unique_ptr<DatacardProcessParametric>(
        new DatacardProcessParametric(name, pdf, *this)));
  datacard_process_.push_back(datacard_process_parametric_.back().get());
  n_processes_++;
  return *this;
}

/*! \brief Sets whether to save data as RooDataSet or RooDataHist
  \param[in] save_data_as_hist whether to save data as RooDataHist
*/
Datacard& Datacard::SaveDataAsHist(bool save_data_as_hist) {
  save_data_as_hist_ = save_data_as_hist;
  return *this;
}

/*! \brief Dummy logic since this method is just needed for interface
  \param[in] tag    -
*/
void Datacard::SetLuminosityTag(const std::string &tag) {
  //dummy logic to avoid unusued variable warning
  std::string temp_tag = tag;
  temp_tag += "";
}

/*! \brief returns the set of processes (not including parametric processes)
*/
std::set<const Process*> Datacard::GetProcesses() const {
  std::set<const Process*> processes;
  for (const std::unique_ptr<Datacard::DatacardProcess> &datacard_process : datacard_process_) {
    processes.insert(datacard_process->process_.get());
  }
  return processes;
}

/*! \brief Returns figure component associated with a particular process
  \param[in] process   process whose figure component to find
*/
Figure::FigureComponent * Datacard::GetComponent(const Process *process) {
  for (const std::unique_ptr<Datacard::DatacardProcess> &datacard_process : datacard_process_) {
    if (datacard_process->process_.get() == process){
      return datacard_process.get();
    }
  }
  throw std::invalid_argument(("Could not find process "+process->name_).c_str());
  return nullptr;
}
