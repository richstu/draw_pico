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
#include "TFile.h"
#include "TH1.h"

#include "core/axis.hpp"
#include "core/baby.hpp"
#include "core/figure.hpp"
#include "core/named_func.hpp"
#include "core/process.hpp"

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
        std::string pdf_name = "pdf_"+proc_name+"_"+channel_name_[ichan];
        RooRealVar norm((pdf_name+"_norm").c_str(),"",data_norm[ichan],
            0,3.0*data_norm[ichan]);
        RooWorkspace ws(("WS_"+proc_name+"_"+channel_name_[ichan]).c_str());
        std::cout << "importing parametric pdf" << std::endl;
        ws.import(*param_process_[iproc_eff][ichan]);
        ws.import(norm);
        ws.Write();
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
        datacard_file << "shapes " << std::left << std::setw(19)
            << param_process_name_[iproc_eff] << std::left << std::setw(19) 
            << channel_name_[ichan] << name_+".root " 
            << "WS_"+param_process_name_[iproc_eff]+"_"+channel_name_[ichan] 
            << ":pdf_"+param_process_name_[iproc_eff]+"_"+channel_name_[ichan] 
            << "\n";
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
          datacard_file << std::left << std::setw(19) << ratio;
          //TODO make a way to make 1.0 systematics just be -
        }
        else {
          datacard_file << std::left << std::setw(19) << "-";
        }
      }
    }
    datacard_file << "\n";
  }
  //no params yets
  datacard_file.close();
  std::cout << "open datacards/"+subdir_mod+name_+".txt" << std::endl;
}

/*! \brief Add parametric process to datacard (ex. background constrained by fit)
  \param[in] name   Name of process
  \param[in] pdf    PDF for process
*/
Datacard& Datacard::AddParametricProcess(const std::string &name, 
                                         std::vector<RooAbsPdf*> &pdf) {
  if (pdf.size() != n_channels_) 
    throw std::invalid_argument(("Wrong number of PDFs for parametric process "+name).c_str());
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    if (pdf[ichan]->GetName() != "pdf_"+name+"_"+channel_name_[ichan]) {
      throw std::invalid_argument("PDF name error. Pleas use pdf_<process>_<channel>");
    }
  }

  n_processes_++;
  param_process_name_.push_back(name);
  //param_process_type_.push_back(type);
  param_process_.push_back(pdf);
  return *this;
}

/*! \brief Replaces an existing process with a parametric shape, while still 
           retaining its normalization
  \param[in] name   Name of process
  \param[in] pdf    PDF for process
*/
Datacard& Datacard::MakeProcessParametric(const std::string &name, 
                                          std::vector<RooAbsPdf*> &pdf) {
  if (pdf.size() != n_channels_) 
    throw std::invalid_argument(("Wrong PDFs for parametric process "+name).c_str());
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    if (pdf[ichan]->GetName() != "pdf_"+name+"_"+channel_name_[ichan]) {
      throw std::invalid_argument("PDF name error. Pleas use pdf_<process>_<channel>");
    }
  }

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
