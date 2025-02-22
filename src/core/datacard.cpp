/*! \class Datacard
 
  \brief A datacard for use with Combine

  Datacard contains machinery to create a Higgs combine datacard

  TODO: add more description here as this class is developed
  TODO: add usage instructions for users
  TODO: add description for datacardprocess subclass

*/

/*! \class Datacard::DatacardProcess

  \brief Container for TH1Ds(->RooDataHist?) associated with a single Process
  
*/

/*! \class Datacard::Systematic

  \brief Container for systematics (nuisance parameters) 
  
  These come in several flavors implemented by various child classes
  two flavours: weight systematics where the nominal weight is 
  replaced by an alternative, and selection systematics where one selection
  is replaced by an alternative
*/

/*! \class Datacard::NormWeightSystematic

  \brief Container for systematics implemented as weight variations
  
  These come in two flavours: weight systematics where the nominal weight is 
  replaced by an alternative, and selection systematics where one selection
  is replaced by an alternative
*/

/*! \class Datacard::SelectionList

  \brief Container for selections to be applied to a particular category
*/

//TODO merge documentation with header

#include "core/datacard.hpp"

#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <sys/stat.h>

#include "RooAbsPdf.h"
#include "RooArgSet.h"
#include "RooArgList.h"
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

//----------------------------------------------------------------------------
//SelectionList
//----------------------------------------------------------------------------

/*!\brief Standard constructor
*/
Datacard::SelectionList::SelectionList(const std::string &channel_name) :
    channel_name_(channel_name) {}

/*!\brief Rename+copy constructor
*/
Datacard::SelectionList::SelectionList(const std::string& channel_name, 
                                       const SelectionList& selection_list) :
  channel_name_(channel_name),
  name_(selection_list.name_),
  selection_(selection_list.selection_) {}

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

//TODO implement sample comparison systematics somehow (ex. Tune)

/*!\brief Constructor for symmetric systematics
 
 \param[in] name             Name of systematic

 \param[in] selection_names  Name of selections to replace OR "weight" OR 
                             "fitvar"
 \param[in] variations       Replacement selections/weights
*/
Datacard::Systematic::Systematic(const std::string &name, 
    const std::vector<std::string> &selection_names, 
    const std::vector<NamedFunc> &variations) :
  is_symmetric(true),
  name_(name) {
  if (selection_names.size() != variations.size()) {
    throw std::invalid_argument("Selection names and replacements must be"
                                "the same size");
  }
  for (unsigned int ivar = 0; ivar < selection_names.size(); ivar++) {
    variation_[selection_names[ivar]] = std::make_shared<NamedFunc>(
        variations[ivar]);
  }
}

/*!\brief Constructor for asymmetric systematics
 
 \param[in] name             Name of systematic

 \param[in] selection_names  Name of selections to replace OR "weight" OR 
                             "fitvar"

 \param[in] variations_up    Replacement selections/weights

 \param[in] variations_dn    Replacement selections/weights
*/
Datacard::Systematic::Systematic(const std::string &name, 
    const std::vector<std::string> &selection_names, 
    const std::vector<NamedFunc> &variations_up,
    const std::vector<NamedFunc> &variations_dn) :
  is_symmetric(false),
  name_(name) {
  if ((selection_names.size() != variations_up.size()) 
      || (selection_names.size() != variations_dn.size())) {
    throw std::invalid_argument("Selection names and replacements must be"
                                "the same size");
  }
  for (unsigned int ivar = 0; ivar < selection_names.size(); ivar++) {
    variation_up_[selection_names[ivar]] = std::make_shared<NamedFunc>(
        variations_up[ivar]);
    variation_dn_[selection_names[ivar]] = std::make_shared<NamedFunc>(
        variations_dn[ivar]);
  }
}

//----------------------------------------------------------------------------
//DatacardProcess
//----------------------------------------------------------------------------

Datacard::DatacardProcess::DatacardProcess(const Figure &figure,
    const std::shared_ptr<Process> &process) :
    FigureComponent(figure, process) {}

Datacard::DatacardProcess::DatacardProcess(const Figure &figure) :
    FigureComponent(figure, std::shared_ptr<Process>(nullptr)) {}

//----------------------------------------------------------------------------
//DatacardProcessNonparametric
//----------------------------------------------------------------------------

/*!\brief Standard constructor
 
  \param[in] figure   Parent figure

  \param[in] process  Process used to fill histogram
*/
Datacard::DatacardProcessNonparametric::DatacardProcessNonparametric(
    const Figure &figure, const std::shared_ptr<Process> &process, 
    const Axis &axis, bool in_datacard) :
    DatacardProcess(figure, process),
    rrv_weight_(RooRealVar("weight","",-50.0,50.0)) {
  if (process_->type_ == Process::Type::data) {
    is_data_ = true;
    is_signal_ = false;
    process_->name_ = "data_obs";
  }
  else if (process_->type_ == Process::Type::signal) {
    is_data_ = false;
    is_signal_ = true;
  } 
  else {
    is_data_ = false;
    is_signal_ = false;
  }
  name_ = process_->name_;
  in_datacard_ = in_datacard;
  replace_with_param_ = false;
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  n_variations_ = 1;
  if (in_datacard_) {
    n_variations_ = datacard->n_variations_;
  }
  float axis_min = axis.Bins().at(0);
  float axis_max = axis.Bins().at(axis.Bins().size()-1);
  for (unsigned ichan = 0; ichan < datacard->n_channels_; ichan++) {
    is_profiled_.push_back(false);
    var_.push_back(RooRealVar(
        (axis.var_.Name()+"_"+datacard->channel_name_[ichan]).c_str(),"",
        axis_min,axis_max));
    var_[ichan].setBins(axis.Nbins());
    dataset_.push_back(std::vector<RooDataSet>());
    for (unsigned isyst = 0; isyst < n_variations_; isyst++) {
      dataset_[ichan].push_back(RooDataSet("","",RooArgSet(var_[ichan],
                                                           rrv_weight_),
                                           RooFit::WeightVar(rrv_weight_)));
    } // loop over variations
  } // loop over channels
}

void Datacard::DatacardProcessNonparametric::RecordEvent(const Baby &baby) {
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  for (unsigned ichan = 0; ichan < datacard->n_channels_; ichan++) {
    for (unsigned isyst = 0; isyst < n_variations_; isyst++) {
      float weight = datacard->weight_[isyst].GetScalar(baby);
      float fit_var = datacard->fit_var_[isyst].GetScalar(baby);
      if (datacard->channel_selection_[isyst][ichan].GetScalar(baby)) {
        var_[ichan].setVal(fit_var);
        rrv_weight_.setVal(weight);
        dataset_[ichan][isyst].add(RooArgSet(var_[ichan]),weight);
      }
    } //loop over systematics/variations
  } //loop over channels
}

/*! \brief Returns the name of the RooWorkspace saved for a given channel
 
  \param[in] channel    index of channel to write
*/
std::string Datacard::DatacardProcessNonparametric::WSName(
    unsigned int channel) {
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  if (is_data_) {
    return "WS_data_obs_"+datacard->channel_name_[channel];
  }
  //is MC
  return "WS_"+name_+"_"+datacard->channel_name_[channel];
}

/*! \brief Returns the name of the RooDataSet/Hist saved for a given channel
 
  \param[in] channel    index of channel to write
*/
std::string Datacard::DatacardProcessNonparametric::DataName(
    unsigned int channel, unsigned int variation) {
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  if (is_data_) {
    return "data_obs_"+datacard->channel_name_[channel];
  }
  return ("mcdata_"+name_+"_"+datacard->channel_name_[channel]+"_"
          +datacard->variation_name_[variation]);
}


/*! \brief Returns the name of the PDF for a given channel
 
  \param[in] channel    index of channel to write
*/
std::string Datacard::DatacardProcessNonparametric::PDFName(
    unsigned int channel) {
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  if (is_data_) {
    return "data_obs_"+datacard->channel_name_[channel];
  }
  return (name_+"_"+datacard->channel_name_[channel]);
}

/*! \brief Returns yield (rate) of events in a given channel
 
  \param[in] channel    index of channel
*/
float Datacard::DatacardProcessNonparametric::Yield(unsigned int channel, 
                                                    unsigned int variation) {
  return dataset_[channel][variation].sumEntries();
}

/*! \brief Writes datasets to a workspace and saves to opened ROOT file
 
  \param[in] channel    index of channel to write
*/
void Datacard::DatacardProcessNonparametric::WriteWorkspace(
    unsigned int channel) {
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  RooWorkspace ws(WSName(channel).c_str());
  for (unsigned int ivar = 0; ivar < n_variations_; ivar++) {
    if (datacard->save_data_as_hist_) {
      //TODO check if weights work correctly for this
      RooDataHist binned_data(DataName(channel,ivar).c_str(),"",
                              RooArgSet(var_[channel]),
                              dataset_[channel][ivar]);
      ws.import(binned_data);
    }
    else {
      dataset_[channel][ivar].SetName(DataName(channel,ivar).c_str());
      ws.import(dataset_[channel][ivar]);
    }
  }
  ws.Write();
}

//----------------------------------------------------------------------------
//DatacardProcessParametric
//----------------------------------------------------------------------------

/*! \brief Add parametric process to datacard (background constrained by fit)
 
  \param[in] name   Name of process
*/
Datacard::DatacardProcessParametric::DatacardProcessParametric(
    const std::string &name, const Figure &figure) :
    DatacardProcess(figure) {

  //initialize properties
  name_ = name;
  is_data_ = false;
  is_signal_ = false; 
  in_datacard_ = true;
}

/*! \brief Wraps a vector of PDFs as a vector of single-entry vectors of PDFs
 
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
std::string Datacard::DatacardProcessParametric::PDFName(
    unsigned int channel) {
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  return "pdf_"+name_+"_"+datacard->channel_name_[channel];
}

/*! \brief Returns the name of the RooWorkspace saved for a given channel
 
  \param[in] channel    index of channel
*/
std::string Datacard::DatacardProcessParametric::WSName(unsigned int channel) {
  //if (figure_==nullptr)
  //  throw std::runtime_error(("WSName called for "+name_
  //                           +" before datacard assignment.").c_str());
  const Datacard* datacard = static_cast<const Datacard*>(&figure_);
  //is MC
  return "WS_"+name_+"_"+datacard->channel_name_[channel];
}

/*! \brief Returns a yield (rate) of 1.0 per combine conventions
 
*/
float Datacard::DatacardProcessParametric::Yield(unsigned int channel, 
                                                 unsigned int variation) {
  channel += 0;
  variation += 0;
  return 1.0;
}

/*! \brief Dummy function for interface
 
  \param[in] channel    index of channel to write
*/
void Datacard::DatacardProcessParametric::WriteWorkspace(
    unsigned int channel) {
  channel += 0;
}

//----------------------------------------------------------------------------
//Datacard
//----------------------------------------------------------------------------

/*! \brief Constructor for shape-based datacard

  \param[in] name         datacard name

  \param[in] channels     selection list for each channel

  \param[in] systematics  list of systematics

  \param[in] processes    standard processes including data, signal, and 
                          MC-based backgrounds

  \param[in] weight       nominal weight

  \param[in] axis         associated histogram axis
*/
//TODO allow different axes in different categories?
Datacard::Datacard(const std::string &name,
                   const std::vector<SelectionList> &channels, 
                   const std::vector<Systematic> &systematics,
                   const std::vector<std::shared_ptr<Process>> &processes,
                   const NamedFunc &weight,
                   const Axis &axis,
                   bool save_as_hist) :
    name_(name),
    n_channels_(channels.size()),
    n_processes_(processes.size()),
    n_systematics_(systematics.size()),
    axis_(axis),
    save_data_as_hist_(save_as_hist) {

  if (!weight.IsScalar()) {
    throw std::invalid_argument(("Weight NamedFunc "+weight.Name()
                                 +" is not a scalar.").c_str());
  }
  if (!axis.var_.IsScalar()) {
    throw std::invalid_argument(("Signal extraction NamedFunc "
                                 +axis.var_.Name()
                                 +" is not a scalar.").c_str());
  }

  for (const SelectionList &channel_map : channels) {
    channel_name_.push_back(channel_map.channel_name_);
  }

  //initialize NamedFunc look-up vectors
  systematics_extended_ = systematics;
  systematics_extended_.insert(systematics_extended_.begin(), 
      Systematic("nominal",{},{}));
  unsigned int ivariation = 0;
  n_variations_ = 0;
  for (const Systematic &systematic : systematics_extended_) {
    std::vector<const std::unordered_map<std::string, 
        std::shared_ptr<NamedFunc>>*> variations;
    if (systematic.is_symmetric) {
      variations.push_back(&systematic.variation_);
      variation_name_.push_back(systematic.name_);
    }
    else {
      variations.push_back(&systematic.variation_up_);
      variations.push_back(&systematic.variation_dn_);
      variation_name_.push_back(systematic.name_+"Up");
      variation_name_.push_back(systematic.name_+"Down");
    }
    for (const std::unordered_map<std::string, std::shared_ptr<NamedFunc>>* 
        variation : variations) {

      channel_selection_.push_back(std::vector<NamedFunc>());
      for (const SelectionList &channel_map : channels) {
        NamedFunc channel_selection(1);
        for (unsigned isel = 0; isel < channel_map.selection_.size(); isel++) {
          std::string selection_name = channel_map.name_[isel];
          if (variation->count(selection_name) != 0) {
            channel_selection = channel_selection 
                                && *(variation->at(selection_name));
          }
          else {
            channel_selection = channel_selection 
                                && channel_map.selection_[isel];
          }
        } //loop over selections in list
        channel_selection_[ivariation].push_back(channel_selection);
      } //loop over channels (selection lists)
        
      if (variation->count("weight") != 0) {
        weight_.push_back(*(variation->at("weight")));
      }
      else {
        weight_.push_back(weight);
      }
      if (variation->count("fitvar") != 0) {
        fit_var_.push_back(*(variation->at("fitvar")));
      }
      else {
        fit_var_.push_back(axis.var_);
      }
        
      ivariation++;
      n_variations_++;
    } //loop over variations (systematics)
  } //loop over systematics
  
  //initialize processes
  int n_data = 0;
  for (const std::shared_ptr<Process> & process : processes) {
    datacard_process_nonparametric_.push_back(
        std::unique_ptr<DatacardProcessNonparametric>(
        new DatacardProcessNonparametric(*this, process, axis)));
    datacard_process_.push_back(datacard_process_nonparametric_.back().get());
    if (datacard_process_.back()->is_data_)
      n_data++;
  }
  if (n_data != 1)
    throw std::invalid_argument("Exactly 1 data process should be provided"
                                "to datacard.");
}

/*! \brief Add processes for histogram generation but not datacard
 
  \param[in] processes  processes to add
 */
Datacard& Datacard::AddHistOnlyProcesses(
    const std::vector<std::shared_ptr<Process>> &processes) {
  //initialize processes
  for (const std::shared_ptr<Process> & process : processes) {
    datacard_process_nonparametric_.push_back(
        std::unique_ptr<DatacardProcessNonparametric>(
        new DatacardProcessNonparametric(*this, process, axis_, false)));
    datacard_process_.push_back(datacard_process_nonparametric_.back().get());
    if (datacard_process_.back()->is_data_)
      throw std::invalid_argument("Must provide data as datacard process.");
  }
  n_processes_ += processes.size();
  return *this;
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
  //These should contain all histograms needed by HtoZG_fitting to generate
  //models
  TFile root_file(("datacards/"+subdir_mod+name_+"_rawdata.root").c_str(),
                  "RECREATE"); 
  //loop over channels
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    //write the workspace for each process/channel to the datacard
    for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
      datacard_process_[iproc]->WriteWorkspace(ichan); 
    }
  }
  root_file.Close();
  std::cout << "open datacards/"+subdir_mod+name_+"_rawdata.root" << std::endl;

  //save datacard txt file
  std::ofstream datacard_file;
  datacard_file.open(("datacards/"+subdir_mod+name_+".txt").c_str(),
                     std::ios::out);
  //header
  datacard_file << "max  " << n_channels_ << " number of categories\n";
  //subtract 1 for data, which is not counted, and 1 for combine conventions
  datacard_file << "jmax " << n_processes_-2 
                << " number of samples minus one\n";
  datacard_file << "kmax " << n_systematics_ 
                << " number of nuisance parameters\n";
  datacard_file << "----------------------------------------------------------"
                   "----------------------------------------------------------"
                   "--------------\n";
  //TODO adjust padding to fit strings
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
      if (datacard_process_[iproc]->in_datacard_) {
        std::string proc_name = datacard_process_[iproc]->name_;

        datacard_file << "shapes " << std::left << std::setw(19)
            << proc_name << std::left << std::setw(19) << channel_name_[ichan] 
            << name_+".root " << datacard_process_[iproc]->WSName(ichan) << ":"
            << datacard_process_[iproc]->PDFName(ichan) << "\n";
      } //is data or parametric
    } //loop over processes
  } //loop over channels
  datacard_file << "----------------------------------------------------------"
                   "----------------------------------------------------------"
                   "--------------\n";
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
  datacard_file << "----------------------------------------------------------"
                   "----------------------------------------------------------"
                   "--------------\n";
  //process rates
  std::ostringstream bin_str, proc_str, index_str, rate_str;
  for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
    int signal_number = -1;
    int background_number = 1;
    for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
      if (datacard_process_[iproc]->in_datacard_ 
          && !datacard_process_[iproc]->is_data_) {
        bin_str << std::left << std::setw(19) << channel_name_[ichan];
        proc_str << std::left << std::setw(19) 
                 << datacard_process_[iproc]->name_;
        if (datacard_process_[iproc]->is_signal_) {
          index_str << std::left << std::setw(19) << signal_number;
          signal_number -= 1;
        }
        else { //background
          index_str << std::left << std::setw(19) << background_number;
          background_number += 1;
        }
        rate_str << std::left << std::setw(19) 
                 << datacard_process_[iproc]->Yield(ichan);
      }
    }
  }
  datacard_file << std::left << std::setw(33) << "bin" << bin_str.str() 
                << "\n";
  datacard_file << std::left << std::setw(33) << "process" << proc_str.str() 
                << "\n";
  datacard_file << std::left << std::setw(33) << "process" << index_str.str() 
                << "\n";
  datacard_file << std::left << std::setw(33) << "rate" << rate_str.str() 
                << "\n";
  datacard_file << "----------------------------------------------------------"
                   "----------------------------------------------------------"
                   "--------------\n";
  //systematics
  unsigned int ivar = 0;
  for (Systematic& systematic : systematics_extended_) {
    if (ivar == 0) {
      ivar++;
      continue;
    }
    if ((systematic.variation_.count("FitVar") == 0)
        && (systematic.variation_up_.count("FitVar") == 0)) {
      //only make lnN constraints for non-shape systematics
      datacard_file << std::left << std::setw(25) << systematic.name_;
      datacard_file << std::left << std::setw(8) << "lnN";
      for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
        for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
          if (datacard_process_[iproc]->in_datacard_ 
              && !datacard_process_[iproc]->is_data_) {
            if (systematic.is_symmetric) {
              float nom_yield = datacard_process_[iproc]->Yield(ichan, 0);
              float syst = 1.0;
              if (nom_yield > 0.0)
                syst = (datacard_process_[iproc]->Yield(ichan, ivar)
                        /nom_yield);
              //systematics less than 0.1% are dropped
              if (fabs(syst-1.0)<1.0e-3)
                datacard_file << std::left << std::setw(19) << "-";
              else
                datacard_file << std::left << std::setw(19) << syst;
            }
            else {
              float nom_yield = datacard_process_[iproc]->Yield(ichan, 0);
              float syst_up = 1.0;
              float syst_dn = 1.0;
              if (nom_yield > 0.0) {
                syst_up = (datacard_process_[iproc]->Yield(ichan, ivar)
                           /nom_yield);
                syst_dn = (datacard_process_[iproc]->Yield(ichan, ivar+1)
                           /nom_yield);
              }
              //systematics less than 0.1% are dropped
              if ((fabs(syst_up-1.0)<1.0e-3) && (fabs(syst_dn-1.0)<1.0e-3))
                datacard_file << std::left << std::setw(19) << "-";
              else {
                std::ostringstream syst_string;
                syst_string << syst_dn << "/" << syst_up;
                datacard_file << std::left << std::setw(19) 
                              << syst_string.str();
              }
            } //asymmetric systematic
          } //process is in syst section
        } //loop over processes
      } //loop over channels
      datacard_file << "\n";
    } //not shape systematic

    if (systematic.is_symmetric) {
      ivar++;
    }
    else {
      ivar += 2;
    }
  }

  //Parameters: deal with this in HtoZG_fitting
  //datacard_file << "\n";
  //for (unsigned ichan = 0; ichan < n_channels_; ichan++) {
  //  for (unsigned iproc = 0; iproc < n_processes_; iproc++) {
  //    if (datacard_process_[iproc]->is_profiled_[ichan]) {
  //      datacard_file << std::left << std::setw(40) 
  //                    << ("pdfindex_"+datacard_process_[iproc]->name_+"_"
  //                        +channel_name_[ichan])
  //                    << std::left << std::left << std::setw(19) 
  //                    << "discrete" << "\n";
  //    }
  //  }
  //}

  datacard_file.close();
  std::cout << "open datacards/"+subdir_mod+name_+".txt" << std::endl;
}

/*! \brief Adds parametric processes
 
  \param[in] process  parametric process to add
*/
Datacard& Datacard::AddParametricProcess(
    const std::string &name) {
  //add
  datacard_process_parametric_.push_back(
        std::unique_ptr<DatacardProcessParametric>(
        new DatacardProcessParametric(name, *this)));
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
  for (const std::unique_ptr<Datacard::DatacardProcessNonparametric> 
      &datacard_process : datacard_process_nonparametric_) {
    processes.insert(datacard_process->process_.get());
  }
  return processes;
}

/*! \brief Returns figure component associated with a particular process
 
  \param[in] process   process whose figure component to find
*/
Figure::FigureComponent * Datacard::GetComponent(const Process *process) {
  for (const std::unique_ptr<Datacard::DatacardProcessNonparametric> 
      &datacard_process : datacard_process_nonparametric_) {
    if (datacard_process->process_.get() == process){
      return datacard_process.get();
    }
  }
  throw std::invalid_argument(("Could not find process "
                               +process->name_).c_str());
  return nullptr;
}
