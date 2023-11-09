#ifndef H_DATACARD
#define H_DATACARD

#include <memory>
#include <set>
#include <string>
#include <vector>

#include "RooAbsPdf.h"
#include "RooMultiPdf.h"
#include "TH1.h"

#include "core/axis.hpp"
#include "core/baby.hpp"
#include "core/figure.hpp"
#include "core/named_func.hpp"
#include "core/process.hpp"

class Datacard final: public Figure{
public:

  class SelectionList {
  public:
    SelectionList(const std::string& name);
    SelectionList() = default;
    SelectionList(const SelectionList &) = default;
    SelectionList& operator=(const SelectionList &) = default;
    SelectionList(SelectionList &&) = default;
    SelectionList& operator=(SelectionList &&) = default;
    SelectionList& AddSelection(const std::string &name, const NamedFunc &selection);
    std::string channel_name_;
    std::vector<std::string> name_;
    std::vector<NamedFunc> selection_;
  private:
  };

  class Systematic {
  public:
    Systematic(const std::string &name, const NamedFunc &alternate_weight);
    Systematic(const std::string &name, const std::string &selection_name, 
               const NamedFunc &alternate_selection);
    Systematic() = default;
    Systematic(const Systematic &) = default;
    Systematic& operator=(const Systematic &) = default;
    Systematic(Systematic &&) = default;
    Systematic& operator=(Systematic &&) = default;
    //for now, hardcoded to only affect signal, later allow systematics
    //affecting specific processes
    bool is_weight_systematic_;
    std::string name_;
    std::string selection_name_;
    NamedFunc content_;
  private:
  };
  //TODO implement asymmetric systematics

  class DatacardProcess final: public Figure::FigureComponent{
  public:
    DatacardProcess(const Figure &figure,
                    const std::shared_ptr<Process> &process,
                    const Axis &axis);
    ~DatacardProcess() = default;
    void RecordEvent(const Baby &baby) final;

    //1 nominal histogram per channel
    std::vector<TH1D> raw_histogram_nom_; 
    //1 histogram per channel per systematic
    std::vector<std::vector<TH1D>> raw_histogram_sys_; 
    bool is_data_;

  private:
    DatacardProcess() = delete;
    DatacardProcess(const DatacardProcess &) = delete;
    DatacardProcess& operator=(const DatacardProcess &) = delete;
    DatacardProcess(DatacardProcess &&) = delete;
    DatacardProcess& operator=(DatacardProcess &&) = delete;
  };

  //cut-and-count constructor, to be implemented
  //Datacard(std::vector<std::unordered_map<NamedFunc>> &channels, 
  //         std::vector<Systematic> &systematics,
  //         std::vector<std::shared_ptr<Process>> &processes);
  Datacard(const std::string &name,
           const std::vector<SelectionList> &channels, 
           const std::vector<Systematic> &systematics,
           const std::vector<std::shared_ptr<Process>> &processes,
           const NamedFunc &weight,
           const Axis &axis);
  Datacard& AddParametricProcess(const std::string &name, 
                                 std::vector<RooAbsPdf*> &pdf); //this will likely not suffice for more complicated procedures like discrete profiling
  Datacard(Datacard &&) = default;
  Datacard& operator=(Datacard &&) = default;
  ~Datacard() = default;

  //functions overwriting parent virtual functions
  void Print(double luminosity,
             const std::string &subdir) final;
  void SetLuminosityTag(const std::string &tag) final;
  std::set<const Process*> GetProcesses() const final;
  std::string GetTag() const final {return "";}
  FigureComponent * GetComponent(const Process *process) final;

  //member data
  std::string name_;
  unsigned int n_channels_;
  unsigned int n_processes_;
  unsigned int n_systematics_;
  NamedFunc nominal_weight_; 
  NamedFunc variable_;
  std::vector<NamedFunc> channel_selection_;
  std::vector<std::string> channel_name_;
  std::vector<std::string> systematic_name_;
  std::vector<std::vector<NamedFunc>> systematic_selection_;
  std::vector<std::vector<NamedFunc>> systematic_weight_;
  std::vector<std::unique_ptr<DatacardProcess>> datacard_process_;
  std::vector<std::string> param_process_name_;
  std::vector<std::vector<std::string>> param_func_name_;
  std::vector<bool> param_process_profile_dec;
  std::vector<std::vector<RooAbsPdf*>> param_process_;
  std::vector<std::vector<RooMultiPdf>> param_profile_process_;
  std::vector<std::vector<RooCategory>> param_profile_ind_process_;


private:

};

#endif //H_DATACARD
