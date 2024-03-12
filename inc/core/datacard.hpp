#ifndef H_DATACARD
#define H_DATACARD

#include <memory>
#include <set>
#include <string>
#include <vector>

#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "TH1.h"

#include "core/axis.hpp"
#include "core/baby.hpp"
#include "core/figure.hpp"
#include "core/named_func.hpp"
#include "core/process.hpp"
#include "core/RooMultiPdf.hpp"

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
    void WriteWorkspace(unsigned int channel) final;
    std::string WSName(unsigned int channel) final;
    std::string PDFName(unsigned int channel) final;
    float Yield(unsigned int channel, unsigned int systematic = 999) final;
    void FitAndFreeze(unsigned int channel);

    std::vector<TH1D> raw_histogram_nom_; //!< nominal histogram per channel 
    std::vector<std::vector<TH1D>> raw_histogram_sys_; //!< histogram per systematic per channel
    std::vector<RooDataSet> raw_dataset_nom_; //!< RooDataSet used for data
    std::vector<RooRealVar> var_;       //!< Signal extraction variable
    RooRealVar weight_;                 //!< Weight variable
    bool replace_with_param_;           //!< If process should be replaced by parametric model
    std::vector<std::shared_ptr<RooAbsPdf>> param_pdf_; //!< Parametric model pdf per channel

  private:
    DatacardProcessNonparametric() = delete;
    DatacardProcessNonparametric(const DatacardProcessNonparametric &) = delete;
    DatacardProcessNonparametric& operator=(const DatacardProcessNonparametric &) = delete;
    DatacardProcessNonparametric(DatacardProcessNonparametric &&) = delete;
    DatacardProcessNonparametric& operator=(DatacardProcessNonparametric &&) = delete;
  };

  /*!\brief Class representing a parametric (i.e. derived purely from a fit
   * to another data set) process in a datacard
  */
  class DatacardProcessParametric final: public DatacardProcess{
  public:
    DatacardProcessParametric(const std::string &name, 
        std::vector<std::shared_ptr<RooAbsPdf>> &pdf, const Figure& figure); 
    DatacardProcessParametric(const std::string &name, 
        std::vector<std::vector<std::shared_ptr<RooAbsPdf>>> &pdf, 
        const Figure& figure);
    ~DatacardProcessParametric() = default;
    void RecordEvent(const Baby &baby) final;
    void WriteWorkspace(unsigned int channel) final; 
    std::string WSName(unsigned int channel) final; 
    std::string PDFName(unsigned int channel) final; 
    float Yield(unsigned int channel, unsigned int systematic = 999) final;

    std::vector<std::vector<std::shared_ptr<RooAbsPdf>>> pdf_; //!<PDF for each channel

  private:
    static std::vector<std::vector<std::shared_ptr<RooAbsPdf>>> 
        MakeDoubleVector(std::vector<std::shared_ptr<RooAbsPdf>> pdfs);
    DatacardProcessParametric() = delete;
    DatacardProcessParametric(const DatacardProcessParametric &) = delete;
    DatacardProcessParametric& operator=(const DatacardProcessParametric &) = delete;
    DatacardProcessParametric(DatacardProcessParametric &&) = delete;
    DatacardProcessParametric& operator=(DatacardProcessParametric &&) = delete;
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
                                 std::vector<std::shared_ptr<RooAbsPdf>> &pdf);
  Datacard& AddParametricProcess(const std::string &name, 
                                 std::vector<std::vector<std::shared_ptr<RooAbsPdf>>> &pdf);
  Datacard& MakeProcessParametric(const std::string &name, 
                                  std::vector<std::shared_ptr<RooAbsPdf>> &pdf);
  Datacard& SaveDataAsHist(bool save_data_as_hist = true);
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
  std::string name_;                         //!<Datacard name
  unsigned int n_channels_;                  //!<Number of channels
  unsigned int n_processes_;                 //!<Number of processes
  unsigned int n_systematics_;               //!<Number of systematics
  NamedFunc nominal_weight_;                 //!<Nominal weight variable
  NamedFunc variable_;                       //!<Variable used for signal extraction
  std::vector<NamedFunc> channel_selection_; //!<Selection for each channel
  std::vector<std::string> channel_name_;    //!<Channel names
  std::vector<std::string> systematic_name_; //!<Name for each systematic
  std::vector<std::vector<NamedFunc>> systematic_selection_; //!<Selection variations for systematics
  std::vector<std::vector<NamedFunc>> systematic_weight_; //!<Weight variations for systematics
  std::vector<DatacardProcess*> datacard_process_; //!<All Processes in datacard
  std::vector<std::unique_ptr<DatacardProcessNonparametric>> datacard_process_nonparametric_; 
      //!<Nonparametric processes in datacard
  std::vector<std::unique_ptr<DatacardProcessParametric>> datacard_process_parametric_; 
      //!<Parametric processes in datacard
  bool save_data_as_hist_;                   //!<Flag indicating to save data as RDataSet or RDataHist

private:

};

#endif //H_DATACARD
