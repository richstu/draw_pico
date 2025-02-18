#ifndef H_DATACARD
#define H_DATACARD

#include <memory>
#include <set>
#include <string>
#include <vector>

#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "TH1.h"

#include "core/axis.hpp"
#include "core/baby.hpp"
#include "core/figure.hpp"
#include "core/named_func.hpp"
#include "core/process.hpp"

/*!\brief Figure subclass that produces combine datacards (statistical models)
*/
class Datacard final: public Figure{
public:

  /*!\brief Class used to store selections applied to sample 
  */
  class SelectionList {
  public:
    SelectionList(const std::string& name);
    SelectionList(const std::string& name, 
                  const SelectionList& selection_list);
    SelectionList() = default;
    SelectionList(const SelectionList &) = default;
    SelectionList& operator=(const SelectionList &) = default;
    SelectionList(SelectionList &&) = default;
    SelectionList& operator=(SelectionList &&) = default;
    SelectionList& AddSelection(const std::string &name, 
                                const NamedFunc &selection);
    std::string channel_name_;
    std::vector<std::string> name_;
    std::vector<NamedFunc> selection_;
  private:
  };

  /*!\brief Class that stores systematic variations, in the selections, the
   * weight applied to each sample, or the fit variable
  */
  class Systematic {
  public:
    Systematic(const std::string &name, 
               const std::vector<std::string> &selection_names, 
               const std::vector<NamedFunc> &variations);

    Systematic(const std::string &name, 
               const std::vector<std::string> &selection_names, 
               const std::vector<NamedFunc> &variations_up,
               const std::vector<NamedFunc> &variations_dn);

    Systematic() = default;
    Systematic(const Systematic &) = default;
    Systematic& operator=(const Systematic &) = default;
    Systematic(Systematic &&) = default;
    Systematic& operator=(Systematic &&) = default;

    bool is_symmetric; //!< systematic type
    std::string name_; //!< systematic name
                       
    std::unordered_map<std::string, std::shared_ptr<NamedFunc>> variation_; 
        //!< map from selection to replace (or "weight") to alternate NamedFunc
    std::unordered_map<std::string, std::shared_ptr<NamedFunc>> variation_up_;
        //!< map from selection to replace (or "weight") to alt up NamedFunc
    std::unordered_map<std::string, std::shared_ptr<NamedFunc>> variation_dn_;
        //!< map from selection to replace (or "weight") to alt dn NamedFunc
  };

  /*!\brief Interface representing processes appearing in the datacard, see 
   * below DatacardProcessNonparametric and DatacardProcessParametric classes
  */
  class DatacardProcess : public Figure::FigureComponent{
  public:
    DatacardProcess(const Figure &figure, 
                    const std::shared_ptr<Process> &process);
    DatacardProcess(const Figure &figure);
    virtual void WriteWorkspace(unsigned int channel) = 0; 
    virtual std::string WSName(unsigned int channel) = 0;
    virtual std::string PDFName(unsigned int channel) = 0;
    virtual float Yield(unsigned int channel, 
                        unsigned int variation = 0) = 0;

    std::string name_;              //!< Process name
    bool is_data_;                  //!< If process is data
    bool is_signal_;                //!< If process is signal
    bool in_datacard_;              //!< If process is to be included in model
    std::vector<bool> is_profiled_; //!< If each channel is profiled
    std::vector<float> data_norm_;  //!< Number of data events for each channel
  };

  /*!\brief Class representing a nonparametric (i.e. derived at some level
   * from data or MC samples) process in a datacard
  */
  class DatacardProcessNonparametric final: public DatacardProcess{
  public:
    DatacardProcessNonparametric(const Figure &figure,
                                 const std::shared_ptr<Process> &process,
                                 const Axis &axis, bool in_datacard=true);
    ~DatacardProcessNonparametric() = default;
    void RecordEvent(const Baby &baby) final;
    void WriteWorkspace(unsigned int channel) final;
    std::string WSName(unsigned int channel) final;
    std::string DataName(unsigned int channel, unsigned int variation);
    std::string PDFName(unsigned int channel) final;
    float Yield(unsigned int channel, unsigned int variation = 0) final;

    unsigned int n_variations_; //!< Number of datasets to store per channel

    std::vector<std::vector<RooDataSet>> dataset_;
        //!< data indexed by channel, then variation
    std::vector<RooRealVar> var_; //!< Signal extraction variable

    RooRealVar rrv_weight_;       //!< Weight variable
    bool replace_with_param_;     //!< Processis replaced by parametric model
    std::vector<std::shared_ptr<RooAbsPdf>> param_pdf_; //!< Parametric model
                                                        //!< by channel

  private:
    DatacardProcessNonparametric() = delete;
    DatacardProcessNonparametric(const DatacardProcessNonparametric &) 
        = delete;
    DatacardProcessNonparametric& operator=(
        const DatacardProcessNonparametric &) = delete;
    DatacardProcessNonparametric(DatacardProcessNonparametric &&) = delete;
    DatacardProcessNonparametric& operator=(DatacardProcessNonparametric &&) 
        = delete;
  };

  /*!\brief Class representing a parametric (i.e. derived purely from a fit
   * to another data set) process in a datacard
  */
  class DatacardProcessParametric final: public DatacardProcess{
  public:
    DatacardProcessParametric(const std::string &name, const Figure& figure); 
    ~DatacardProcessParametric() = default;
    void RecordEvent(const Baby &baby) final;
    void WriteWorkspace(unsigned int channel) final; 
    std::string WSName(unsigned int channel) final; 
    std::string PDFName(unsigned int channel) final; 
    float Yield(unsigned int channel, unsigned int variation = 0) final;

    std::vector<std::vector<std::shared_ptr<RooAbsPdf>>> pdf_; //!<channel PDF

  private:
    static std::vector<std::vector<std::shared_ptr<RooAbsPdf>>> 
        MakeDoubleVector(std::vector<std::shared_ptr<RooAbsPdf>> pdfs);
    DatacardProcessParametric() = delete;
    DatacardProcessParametric(const DatacardProcessParametric &) = delete;
    DatacardProcessParametric& operator=(const DatacardProcessParametric &) 
        = delete;
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
           const Axis &axis,
           bool save_as_hist = false);
  Datacard& AddHistOnlyProcesses(
      const std::vector<std::shared_ptr<Process>> &processes);
  Datacard& AddParametricProcess(const std::string &name);
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
  unsigned int n_variations_;                //!<Internal number of systematics
  Axis axis_; //!< Axis
                                             
  std::vector<std::vector<NamedFunc>> channel_selection_;
      //!< Selections indexed by systematic, then channel
  std::vector<NamedFunc> fit_var_;           //!< Fitting var by systematic
  std::vector<NamedFunc> weight_;            //!< Weight indexed by systematic
                                             
  std::vector<Systematic> systematics_extended_; //!< systematics + nominal
  std::vector<std::string> channel_name_;    //!<Channel names
  std::vector<std::string> variation_name_;  //!<Name for each variation
  std::vector<DatacardProcess*> datacard_process_; //!<All Processes
  std::vector<std::unique_ptr<DatacardProcessNonparametric>> 
      datacard_process_nonparametric_; //!<Nonparametric processes in datacard
  std::vector<std::unique_ptr<DatacardProcessParametric>> 
      datacard_process_parametric_; //!<Parametric processes in datacard
  bool save_data_as_hist_;                   //!<Save data as RooDataHist

private:

};

#endif //H_DATACARD
