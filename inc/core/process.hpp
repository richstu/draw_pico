#ifndef H_PROCESS
#define H_PROCESS

#include <memory>
#include <string>
#include <map>
#include <set>
#include <mutex>

#include "core/baby.hpp"
#include "core/named_func.hpp"
#include "core/utilities.hpp"

class Process : public TAttFill, public TAttLine, public TAttMarker{
public:
  enum class Type{data, background, signal};

  template<typename BabyType>
  static std::shared_ptr<Process> MakeShared(const std::string &name,
                                             Type type,
                                             int color,
                                             const std::set<std::string> &files,
                                             const NamedFunc &cut = true){
    return std::shared_ptr<Process>(new Process(static_cast<BabyType*>(nullptr),
                                                name, type, color, {files}, cut));
  }
  template<typename BabyType>
  static std::shared_ptr<Process> MakeShared(const std::string &name,
                                             Type type,
                                             int color,
                                             const std::vector<std::set<std::string>> &files,
                                             const NamedFunc &cut = true){
    return std::shared_ptr<Process>(new Process(static_cast<BabyType*>(nullptr),
                                                name, type, color, files, cut));
  }

  std::string name_;
  Type type_;
  NamedFunc cut_;
  int color_;

  std::set<Baby*> Babies() const;

  ~Process();

private:
  template<typename BabyType>
  Process(BabyType * dummy_baby,
          const std::string &name,
          Type type,
          int color,
          const std::vector<std::set<std::string>> &files,
          const NamedFunc &cut);

  Process() = delete;
  Process(const Process &) = delete;
  Process& operator=(const Process &) = delete;
  Process(Process &&) = delete;
  Process& operator=(Process &&) = delete;

  static std::set<std::unique_ptr<Baby> > baby_pool_;
  static std::mutex mutex_;
};

template<typename BabyType>
Process::Process(BabyType * /*dummy_baby*/,
                 const std::string &name,
                 Type type,
                 int color,
                 const std::vector<std::set<std::string>> &files,
                 const NamedFunc &cut):
  TAttFill(color, type == Type::background ? 1001 : 0),
  TAttLine(type == Type::background ? 1 : color,
           1,
           type == Type::background ? 1 :
           type == Type::signal ? 5
           : 3),
  TAttMarker(color, 20, 1.2),
  name_(name),
  type_(type),
  cut_(cut),
  color_(color){
  std::lock_guard<std::mutex> lock(mutex_);
  for (const auto &baby_files: files) {
    std::set<std::string> baby_full_files;
    for(const auto &file: baby_files){
      const auto &full_files = Glob(file);
      baby_full_files.insert(full_files.begin(), full_files.end());
    }
    bool found = false;
    for(auto &baby_p: baby_pool_){
      auto &baby = *baby_p;
      if(typeid(baby) != typeid(BabyType)) continue;
      if(baby_full_files == baby.FileNames()){
        baby.processes_.insert(this);
        found = true;
        break;
      }
    }
    if(!found){
      baby_pool_.emplace(static_cast<Baby*>(new BabyType(baby_full_files,
                                                         std::set<const Process*>{this})));
    }
  }
  }

#endif
