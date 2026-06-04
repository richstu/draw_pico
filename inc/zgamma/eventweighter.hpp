//draw_pico version of nano2pico event weight for on-the-fly corrections

#ifndef H_EVENTWEIGHTER
#define H_EVENTWEIGHTER

#include "core/baby.hpp"
#include "core/correction.hpp"
#include "core/named_func.hpp"

#include <string>
#include <memory>
#include <vector>

class ElectronEventWeighter{
public:
  ElectronEventWeighter(std::string year);

  std::vector<double> ElectronSF(const Baby &b);

  std::vector<double> ElectronMinisoSF(const Baby &b);

  std::vector<double> ElectronCombinedSF(const Baby &b);

private:
  std::string in_file_electron_;
  std::string in_file_electron_iso0p10_;
  std::string in_file_electron_iso0p15_;
  std::string in_file_electron_reco_;
  std::unique_ptr<correction::CorrectionSet> cs_electron_;
  std::unique_ptr<correction::CorrectionSet> cs_electron_reco_;
  std::unique_ptr<correction::CorrectionSet> cs_electron_bpixhole_;
  std::unique_ptr<correction::CorrectionSet> cs_el_iso0p10_;
  std::unique_ptr<correction::CorrectionSet> cs_el_iso0p15_;
  std::unique_ptr<correction::CorrectionSet> cs_el_hole_iso0p10_;
  std::unique_ptr<correction::CorrectionSet> cs_el_hole_iso0p15_;
  correction::Correction::Ref map_electron_id_pass_;
  correction::Correction::Ref map_electron_id_pass_unc_;
  correction::Correction::Ref map_electron_id_fail_;
  correction::Correction::Ref map_electron_id_fail_unc_;
  correction::Correction::Ref map_electron_hole_id_pass_;
  correction::Correction::Ref map_electron_hole_id_pass_unc_;
  correction::Correction::Ref map_electron_hole_id_fail_;
  correction::Correction::Ref map_electron_hole_id_fail_unc_;
  correction::Correction::Ref map_electron_reco_pass_;
  correction::Correction::Ref map_electron_reco_pass_unc_;
  correction::Correction::Ref map_electron_reco_fail_;
  correction::Correction::Ref map_electron_reco_fail_unc_;
  bool post_bpix_;
};

extern const NamedFunc sys_el;
extern const NamedFunc w_el;
extern const NamedFunc sys_el_up;
extern const NamedFunc sys_el_dn;

float get_w_fakephoton(std::string run, int ph_isjet, float ph_pt, float ph_eta);

#endif
