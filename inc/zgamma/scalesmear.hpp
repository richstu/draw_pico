// Joseph's code implementing on-the-fly energy error scale and smear for run 3

#ifndef H_SCALESMEAR
#define H_SCALESMEAR

#include "core/named_func.hpp"
#include "core/correction.hpp"

#include <memory>
#include <string>

extern std::unique_ptr<correction::CorrectionSet> cs_scale_syst_2022;
extern std::unique_ptr<correction::CorrectionSet> cs_scale_syst_2022EE;
extern std::unique_ptr<correction::CorrectionSet> cs_scale_syst_2023;
extern std::unique_ptr<correction::CorrectionSet> cs_scale_syst_2023BPix;

extern correction::CompoundCorrection::Ref map_scale_2022;
extern correction::Correction::Ref map_smearing_2022;
extern correction::CompoundCorrection::Ref map_scale_2022EE;
extern correction::Correction::Ref map_smearing_2022EE;
extern correction::CompoundCorrection::Ref map_scale_2023;
extern correction::Correction::Ref map_smearing_2023;
extern correction::CompoundCorrection::Ref map_scale_2023BPix;
extern correction::Correction::Ref map_smearing_2023BPix;

void cs_setter();
extern NamedFunc corrected_energy_lead;

NamedFunc assign_variation_lead_photon_relpterr_corrected(
    const NamedFunc &ph_sig, const NamedFunc &ph_pt, const std::string &name);

#endif
