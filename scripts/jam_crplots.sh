#plot_folder="./plots/an_cr-f -pull_plots"
plot_folder="./plots/lassen_v0_ALL_remake"
#plot_folder="./plots/kingscanyon_v1_sideband"
isrun3=${1:-""}
r3_str=""
r3_label="-run2"

if [ "$1" == "run3" ]; then
  r3_str="_run3"
  r3_label="-run3"

  plot_folder="./plots/lassen_v0_ALL_remake_run3"
fi 

cr_regions=("sideband" "ptyrat_l15o110" "Zfsrpeak" "lowIDMVA")
flavor=("ee" "mumu" "ll")
categories=("ggF" "VBF" "ttH_lep" "WH_3l" "ttH_had" "ZH_MET")

for cr in ${cr_regions[@]}
do
  for lep in ${flavor[@]}
  do
    control_region_mlly_plots="\
      ${plot_folder}/an_controlregions_cat_ggF_${lep}_${cr}${r3_str}_ggF__lly_m__mlly__wgt__lumi_lin.pdf \
      ${plot_folder}/an_controlregions_cat_VBF_${lep}_${cr}${r3_str}_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
      ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${cr}${r3_str}_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
      ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${cr}${r3_str}_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
      ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${cr}${r3_str}_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
      ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${cr}${r3_str}_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
      "

    ./scripts/pdf_combine.py -i $control_region_mlly_plots   -o plots/${cr}-${lep}${r3_label}-mlly-plots.pdf -x 3 -y 2 -f -p
  done
done


for cat in ${categories[@]}
do
  for lep in ${flavor[@]}
  do
    control_region_mlly_plots="\
    ${plot_folder}/an_controlregions_cat_${cat}_${lep}_sideband${r3_str}_${cat}__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_${cat}_${lep}_ptyrat_l15o110${r3_str}_${cat}__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_${cat}_${lep}_Zfsrpeak${r3_str}_${cat}__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_${cat}_${lep}_lowIDMVA${r3_str}_${cat}__lly_m__mlly__wgt__lumi_lin.pdf \
    "
    ./scripts/pdf_combine.py -i $control_region_mlly_plots -o plots/${cat}-controlregion-${lep}${r3_label}-mlly-plots.pdf -x 2 -y 2 -f -p

  done
done


########################################################################################################################
#ggF only control region plots
ggF_cr_regions=("sideband" "ptyrat_l15o110" "Zfsrpeak" "lowIDMVA")

for control_region in ${ggF_cr_regions[@]}
do
  for lep in ${flavor[@]}
  do

    ggF_control_regions_plots1_ee="\
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__phi_el1__el_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__phi_el2__el_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    ggF_control_regions_plots1_mumu="\

    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__phi_mu1__mu_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__phi_mu2__mu_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    ggF_control_regions_plots1_ll="\
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__phi_el1__el_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__phi_el2__el_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__phi_mu1__mu_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__phi_mu2__mu_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    ggF_control_regions_plots2="\
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__ll_m__ll_m0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__ll_rap__rap_ll__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__lly_pt__pT_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__lly_eta__eta_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__lly_rap__rap_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__lly_phi__phi_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF_cTHETA__lly_cosTheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF_ctheta__lly_costheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF_phi__lly_angphi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF__gamma_relunc__photon_energyErr0dphoton_pt0__wgt__lumi_lin.pdf \
    "

    ggF_control_regions_plots3="\
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_nel__nel__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_nmu__nmu__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_nll__nll__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_nphoton__nphoton__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_nllphoton__nllphoton__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_njet__njet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_nbdfm__nbdfm__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF_met__met__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF_j1_pt__j1_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF_j1_eta__j1_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF_j1_phi__j1_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_ggF_j1_m__j1_m__wgt__lumi_lin.pdf \
    "
    ggF_control_regions_plots4="\
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_ll_pt__ll_refit_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_ll_eta__ll_refit_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_ll_phi__ll_refit_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_ll_m__ll_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_lly_pt__llphoton_refit_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_lly_eta__rap_lly_refit__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_lly_phi__llphoton_refit_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_lly_m__llphoton_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_lly_ptom__llphoton_refit_ptdllphoton_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_cTHETA__llphoton_refit_cosTheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_ctheta__llphoton_refit_costheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_ll_y_dR__llphoton_refit_dr__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_ll_y_dphi__llphoton_refit_dphi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_pTt__llphoton_refit_pTt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_phi__llphoton_refit_psi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}${r3_str}_refit_y_ptom__photon_pt0dllphoton_refit_m__wgt__lumi_lin.pdf \
    "

    if [ "$lep" == "ee" ]; then
      ./scripts/pdf_combine.py -i $ggF_control_regions_plots1_ee -o plots/app-cr-ggf-${control_region}${r3_label}-${lep}-1.pdf -x 6 -y 2 -f -p
    elif [ "$lep" = "mumu" ]; then
      ./scripts/pdf_combine.py -i $ggF_control_regions_plots1_mumu -o plots/app-cr-ggf-${control_region}${r3_label}-${lep}-1.pdf -x 6 -y 2 -f -p
    else
      ./scripts/pdf_combine.py -i $ggF_control_regions_plots1_ll -o plots/app-cr-ggf-${control_region}${r3_label}-${lep}-1.pdf -x 6 -y 3 -f -p
    fi

    ./scripts/pdf_combine.py -i $ggF_control_regions_plots2 -o plots/app-cr-ggf-${control_region}${r3_label}-${lep}-2.pdf -x 5 -y 4 -f -p
    ./scripts/pdf_combine.py -i $ggF_control_regions_plots3 -o plots/app-cr-ggf-${control_region}${r3_label}-${lep}-3.pdf -x 4 -y 3 -f -p
    ./scripts/pdf_combine.py -i $ggF_control_regions_plots4 -o plots/app-cr-ggf-${control_region}${r3_label}-${lep}-4.pdf -x 4 -y 4 -f -p

   done
done



##########################################################################################################################

########################################################################################################################
VBF_cr_regions=("sideband" "ptyrat_l15o110" "Zfsrpeak" "lowIDMVA" "Zfsrpeak_NOT_dijet_m_g800" "Zfsrpeakdijet_m_g800" "sideband_NOT_dijet_m_g800" "sidebanddijet_m_g800")

for control_region in ${VBF_cr_regions[@]}
do
  for lep in ${flavor[@]}
  do

    VBF_control_regions_plots1_ee="\
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__phi_el1__el_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__phi_el2__el_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    VBF_control_regions_plots1_mumu="\
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__phi_mu1__mu_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__phi_mu2__mu_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    VBF_control_regions_plots1_ll="\
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__phi_el1__el_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__phi_el2__el_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__phi_mu1__mu_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__phi_mu2__mu_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "


    VBF_control_regions_plots2="\
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_j1_eta__j1_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_j1_phi__j1_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_j1_pt__j1_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_j2_eta__j2_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_j2_phi__j2_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_j2_pt__j2_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_deta__dijet_deta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_dphi__dijet_dphi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_dijet_m__dijet_m__wgt__lumi_lin.pdf \
    "

    VBF_control_regions_plots3="\
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__ll_m__ll_m0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__ll_rap__rap_ll__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__lly_pt__pT_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__lly_eta__eta_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__lly_rap__rap_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__lly_phi__phi_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF__gamma_relunc__photon_energyErr0dphoton_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_jjlly_dphi__dphi_lly_dijet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_sysbal__balance_lly_dijet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_yjdr__mindR_y_jet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_VBF_zep__zeppenfeld_y_jet__wgt__lumi_lin.pdf \
    "

    VBF_control_regions_plots4="\
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_nel__nel__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_nmu__nmu__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_nll__nll__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_nphoton__nphoton__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_nllphoton__nllphoton__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_njet__njet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_nbdfm__nbdfm__wgt__lumi_lin.pdf \
    "

    VBF_control_regions_plots5="\
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_ll_pt__ll_refit_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_ll_eta__ll_refit_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_ll_phi__ll_refit_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_ll_m__ll_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_lly_pt__llphoton_refit_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_lly_eta__rap_lly_refit__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_lly_phi__llphoton_refit_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_lly_m__llphoton_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_lly_ptom__llphoton_refit_ptdllphoton_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_cTHETA__llphoton_refit_cosTheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_ctheta__llphoton_refit_costheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_ll_y_dR__llphoton_refit_dr__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_ll_y_dphi__llphoton_refit_dphi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_pTt__llphoton_refit_pTt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_phi__llphoton_refit_psi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}${r3_str}_refit_y_ptom__photon_pt0dllphoton_refit_m__wgt__lumi_lin.pdf \
    "

    if [ "$lep" = "ee" ]; then
      ./scripts/pdf_combine.py -i $VBF_control_regions_plots1_ee -o plots/app-cr-vbf-${control_region}${r3_label}-${lep}-1.pdf -x 3 -y 4 -f -p
    elif [ "$lep" = "mumu" ]; then
      ./scripts/pdf_combine.py -i $VBF_control_regions_plots1_mumu -o plots/app-cr-vbf-${control_region}${r3_label}-${lep}-1.pdf -x 3 -y 4 -f -p
    else
      ./scripts/pdf_combine.py -i $VBF_control_regions_plots1_ll -o plots/app-cr-vbf-${control_region}${r3_label}-${lep}-1.pdf -x 6 -y 3 -f -p
    fi

    ./scripts/pdf_combine.py -i $VBF_control_regions_plots2 -o plots/app-cr-vbf-${control_region}${r3_label}-${lep}-2.pdf -x 3 -y 3 -f -p
    ./scripts/pdf_combine.py -i $VBF_control_regions_plots3 -o plots/app-cr-vbf-${control_region}${r3_label}-${lep}-3.pdf -x 5 -y 4 -f -p
    ./scripts/pdf_combine.py -i $VBF_control_regions_plots4 -o plots/app-cr-vbf-${control_region}${r3_label}-${lep}-4.pdf -x 4 -y 2 -f -p
    ./scripts/pdf_combine.py -i $VBF_control_regions_plots5 -o plots/app-cr-vbf-${control_region}${r3_label}-${lep}-5.pdf -x 4 -y 4 -f -p

  done
done



#VBF only control region plots
VBF_controlregions_meey_plots="\
${plot_folder}/an_controlregions_cat_VBF_ee_sideband_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_VBF_ee_ptyrat_l15o110_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_VBF_ee_Zfsrpeak_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_VBF_ee_lowIDMVA_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
"

VBF_controlregions_mmumuy_plots="\
${plot_folder}/an_controlregions_cat_VBF_mumu_sideband_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_VBF_mumu_ptyrat_l15o110_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_VBF_mumu_Zfsrpeak_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_VBF_mumu_lowIDMVA_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
"

VBF_controlregions_mlly_plots="\
${plot_folder}/an_controlregions_cat_VBF_ll_sideband_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_VBF_ll_ptyrat_l15o110_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_VBF_ll_Zfsrpeak_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_VBF_ll_lowIDMVA_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
"

#./scripts/pdf_combine.py -i $VBF_controlregion_meey_plots   -o plots/VBF_controlregion_ee_mlly_plots.pdf -x 3 -y 2 -f -p
#./scripts/pdf_combine.py -i $VBF_controlregion_mmumuy_plots -o plots/VBF_controlregion_mumu_mlly_plots.pdf -x 3 -y 2 -f -p
#./scripts/pdf_combine.py -i $VBF_controlregion_mlly_plots   -o plots/VBF_controlregion_ll_mlly_plots.pdf -x 3 -y 2 -f -p
##########################################################################################################################

########################################################################################################################
WH_3l_cr_regions=("sideband" "ptyrat_l15o110" "Zfsrpeak" "lowIDMVA")

for control_region in ${WH_3l_cr_regions[@]}
do
  for lep in ${flavor[@]}
  do

    WH_3l_control_regions_plots1_ee="\
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__phi_el1__el_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__phi_el2__el_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    WH_3l_control_regions_plots1_mumu="\
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__phi_mu1__mu_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__phi_mu2__mu_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    WH_3l_control_regions_plots1_ll="\
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__phi_el1__el_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__phi_el2__el_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__phi_mu1__mu_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__phi_mu2__mu_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    WH_3l_control_regions_plots2="\
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__ll_m__ll_m0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__ll_rap__rap_ll__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__lly_pt__pT_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__lly_eta__eta_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__lly_rap__rap_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__lly_phi__phi_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__gamma_relunc__photon_energyErr0dphoton_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    "

    WH_3l_control_regions_plots3="\
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_e3_pt__l3_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_e3_eta__l3_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_e3_phi__l3_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_mu3_pt__l3_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_mu3_eta__l3_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_mu3_phi__l3_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_e3idmva__l3_idmva__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_lep_l3_Imini__l3_mini__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_max_Imini__max_Imini__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_Zl3_DR__l3_Z_DR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_MET__met__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_dphi_Hmet__dphi_Hmet__wgt__lumi_lin.pdf \
    "

    WH_3l_control_regions_plots4="\
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_nel__nel__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_nmu__nmu__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_nll__nll__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_nphoton__nphoton__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_nllphoton__nllphoton__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_njet__njet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_nbdfm__nbdfl__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_nbdfm__nbdfm__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_j1_pt__j1_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_j1_eta__j1_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_j1_phi__j1_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_j1_m__j1_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_j2_pt__j2_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_j2_eta__j2_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_j2_phi__j2_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_j2_m__j2_m__wgt__lumi_lin.pdf \

    "

    WH_3l_control_regions_plots5="\
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_ll_pt__ll_refit_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_ll_eta__ll_refit_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_ll_phi__ll_refit_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_ll_m__ll_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_lly_pt__llphoton_refit_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_lly_eta__rap_lly_refit__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_lly_phi__llphoton_refit_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_lly_m__llphoton_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_lly_ptom__llphoton_refit_ptdllphoton_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_cTHETA__llphoton_refit_cosTheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_ctheta__llphoton_refit_costheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_ll_y_dR__llphoton_refit_dr__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_ll_y_dphi__llphoton_refit_dphi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_pTt__llphoton_refit_pTt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_phi__llphoton_refit_psi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_refit_y_ptom__photon_pt0dllphoton_refit_m__wgt__lumi_lin.pdf \
    "


    #${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_e3_phi__l3_phi__wgt__lumi_lin.pdf \
    #${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_mu3_phi__l3_phi__wgt__lumi_lin.pdf \
    #${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_e3_eta__l3_eta__wgt__lumi_lin.pdf \
    #${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}${r3_str}_WH_3l_e3_eta__l3_eta__wgt__lumi_lin.pdf \

    if [ "$lep" = "ee" ]; then
      ./scripts/pdf_combine.py -i $WH_3l_control_regions_plots1_ee -o plots/app-cr-vh3l-${control_region}${r3_label}-${lep}-1.pdf -x 3 -y 4 -f -p
    elif [ "$lep" = "mumu" ]; then
      ./scripts/pdf_combine.py -i $WH_3l_control_regions_plots1_mumu -o plots/app-cr-vh3l-${control_region}${r3_label}-${lep}-1.pdf -x 3 -y 4 -f -p
    else
      ./scripts/pdf_combine.py -i $WH_3l_control_regions_plots1_ll -o plots/app-cr-vh3l-${control_region}${r3_label}-${lep}-1.pdf -x 6 -y 3 -f -p
    fi

    ./scripts/pdf_combine.py -i $WH_3l_control_regions_plots2 -o plots/app-cr-vh3l-${control_region}${r3_label}-${lep}-2.pdf -x 5 -y 3 -f -p
    ./scripts/pdf_combine.py -i $WH_3l_control_regions_plots3 -o plots/app-cr-vh3l-${control_region}${r3_label}-${lep}-3.pdf -x 6 -y 2 -f -p
    ./scripts/pdf_combine.py -i $WH_3l_control_regions_plots4 -o plots/app-cr-vh3l-${control_region}${r3_label}-${lep}-4.pdf -x 4 -y 2 -f -p
    ./scripts/pdf_combine.py -i $WH_3l_control_regions_plots5 -o plots/app-cr-vh3l-${control_region}${r3_label}-${lep}-5.pdf -x 4 -y 4 -f -p

  done
done
#WH_3l only control region plots
WH_3l_controlregions_meey_plots="\
${plot_folder}/an_controlregions_cat_WH_3l_ee_sideband_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_WH_3l_ee_ptyrat_l15o110_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_WH_3l_ee_Zfsrpeak_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_WH_3l_ee_lowIDMVA_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
"

WH_3l_controlregions_mmumuy_plots="\
${plot_folder}/an_controlregions_cat_WH_3l_mumu_sideband_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_WH_3l_mumu_ptyrat_l15o110_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_WH_3l_mumu_Zfsrpeak_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_WH_3l_mumu_lowIDMVA_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
"

WH_3l_controlregions_mlly_plots="\
${plot_folder}/an_controlregions_cat_WH_3l_ll_sideband_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_WH_3l_ll_ptyrat_l15o110_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_WH_3l_ll_Zfsrpeak_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_WH_3l_ll_lowIDMVA_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
"

#./scripts/pdf_combine.py -i $WH_3l_controlregion_meey_plots   -o plots/WH_3l_controlregion_ee_mlly_plots.pdf -x 3 -y 2 -f -p
#./scripts/pdf_combine.py -i $WH_3l_controlregion_mmumuy_plots -o plots/WH_3l_controlregion_mumu_mlly_plots.pdf -x 3 -y 2 -f -p
#./scripts/pdf_combine.py -i $WH_3l_controlregion_mlly_plots   -o plots/WH_3l_controlregion_ll_mlly_plots.pdf -x 3 -y 2 -f -p
##########################################################################################################################

########################################################################################################################
#ZH_MET only control region plots

ZH_MET_cr_regions=("sideband" "ptyrat_l15o110" "Zfsrpeak" "lowIDMVA")

for control_region in ${ZH_MET_cr_regions[@]}
do
  for lep in ${flavor[@]}
  do

    ZH_MET_control_regions_plots1_ee="\
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__phi_el1__el_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__phi_el2__el_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    ZH_MET_control_regions_plots1_mumu="\
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__phi_mu1__mu_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__phi_mu2__mu_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    ZH_MET_control_regions_plots1_ll="\
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__phi_el1__el_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__phi_el2__el_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__phi_mu1__mu_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__phi_mu2__mu_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    ZH_MET_control_regions_plots2="\
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__ll_m__ll_m0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__ll_rap__rap_ll__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__lly_pt__pT_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__lly_eta__eta_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__lly_rap__rap_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__lly_phi__phi_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__gamma_relunc__photon_energyErr0dphoton_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    "

    ZH_MET_control_regions_plots3="\
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET_njet__njet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET_nbdfl__nbdfl__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET_nbdfm__nbdfm__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET_MET__met__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET_dphi_Hmet__dphi_Hmet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET_max_Imini__max_Imini__wgt__lumi_lin.pdf \
    "

    ZH_MET_control_regions_plots4="\
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_nel__nel__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_nmu__nmu__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_nll__nll__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_nphoton__nphoton__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_nllphoton__nllphoton__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_njet__njet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_nbdfm__nbdfl__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_nbdfm__nbdfm__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET_j1_pt__j1_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET_j1_eta__j1_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET_j1_phi__j1_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_ZH_MET_j1_m__j1_m__wgt__lumi_lin.pdf \
    "

    ZH_MET_control_regions_plots5="\
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_ll_pt__ll_refit_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_ll_eta__ll_refit_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_ll_phi__ll_refit_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_ll_m__ll_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_lly_pt__llphoton_refit_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_lly_eta__rap_lly_refit__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_lly_phi__llphoton_refit_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_lly_m__llphoton_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_lly_ptom__llphoton_refit_ptdllphoton_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_cTHETA__llphoton_refit_cosTheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_ctheta__llphoton_refit_costheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_ll_y_dR__llphoton_refit_dr__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_ll_y_dphi__llphoton_refit_dphi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_pTt__llphoton_refit_pTt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_phi__llphoton_refit_psi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}${r3_str}_refit_y_ptom__photon_pt0dllphoton_refit_m__wgt__lumi_lin.pdf \
    "

    if [ "$lep" = "ee" ]; then
      ./scripts/pdf_combine.py -i $ZH_MET_control_regions_plots1_ee -o plots/app-cr-zhptmiss-${control_region}${r3_label}-${lep}-1.pdf -x 3 -y 4 -f -p
    elif [ "$lep" = "mumu" ]; then
      ./scripts/pdf_combine.py -i $ZH_MET_control_regions_plots1_mumu -o plots/app-cr-zhptmiss-${control_region}${r3_label}-${lep}-1.pdf -x 3 -y 4 -f -p
    else
      ./scripts/pdf_combine.py -i $ZH_MET_control_regions_plots1_ll -o plots/app-cr-zhptmiss-${control_region}${r3_label}-${lep}-1.pdf -x 6 -y 3 -f -p
    fi

    ./scripts/pdf_combine.py -i $ZH_MET_control_regions_plots2 -o plots/app-cr-zhptmiss-${control_region}${r3_label}-${lep}-2.pdf -x 5 -y 3 -f -p
    ./scripts/pdf_combine.py -i $ZH_MET_control_regions_plots3 -o plots/app-cr-zhptmiss-${control_region}${r3_label}-${lep}-3.pdf -x 5 -y 2 -f -p
    ./scripts/pdf_combine.py -i $ZH_MET_control_regions_plots4 -o plots/app-cr-zhptmiss-${control_region}${r3_label}-${lep}-4.pdf -x 4 -y 2 -f -p
    ./scripts/pdf_combine.py -i $ZH_MET_control_regions_plots5 -o plots/app-cr-zhptmiss-${control_region}${r3_label}-${lep}-5.pdf -x 4 -y 4 -f -p

  done
done

ZH_MET_controlregions_meey_plots="\
${plot_folder}/an_controlregions_cat_ZH_MET_ee_sideband_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ZH_MET_ee_ptyrat_l15o110_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ZH_MET_ee_Zfsrpeak_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ZH_MET_ee_lowIDMVA_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
"

ZH_MET_controlregions_mmumuy_plots="\
${plot_folder}/an_controlregions_cat_ZH_MET_mumu_sideband_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ZH_MET_mumu_ptyrat_l15o110_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ZH_MET_mumu_Zfsrpeak_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ZH_MET_mumu_lowIDMVA_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
"

ZH_MET_controlregions_mlly_plots="\
${plot_folder}/an_controlregions_cat_ZH_MET_ll_sideband_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ZH_MET_ll_ptyrat_l15o110_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ZH_MET_ll_Zfsrpeak_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ZH_MET_ll_lowIDMVA_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
"

#./scripts/pdf_combine.py -i $ZH_MET_controlregion_meey_plots   -o plots/ZH_MET_controlregion_ee_mlly_plots.pdf -x 3 -y 2 -f -p
#./scripts/pdf_combine.py -i $ZH_MET_controlregion_mmumuy_plots -o plots/ZH_MET_controlregion_mumu_mlly_plots.pdf -x 3 -y 2 -f -p
#./scripts/pdf_combine.py -i $ZH_MET_controlregion_mlly_plots   -o plots/ZH_MET_controlregion_ll_mlly_plots.pdf -x 3 -y 2 -f -p
##########################################################################################################################

########################################################################################################################
#ttH_had only control region plots
ttH_had_cr_regions=("sideband" "ptyrat_l15o110" "Zfsrpeak" "lowIDMVA" "Zfsrpeak_NOT_metg50" "Zfsrpeakmetg50")

for control_region in ${ttH_had_cr_regions[@]}
do
  for lep in ${flavor[@]}
  do

    ttH_had_control_regions_plots1_ee="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__phi_el1__el_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__phi_el2__el_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    ttH_had_control_regions_plots1_mumu="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__phi_mu1__mu_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__phi_mu2__mu_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    ttH_had_control_regions_plots1_ll="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__phi_el1__el_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__phi_el2__el_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__phi_mu1__mu_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__phi_mu2__mu_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    ttH_had_control_regions_plots2="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__ll_m__ll_m0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__ll_rap__rap_ll__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__lly_pt__pT_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__lly_eta__eta_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__lly_rap__rap_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__lly_phi__phi_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__gamma_relunc__photon_energyErr0dphoton_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    "

    ttH_had_control_regions_plots3="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_Njet__njet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_HT__ht__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_MET__met__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_max_Imini__max_Imini__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_ptj_m4j__ptj_m4j__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_ptj_m5j__ptj_m5j__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    "
    #${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_deepflav_BASE_MAX__max_deepflav__wgt__lumi_lin.pdf \ 

    ttH_had_control_regions_plots4="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_nel__nel__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_nmu__nmu__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_nll__nll__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_nphoton__nphoton__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_nllphoton__nllphoton__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_njet__njet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_nbdfm__nbdfm__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_nbdfm__nbdft__wgt__lumi_lin.pdf \
    "

    ttH_had_control_regions_plots5="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_ll_pt__ll_refit_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_ll_eta__ll_refit_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_ll_phi__ll_refit_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_ll_m__ll_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_lly_pt__llphoton_refit_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_lly_eta__rap_lly_refit__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_lly_phi__llphoton_refit_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_lly_m__llphoton_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_lly_ptom__llphoton_refit_ptdllphoton_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_cTHETA__llphoton_refit_cosTheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_ctheta__llphoton_refit_costheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_ll_y_dR__llphoton_refit_dr__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_ll_y_dphi__llphoton_refit_dphi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_pTt__llphoton_refit_pTt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_phi__llphoton_refit_psi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_refit_y_ptom__photon_pt0dllphoton_refit_m__wgt__lumi_lin.pdf \
    "

    ttH_had_control_regions_plots6="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j1_pt__j1_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j1_eta__j1_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j1_phi__j1_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j1_m__j1_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j2_pt__j2_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j2_eta__j2_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j2_phi__j2_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j2_m__j2_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j3_pt__j3_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j3_eta__j3_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j3_phi__j3_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j3_m__j3_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j4_pt__j4_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j4_eta__j4_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j4_phi__j4_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j4_m__j4_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j5_pt__j5_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j5_eta__j5_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j5_phi__j5_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}${r3_str}_ttH_had_j5_m__j5_m__wgt__lumi_lin.pdf \
    "

    if [ "$lep" = "ee" ]; then
      ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots1_ee -o plots/app-cr-tthhad-${control_region}${r3_label}-${lep}-1.pdf -x 3 -y 4 -f -p
    elif [ "$lep" = "mumu" ]; then
      ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots1_mumu -o plots/app-cr-tthhad-${control_region}${r3_label}-${lep}-1.pdf -x 3 -y 4 -f -p
    else
      ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots1_ll -o plots/app-cr-tthhad-${control_region}${r3_label}-${lep}-1.pdf -x 6 -y 3 -f -p
    fi

    ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots2 -o plots/app-cr-tthhad-${control_region}${r3_label}-${lep}-2.pdf -x 5 -y 3 -f -p
    ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots3 -o plots/app-cr-tthhad-${control_region}${r3_label}-${lep}-3.pdf -x 5 -y 2 -f -p
    ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots4 -o plots/app-cr-tthhad-${control_region}${r3_label}-${lep}-4.pdf -x 5 -y 2 -f -p
    ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots5 -o plots/app-cr-tthhad-${control_region}${r3_label}-${lep}-5.pdf -x 4 -y 4 -f -p
    ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots6 -o plots/app-cr-tthhad-${control_region}${r3_label}-${lep}-6.pdf -x 4 -y 5 -f -p

  done
done

ttH_had_controlregions_meey_plots="\
${plot_folder}/an_controlregions_cat_ttH_had_ee_sideband_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_had_ee_ptyrat_l15o110_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_had_ee_Zfsrpeak_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_had_ee_lowIDMVA_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
"

ttH_had_controlregions_mmumuy_plots="\
${plot_folder}/an_controlregions_cat_ttH_had_mumu_sideband_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_had_mumu_ptyrat_l15o110_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_had_mumu_Zfsrpeak_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_had_mumu_lowIDMVA_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
"

ttH_had_controlregions_mlly_plots="\
${plot_folder}/an_controlregions_cat_ttH_had_ll_sideband_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_had_ll_ptyrat_l15o110_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_had_ll_Zfsrpeak_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_had_ll_lowIDMVA_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
"

#./scripts/pdf_combine.py -i $ttH_had_controlregion_meey_plots   -o plots/ttH_had_controlregion_ee_mlly_plots.pdf -x 3 -y 2 -f -p
#./scripts/pdf_combine.py -i $ttH_had_controlregion_mmumuy_plots -o plots/ttH_had_controlregion_mumu_mlly_plots.pdf -x 3 -y 2 -f -p
#./scripts/pdf_combine.py -i $ttH_had_controlregion_mlly_plots   -o plots/ttH_had_controlregion_ll_mlly_plots.pdf -x 3 -y 2 -f -p
##########################################################################################################################

########################################################################################################################
#ttH_lep only control region plots
ttH_lep_cr_regions=("sideband" "ptyrat_l15o110" "Zfsrpeak" "lowIDMVA" "Zfsrpeak_NOT_metg50" "Zfsrpeakmetg50" "sideband_NOT_metg50" "sidebandmetg50")

for control_region in ${ttH_lep_cr_regions[@]}
do
  for lep in ${flavor[@]}
  do

    ttH_lep_control_regions_plots1_ee="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__phi_el1__el_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__phi_el2__el_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    ttH_lep_control_regions_plots1_mumu="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__phi_mu1__mu_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__phi_mu2__mu_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "

    ttH_lep_control_regions_plots1_ll="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__phi_el1__el_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__phi_el2__el_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__phi_mu1__mu_phill_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__phi_mu2__mu_phill_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    "
    ttH_lep_control_regions_plots2="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__ll_m__ll_m0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__ll_rap__rap_ll__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__lly_pt__pT_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__lly_eta__eta_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__lly_rap__rap_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__lly_phi__phi_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__gamma_relunc__photon_energyErr0dphoton_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    "

    ttH_lep_control_regions_plots3="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_e3_pt__l3_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_e3_eta__l3_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_e3_phi__l3_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_mu3_pt__l3_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_mu3_eta__l3_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_mu3_phi__l3_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_l3_pt__l3_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_l3_eta__l3_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_l3_phi__l3_phi__wgt__lumi_lin.pdf \
    "

    ttH_lep_control_regions_plots4="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_Nlep__nlep__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_Njet__njet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_Nbdfm__nbdfm__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_MET__met__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_HT__ht__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_e3idmva__l3_idmva__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_l3_Imini__l3_mini__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_max_Imini__max_Imini__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_H_Z_dR__H_Z_dR__wgt__lumi_lin.pdf \
    "

    ttH_lep_control_regions_plots5="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_nel__nel__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_nmu__nmu__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_nll__nll__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_nphoton__nphoton__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_nllphoton__nllphoton__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_njet__njet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_nbdfm__nbdfm__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_nbdfm__nbdft__wgt__lumi_lin.pdf \
    "

    ttH_lep_control_regions_plots6="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_ll_pt__ll_refit_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_ll_eta__ll_refit_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_ll_phi__ll_refit_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_ll_m__ll_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_lly_pt__llphoton_refit_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_lly_eta__rap_lly_refit__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_lly_phi__llphoton_refit_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_lly_m__llphoton_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_lly_ptom__llphoton_refit_ptdllphoton_refit_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_cTHETA__llphoton_refit_cosTheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_ctheta__llphoton_refit_costheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_ll_y_dR__llphoton_refit_dr__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_ll_y_dphi__llphoton_refit_dphi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_pTt__llphoton_refit_pTt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_phi__llphoton_refit_psi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_refit_y_ptom__photon_pt0dllphoton_refit_m__wgt__lumi_lin.pdf \
    "

    ttH_lep_control_regions_plots7="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_j1_pt__j1_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_j1_eta__j1_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_j1_phi__j1_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_j1_m__j1_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_j2_pt__j2_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_j2_eta__j2_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_j2_phi__j2_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_j2_m__j2_m__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_j3_pt__j3_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_j3_eta__j3_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_j3_phi__j3_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}${r3_str}_ttH_lep_j3_m__j3_m__wgt__lumi_lin.pdf \
    "

    if [ "$lep" = "ee" ]; then
      ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots1_ee -o plots/app-cr-tthlep-${control_region}${r3_label}-${lep}-1.pdf -x 3 -y 4 -f -p
    elif [ "$lep" = "mumu" ]; then
      ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots1_mumu -o plots/app-cr-tthlep-${control_region}${r3_label}-${lep}-1.pdf -x 3 -y 4 -f -p
    else
      ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots1_ll -o plots/app-cr-tthlep-${control_region}${r3_label}-${lep}-1.pdf -x 6 -y 3 -f -p
    fi

    ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots2 -o plots/app-cr-tthlep-${control_region}${r3_label}-${lep}-2.pdf -x 5 -y 3 -f -p
    ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots3 -o plots/app-cr-tthlep-${control_region}${r3_label}-${lep}-3.pdf -x 3 -y 3 -f -p
    ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots4 -o plots/app-cr-tthlep-${control_region}${r3_label}-${lep}-4.pdf -x 4 -y 3 -f -p
    ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots5 -o plots/app-cr-tthlep-${control_region}${r3_label}-${lep}-5.pdf -x 4 -y 2 -f -p
    ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots6 -o plots/app-cr-tthlep-${control_region}${r3_label}-${lep}-6.pdf -x 4 -y 4 -f -p
    ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots7 -o plots/app-cr-tthlep-${control_region}${r3_label}-${lep}-7.pdf -x 4 -y 3 -f -p

  done
done



ttH_lep_controlregions_meey_plots="\
${plot_folder}/an_controlregions_cat_ttH_lep_ee_sideband_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_lep_ee_ptyrat_l15o110_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_lep_ee_Zfsrpeak_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_lep_ee_lowIDMVA_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
"

ttH_lep_controlregions_mmumuy_plots="\
${plot_folder}/an_controlregions_cat_ttH_lep_mumu_sideband_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_lep_mumu_ptyrat_l15o110_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_lep_mumu_Zfsrpeak_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_lep_mumu_lowIDMVA_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
"

ttH_lep_controlregions_mlly_plots="\
${plot_folder}/an_controlregions_cat_ttH_lep_ll_sideband_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_lep_ll_ptyrat_l15o110_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_lep_ll_Zfsrpeak_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
${plot_folder}/an_controlregions_cat_ttH_lep_ll_lowIDMVA_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
"

#./scripts/pdf_combine.py -i $ttH_lep_controlregion_meey_plots   -o plots/ttH_lep_controlregion_ee_mlly_plots.pdf -x 3 -y 2 -f -p
#./scripts/pdf_combine.py -i $ttH_lep_controlregion_mmumuy_plots -o plots/ttH_lep_controlregion_mumu_mlly_plots.pdf -x 3 -y 2 -f -p
#./scripts/pdf_combine.py -i $ttH_lep_controlregion_mlly_plots   -o plots/ttH_lep_controlregion_ll_mlly_plots.pdf -x 3 -y 2 -f -p
##########################################################################################################################


