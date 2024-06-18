plot_folder="./plots/an_cr_full_plots"

cr_regions=("sideband" "ptyrat_l15o110" "Zfsrpeak" "lowIDMVA")
flavor=("ee" "mumu" "ll")
categories=("ggF" "VBF" "WH_3l" "ZH_MET" "ttH_had" "ttH_lep")

for cr in ${cr_regions[@]}
do
  for lep in ${flavor[@]}
  do
    control_region_mlly_plots="\
      ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${cr}_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
      ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${cr}_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
      ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${cr}_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
      ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${cr}_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
      ${plot_folder}/an_controlregions_cat_VBF_${lep}_${cr}_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
      ${plot_folder}/an_controlregions_cat_ggF_${lep}_${cr}_ggF__lly_m__mlly__wgt__lumi_lin.pdf \
    "

    ./scripts/pdf_combine.py -i $control_region_mlly_plots   -o plots/${cr}_${lep}_mlly_plots.pdf -x 3 -y 2 -f

  done
done


for cat in ${categories[@]}
do
  for lep in ${flavor[@]}
  do
    control_region_mlly_plots="\
    ${plot_folder}/an_controlregions_cat_${cat}_${lep}_sideband_${cat}__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_${cat}_${lep}_ptyrat_l15o110_${cat}__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_${cat}_${lep}_Zfsrpeak_${cat}__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_${cat}_${lep}_lowIDMVA_${cat}__lly_m__mlly__wgt__lumi_lin.pdf \
    "
    ./scripts/pdf_combine.py -i $control_region_mlly_plots -o plots/${cat}_controlregion_${lep}_mlly_plots.pdf -x 2 -y 2 -f

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
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__10_pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__10_pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__11_eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__11_eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_m__ll_m0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    "

    ggF_control_regions_plots1_mumu="\
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__10_pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__10_pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__11_eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__11_eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_m__ll_m0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    "

    ggF_control_regions_plots1_ll="\
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__10_pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__10_pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__11_eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__11_eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__10_pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__10_pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__11_eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__11_eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_m__ll_m0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    "

    ggF_control_regions_plots2="\
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_cTHETA__llphoton_cosTheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ctheta__llphoton_costheta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_phi__llphoton_psi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__gamma_relunc__photon_energyErr0dphoton_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__lly_pt__pT_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__lly_eta__eta_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__lly_eta__rap_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ggF_${lep}_${control_region}_ggF__lly_phi__phi_lly__wgt__lumi_lin.pdf \
    "
 
    if [["${lep}"=="ee"]]; then
      ./scripts/pdf_combine.py -i $ggF_control_regions_plots1_ee -o plots/ggF_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 3 -f
    elif[["${lep}"=="mumu"]]; then
      ./scripts/pdf_combine.py -i $ggF_control_regions_plots1_mumu -o plots/ggF_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 3 -f
    else
      ./scripts/pdf_combine.py -i $ggF_control_regions_plots1_ll -o plots/ggF_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 4 -f
    fi

    ./scripts/pdf_combine.py -i $ggF_control_regions_plots2 -o plots/ggF_cr_${control_region}_${lep}_plots_part2.pdf -x 4 -y 4 -f
  done
done



##########################################################################################################################

########################################################################################################################
VBF_cr_regions=("sideband")
#"ptyrat_l15o110" "Zfsrpeak" "lowIDMVA")


for control_region in ${VBF_cr_regions[@]}
do
  for lep in ${flavor[@]}
  do

    VBF_control_regions_plots1_ee="\
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__10_pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__10_pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__11_eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__11_eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_m__ll_m0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    "

    VBF_control_regions_plots1_mumu="\
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__10_pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__10_pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__11_eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__11_eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_m__ll_m0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    "

    VBF_control_regions_plots1_ll="\
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__10_pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__10_pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__11_eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__11_eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__10_pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__10_pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__11_eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__11_eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_m__ll_m0__wgt__lumi_lin.pdf \
    "

    VBF_control_regions_plots2="\
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_j1_eta__j1_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_j1_phi__j1_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_j1_pt__j1_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_j2_eta__j2_eta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_j2_phi__j2_phi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_j2_pt__j2_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_deta__dijet_deta__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_dphi__dijet_dphi__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_dijet_m__dijet_m__wgt__lumi_lin.pdf \
    "

    VBF_control_regions_plots3="\
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__lly_pt__pT_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__lly_eta__eta_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__lly_eta__rap_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__lly_phi__phi_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__gamma_relunc__photon_energyErr0dphoton_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_jjlly_dphi__llphoton_dijet_dphi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_sysbal__llphoton_dijet_balance0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_yjdr__photon_jet_mindr0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_VBF_${lep}_${control_region}_VBF_zep__photon_zeppenfeld0__wgt__lumi_lin.pdf \
    "


    if [["${lep}"=="ee"]]; then
      ./scripts/pdf_combine.py -i $VBF_control_regions_plots1_ee -o plots/VBF_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 3 -f
    elif[["${lep}"=="mumu"]]; then
      ./scripts/pdf_combine.py -i $VBF_control_regions_plots1_mumu -o plots/VBF_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 3 -f
    else
      ./scripts/pdf_combine.py -i $VBF_control_regions_plots1_ll -o plots/VBF_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 4 -f
    fi

    ./scripts/pdf_combine.py -i $VBF_control_regions_plots2 -o plots/VBF_cr_${control_region}_${lep}_plots_part2.pdf -x 3 -y 3 -f
    ./scripts/pdf_combine.py -i $VBF_control_regions_plots3 -o plots/VBF_cr_${control_region}_${lep}_plots_part3.pdf -x 4 -y 4 -f

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

#./scripts/pdf_combine.py -i $VBF_controlregion_meey_plots   -o plots/VBF_controlregion_ee_mlly_plots.pdf -x 3 -y 2 -f
#./scripts/pdf_combine.py -i $VBF_controlregion_mmumuy_plots -o plots/VBF_controlregion_mumu_mlly_plots.pdf -x 3 -y 2 -f
#./scripts/pdf_combine.py -i $VBF_controlregion_mlly_plots   -o plots/VBF_controlregion_ll_mlly_plots.pdf -x 3 -y 2 -f
##########################################################################################################################

########################################################################################################################
WH_3l_cr_regions=("sideband")
#"ptyrat_l15o110" "Zfsrpeak" "lowIDMVA")

for control_region in ${WH_3l_cr_regions[@]}
do
  for lep in ${flavor[@]}
  do

    WH_3l_control_regions_plots1_ee="\
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__10_pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__10_pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__11_eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__11_eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l_e3_pt__l3_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l_mu3_pt__l3_pt__wgt__lumi_lin.pdf \
    "

    WH_3l_control_regions_plots1_mumu="\
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__10_pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__10_pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__11_eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__11_eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l_e3_pt__l3_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l_mu3_pt__l3_pt__wgt__lumi_lin.pdf \
    "

    WH_3l_control_regions_plots1_ll="\
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__10_pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__10_pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__11_eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__11_eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__10_pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__10_pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__11_eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__11_eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l_e3_pt__l3_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l_mu3_pt__l3_pt__wgt__lumi_lin.pdf \
    "

    WH_3l_control_regions_plots2="\
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__ll_m__ll_m0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__lly_pt__pT_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__lly_eta__eta_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__lly_eta__rap_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__lly_phi__phi_lly__wgt__lumi_lin.pdf \
    "

    WH_3l_control_regions_plots3="\
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l_Zl3_DR__l3_Z_DR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l_MET__met__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__gamma_relunc__photon_energyErr0dphoton_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l_dphi_Hmet__dphi_Hmet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l_e3idmva__l3_idmva__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l_lep_l3_Imini__l3_mini__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_WH_3l_${lep}_${control_region}_WH_3l_max_Imini__max_Imini__wgt__lumi_lin.pdf \
    "


    if [["${lep}"=="ee"]]; then
      ./scripts/pdf_combine.py -i $WH_3l_control_regions_plots1_ee -o plots/WH_3l_cr_${control_region}_${lep}_plots_part1.pdf -x 3 -y 2 -f
    elif[["${lep}"=="mumu"]]; then
      ./scripts/pdf_combine.py -i $WH_3l_control_regions_plots1_mumu -o plots/WH_3l_cr_${control_region}_${lep}_plots_part1.pdf -x 3 -y 2 -f
    else
      ./scripts/pdf_combine.py -i $WH_3l_control_regions_plots1_ll -o plots/WH_3l_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 3 -f
    fi

    ./scripts/pdf_combine.py -i $WH_3l_control_regions_plots2 -o plots/WH_3l_cr_${control_region}_${lep}_plots_part2.pdf -x 3 -y 3 -f
    ./scripts/pdf_combine.py -i $WH_3l_control_regions_plots3 -o plots/WH_3l_cr_${control_region}_${lep}_plots_part3.pdf -x 4 -y 4 -f

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

#./scripts/pdf_combine.py -i $WH_3l_controlregion_meey_plots   -o plots/WH_3l_controlregion_ee_mlly_plots.pdf -x 3 -y 2 -f
#./scripts/pdf_combine.py -i $WH_3l_controlregion_mmumuy_plots -o plots/WH_3l_controlregion_mumu_mlly_plots.pdf -x 3 -y 2 -f
#./scripts/pdf_combine.py -i $WH_3l_controlregion_mlly_plots   -o plots/WH_3l_controlregion_ll_mlly_plots.pdf -x 3 -y 2 -f
##########################################################################################################################

########################################################################################################################
#ZH_MET only control region plots

ZH_MET_cr_regions=("sideband")
#"ptyrat_l15o110" "Zfsrpeak" "lowIDMVA")

for control_region in ${ZH_MET_cr_regions[@]}
do
  for lep in ${flavor[@]}
  do

    ZH_MET_control_regions_plots1_ee="\
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__10_pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__10_pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__11_eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__11_eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_m__ll_m0__wgt__lumi_lin.pdf \
    "

    ZH_MET_control_regions_plots1_mumu="\
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__10_pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__10_pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__11_eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__11_eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_m__ll_m0__wgt__lumi_lin.pdf \
    "

    ZH_MET_control_regions_plots1_ll="\
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__10_pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__10_pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__11_eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__11_eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__10_pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__10_pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__11_eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__11_eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_m__ll_m0__wgt__lumi_lin.pdf \
    "

    ZH_MET_control_regions_plots2="\
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__gamma_relunc__photon_energyErr0dphoton_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__lly_pt__pT_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__lly_eta__eta_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__lly_phi__phi_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__lly_eta__rap_lly__wgt__lumi_lin.pdf \
    "

    ZH_MET_control_regions_plots3="\
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET_MET__met__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET_dphi_Hmet__dphi_Hmet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET_max_Imini__max_Imini__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET_nbdfl__nbdfl__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET_nbdfm__nbdfm__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET_njet__njet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ZH_MET_${lep}_${control_region}_ZH_MET_pTt__llgamma_pTt__wgt__lumi_lin.pdf \
    "

    if [["${lep}"=="ee"]]; then
      ./scripts/pdf_combine.py -i $ZH_MET_control_regions_plots1_ee -o plots/ZH_MET_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 3 -f
    elif[["${lep}"=="mumu"]]; then
      ./scripts/pdf_combine.py -i $ZH_MET_control_regions_plots1_mumu -o plots/ZH_MET_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 3 -f
    else
      ./scripts/pdf_combine.py -i $ZH_MET_control_regions_plots1_ll -o plots/ZH_MET_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 3 -f
    fi

    ./scripts/pdf_combine.py -i $ZH_MET_control_regions_plots2 -o plots/ZH_MET_cr_${control_region}_${lep}_plots_part2.pdf -x 4 -y 3 -f
    ./scripts/pdf_combine.py -i $ZH_MET_control_regions_plots3 -o plots/ZH_MET_cr_${control_region}_${lep}_plots_part3.pdf -x 4 -y 3 -f


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

#./scripts/pdf_combine.py -i $ZH_MET_controlregion_meey_plots   -o plots/ZH_MET_controlregion_ee_mlly_plots.pdf -x 3 -y 2 -f
#./scripts/pdf_combine.py -i $ZH_MET_controlregion_mmumuy_plots -o plots/ZH_MET_controlregion_mumu_mlly_plots.pdf -x 3 -y 2 -f
#./scripts/pdf_combine.py -i $ZH_MET_controlregion_mlly_plots   -o plots/ZH_MET_controlregion_ll_mlly_plots.pdf -x 3 -y 2 -f
##########################################################################################################################

########################################################################################################################
#ttH_had only control region plots
ttH_had_cr_regions=("sideband")
#"ptyrat_l15o110" "Zfsrpeak" "lowIDMVA")

for control_region in ${ttH_had_cr_regions[@]}
do
  for lep in ${flavor[@]}
  do

    ttH_had_control_regions_plots1_ee="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__10_pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__10_pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__11_eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__11_eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_m__ll_m0__wgt__lumi_lin.pdf \
    "

    ttH_had_control_regions_plots1_mumu="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__10_pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__10_pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__11_eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__11_eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_m__ll_m0__wgt__lumi_lin.pdf \
    "

    ttH_had_control_regions_plots1_ll="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__10_pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__10_pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__11_eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__11_eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__10_pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__10_pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__11_eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__11_eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_m__ll_m0__wgt__lumi_lin.pdf \
    "

    ttH_had_control_regions_plots2="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__gamma_relunc__photon_energyErr0dphoton_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__lly_pt__pT_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__lly_eta__eta_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__lly_eta__rap_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__lly_phi__phi_lly__wgt__lumi_lin.pdf \
    "

    ttH_had_control_regions_plots3="\
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had_HT__ht__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had_MET__met__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had_Njet__njet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had_max_Imini__max_Imini__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had_deepflav_BASE_MAX__max_deepflav__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had_ptj_m4j__ptj_m4j__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_had_${lep}_${control_region}_ttH_had_ptj_m5j__ptj_m5j__wgt__lumi_lin.pdf \
    "

    if [["${lep}"=="ee"]]; then
      ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots1_ee -o plots/ttH_had_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 2 -f
    elif[["${lep}"=="mumu"]]; then
      ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots1_mumu -o plots/ttH_had_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 2 -f
    else
      ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots1_ll -o plots/ttH_had_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 3 -f
    fi

    ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots2 -o plots/ttH_had_cr_${control_region}_${lep}_plots_part2.pdf -x 4 -y 3 -f
    ./scripts/pdf_combine.py -i $ttH_had_control_regions_plots3 -o plots/ttH_had_cr_${control_region}_${lep}_plots_part3.pdf -x 4 -y 3 -f

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

#./scripts/pdf_combine.py -i $ttH_had_controlregion_meey_plots   -o plots/ttH_had_controlregion_ee_mlly_plots.pdf -x 3 -y 2 -f
#./scripts/pdf_combine.py -i $ttH_had_controlregion_mmumuy_plots -o plots/ttH_had_controlregion_mumu_mlly_plots.pdf -x 3 -y 2 -f
#./scripts/pdf_combine.py -i $ttH_had_controlregion_mlly_plots   -o plots/ttH_had_controlregion_ll_mlly_plots.pdf -x 3 -y 2 -f
##########################################################################################################################

########################################################################################################################
#ttH_lep only control region plots
ttH_lep_cr_regions=("sideband")
#"ptyrat_l15o110" "Zfsrpeak" "lowIDMVA")

for control_region in ${ttH_lep_cr_regions[@]}
do
  for lep in ${flavor[@]}
  do
    ttH_lep_control_regions_plots1_ee="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__10_pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__10_pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__11_eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__11_eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_m__ll_m0__wgt__lumi_lin.pdf \
    "

    ttH_lep_control_regions_plots1_mumu="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__10_pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__10_pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__11_eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__11_eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_m__ll_m0__wgt__lumi_lin.pdf \
    "

    ttH_lep_control_regions_plots1_ll="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__10_pt_el1__el_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__10_pt_el2__el_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__11_eta_el1__el_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__11_eta_el2__el_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__10_pt_mu1__mu_ptll_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__10_pt_mu2__mu_ptll_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__11_eta_mu1__mu_etall_i10__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__11_eta_mu2__mu_etall_i20__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_pt__ll_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_eta__ll_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_phi__ll_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_m__ll_m0__wgt__lumi_lin.pdf \
    "


    ttH_lep_control_regions_plots2="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__lly_m__mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__gamma_pT__photon_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__gamma_eta__photon_eta0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__gamma_phi__photon_phi0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__gamma_maxDR__photon_drmax0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__gamma_minDR__photon_drmin0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__gamma_IDMVA__photon_idmva0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__gamma_relunc__photon_energyErr0dphoton_pt0__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__lly_pt__pT_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__lly_eta__eta_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__lly_eta__rap_lly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__lly_phi__phi_lly__wgt__lumi_lin.pdf \
    "

    ttH_lep_control_regions_plots3="\
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__ll_y_dR__lly_dR__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__lly_ptom__pTlly_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep__y_ptom__pTy_mlly__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep_e3_pt__l3_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep_e3idmva__l3_idmva__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep_mu3_pt__l3_pt__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep_l3_Imini__l3_mini__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep_max_Imini__max_Imini__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep_HT__ht__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep_MET__met__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep_Nbdfm__nbdfm__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep_Njet__njet__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep_Nlep__nlep__wgt__lumi_lin.pdf \
    ${plot_folder}/an_controlregions_cat_ttH_lep_${lep}_${control_region}_ttH_lep_H_Z_dR__H_Z_dR__wgt__lumi_lin.pdf \
    "

    if [["${lep}"=="ee"]]; then
      ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots1 -o plots/ttH_lep_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 2 -f
    elif[["${lep}"=="mumu"]]; then
      ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots1 -o plots/ttH_lep_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 2 -f
    else
      ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots1 -o plots/ttH_lep_cr_${control_region}_${lep}_plots_part1.pdf -x 4 -y 3 -f
    fi

    ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots2 -o plots/ttH_lep_cr_${control_region}_${lep}_plots_part2.pdf -x 4 -y 3 -f
    ./scripts/pdf_combine.py -i $ttH_lep_control_regions_plots3 -o plots/ttH_lep_cr_${control_region}_${lep}_plots_part3.pdf -x 4 -y 4 -f

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

#./scripts/pdf_combine.py -i $ttH_lep_controlregion_meey_plots   -o plots/ttH_lep_controlregion_ee_mlly_plots.pdf -x 3 -y 2 -f
#./scripts/pdf_combine.py -i $ttH_lep_controlregion_mmumuy_plots -o plots/ttH_lep_controlregion_mumu_mlly_plots.pdf -x 3 -y 2 -f
#./scripts/pdf_combine.py -i $ttH_lep_controlregion_mlly_plots   -o plots/ttH_lep_controlregion_ll_mlly_plots.pdf -x 3 -y 2 -f
##########################################################################################################################


