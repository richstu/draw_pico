#quickly get scale of systematic unceratinties

from statistics import mean, stdev

SYSTS = ['pdf_Higgs_gg', 'pdf_Higgs_qqbar', 'pdf_Higgs_ttH', 'QCD_scale_ggH', 'QCD_scale_qqH', 'QCD_scale_VH', 'QCD_scale_ttH', 'QCD_scale_diff_ggH', 'QCD_scale_diff_qqH', 'QCD_scale_diff_VH', 'QCD_scale_diff_ttH', 'BR_hzg', 'BR_hmm', 'param_alphas', 'param_mq', 'ps_isr', 'ps_fsr', 'underlying_event', 'lumi_13TeV', 'lumi_2022', 'lumi_2023', 'CMS_pileup', 'CMS_l1_ecal_prefiring', 'CMS_eff_e', 'CMS_eff_m', 'CMS_eff_g', 'CMS_trigger_e', 'CMS_trigger_m', 'CMS_quality_g', 'CMS_btag_heavy', 'CMS_btag_light', 'CMS_btag_heavy_2016APV', 'CMS_btag_light_2016APV', 'CMS_btag_heavy_2016', 'CMS_btag_light_2016', 'CMS_btag_heavy_2017', 'CMS_btag_light_2017', 'CMS_btag_heavy_2018', 'CMS_btag_light_2018', 'CMS_btag_heavy_2022', 'CMS_btag_light_2022', 'CMS_btag_heavy_2022EE', 'CMS_btag_light_2022EE', 'CMS_btag_heavy_2023', 'CMS_btag_light_2023', 'CMS_btag_heavy_2023BPix', 'CMS_btag_light_2023BPix', 'CMS_scale_e', 'CMS_res_e', 'CMS_res_m', 'CMS_scale_g', 'CMS_res_g', 'CMS_scale_j_2016APV', 'CMS_res_j_2016APV', 'CMS_scale_j_2016', 'CMS_res_j_2016', 'CMS_scale_j_2017', 'CMS_res_j_2017', 'CMS_scale_j_2018', 'CMS_res_j_2018', 'CMS_scale_j_2022', 'CMS_res_j_2022', 'CMS_scale_j_2022EE', 'CMS_res_j_2022EE', 'CMS_scale_j_2023', 'CMS_res_j_2023', 'CMS_scale_j_2023BPix', 'CMS_res_j_2023BPix']
MCSTAT_SYSTS = ['MCstat_Htozg_el_cat_ggf4', 'MCstat_Htozg_mu_cat_ggf4', 'MCstat_Htomm_cat_ggf4', 'MCstat_Htozg_el_cat_ggf3', 'MCstat_Htozg_mu_cat_ggf3', 'MCstat_Htomm_cat_ggf3', 'MCstat_Htozg_el_cat_ggf2', 'MCstat_Htozg_mu_cat_ggf2', 'MCstat_Htomm_cat_ggf2', 'MCstat_Htozg_el_cat_ggf1', 'MCstat_Htozg_mu_cat_ggf1', 'MCstat_Htomm_cat_ggf1', 'MCstat_Htozg_el_cat_vbf4', 'MCstat_Htozg_mu_cat_vbf4', 'MCstat_Htomm_cat_vbf4', 'MCstat_Htozg_el_cat_vbf3', 'MCstat_Htozg_mu_cat_vbf3', 'MCstat_Htomm_cat_vbf3', 'MCstat_Htozg_el_cat_vbf2', 'MCstat_Htozg_mu_cat_vbf2', 'MCstat_Htomm_cat_vbf2', 'MCstat_Htozg_el_cat_vbf1', 'MCstat_Htozg_mu_cat_vbf1', 'MCstat_Htomm_cat_vbf1', 'MCstat_Htozg_el_cat_vh3l', 'MCstat_Htozg_mu_cat_vh3l', 'MCstat_Htomm_cat_vh3l', 'MCstat_Htozg_el_cat_vhmet', 'MCstat_Htozg_mu_cat_vhmet', 'MCstat_Htomm_cat_vhmet', 'MCstat_Htozg_el_cat_tthhad', 'MCstat_Htozg_mu_cat_tthhad', 'MCstat_Htomm_cat_tthhad', 'MCstat_Htozg_el_cat_tthlep', 'MCstat_Htozg_mu_cat_tthlep', 'MCstat_Htomm_cat_tthlep', 'MCstat_Htozg_el_cat_untagged', 'MCstat_Htozg_mu_cat_untagged', 'MCstat_Htomm_cat_untagged']

DATACARD_NAME = 'datacards/hzg_datacard_v1p0.txt'
N_CATEGORIES = 13

systs_sorted = []

with open(DATACARD_NAME, 'r') as datacard:
  datacard_content = datacard.read().split('\n')
  for line in datacard_content:
    line_split = line.split()
    if len(line_split)>0:
      if (line_split[0] in SYSTS) and not ('param' in line):
        print(f'{line_split[0]}: ',end='')
        syst_values = []
        for i in range(N_CATEGORIES):
          for j in range(2):
            syst_value = line_split[i*4+j+2]
            #syst_value = line_split[i*3+j+2] #no hto mumu
            if '/' in syst_value:
              syst_value_up = float(syst_value.split('/')[0])
              syst_value_dn = float(syst_value.split('/')[1])
              if (syst_value_up < 1.0 and syst_value_up > 0.0):
                syst_value_up = 1.0/syst_value_up
              elif syst_value_up == 0.0:
                syst_value_up = 9.99
              if (syst_value_dn < 1.0 and syst_value_dn > 0.0):
                syst_value_dn = 1.0/syst_value_dn
              elif syst_value_dn == 0.0:
                syst_value_dn = 9.99
              syst_values.append(syst_value_up)
              syst_values.append(syst_value_dn)
            else:
              if syst_value == '-':
                syst_value = 1.0
              else:
                syst_value = float(syst_value)
              if (syst_value < 1.0):
                syst_value = 1.0/syst_value
              syst_values.append(syst_value)
        print(f'{min(syst_values):.4f} {max(syst_values):.4f} ',end='')
        print(f'{mean(syst_values):.4f} {stdev(syst_values):.4f} ')
        systs_sorted.append((line_split[0], max(syst_values)))

systs_sorted = sorted(systs_sorted, key=lambda a: a[1])
print_sorted = True
if (print_sorted):
  for syst in systs_sorted:
    print(f'{syst[0]}: {syst[1]}')
