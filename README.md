draw_pico
========

**draw_pico** is a user-level analysis framework for making 1-d histograms, 2-d histograms, efficiency plots, tables, events scans, and statistical models (in the form of [combine](https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit) datacards) from high energy physics n-tuples. Currently supported n-tuple formats include [picos](https://github.com/richstu/nano2pico) and [NanoAOD](https://github.com/cms-sw/cmssw/tree/master/PhysicsTools/NanoAOD), though new formats can be easily added with a specification text file. It has been used for the the CMS 1 lepton+jets+MET (MJ) SUSY analysis, the HH(4b)+MET analysis, and the CMS H to Z gamma analysis.

Repository inherited from: [https://github.com/ald77/ra4_draw](https://github.com/ald77/ra4_draw)

For examples of functionality see: [src/core/test.cxx](src/core/test.cxx)

## Setup

Prerequisites: (Note, on UCSB servers, `source set_env.sh` can be used to fulfill these prerequisites, which is described below in the UCSB server section.)
- ROOT 
- (Optional) Scons for building the library. Can also use the `compile.py` script instead.

Please run the below command to clone the repository

~~~~bash
git clone git@github.com:richstu/draw_pico.git
cd draw_pico
~~~~

Please run the below command to build the libraries.

~~~~bash
./compile.py
~~~~

or

~~~~bash
scons
~~~~

### Setup on UCSB servers

If you use the USCB servers,
you can use the below commands for the setup.
This also sets up the environment for the batch system at UCSB

~~~~bash
git clone --recurse-submodules git@github.com:richstu/draw_pico.git
cd draw_pico
source set_env.sh
./compile.py or scons
~~~~

An alternative method is shown below:

~~~~bash
git clone git@github.com:richstu/draw_pico.git
cd draw_pico
git submodule update --init --remote --recursive
source set_env.sh
./compile.py or scons
~~~~

## Higgsino useful commands

### To make datacards and get limits:

To plot overlap between the boosted and resolved analysis:

~~~~bash
./run/higgsino/plot_regions.exe #boosted overlap in resolved regions
./run/higgsino/plot_regions.exe --boo #resolved overlap in boosted regions
~~~~

Write all the datacards:

~~~~bash
./run/higgsino/write_kappa_datacards_priority.exe -o datacards/tchihh_onedim/ -m CN -1 -r 1 --unblind --unblind_signalregion
./run/higgsino/write_kappa_datacards_priority.exe -o datacards/tchihh_twodim/ -m N1N2 -2 -r 1 --unblind --unblind_signalregion
./run/higgsino/write_kappa_datacards_priority.exe -o datacards/t5hh_onedim/ -m T5HH -f -1 -r 1 --unblind --unblind_signalregion
./run/higgsino/write_kappa_datacards_priority.exe -o datacards/t5hh_twodim/ -m T5HH -2 -r 1 --unblind --unblind_signalregion
~~~~

For comparison of signal yields to other tables, remember to use `--recomet` such that the yield is not averaged with the one obtained using GenMET. To make a datacard for the boosted case, use `-t boosted`. Use option `--unblind` to include data. To run on a particular point add, e.g. `-p "700_1"`. For computers with smaller amounts of memory, there may not be sufficient memory to generate all the 2D T5HH datacards at once. In this case, you can split the job in half with the `-s` flag, which sets a lower/upper bound on the gluino mass. 

~~~~bash
./run/higgsino/write_kappa_datacards_priority.exe -o datacards/t5hh_twodim/ -m T5HH -2 -r 1 --unblind --unblind_signalregion -s "-2100"
./run/higgsino/write_kappa_datacards_priority.exe -o datacards/t5hh_twodim/ -m T5HH -2 -r 1 --unblind --unblind_signalregion -s 2100
~~~~

Then to get a limit interactively:

~~~~bash
./run/higgsino/scan_point.exe -f datacards/tchihh_onedim/datacard-TChiHH_mChi-700_mLSP-0_Tune_2016,2017,2018_priority1_resolved.txt
~~~~

To get limits for the full scan in the batch, use the following script to generate batch commands.

~~~~bash
./scripts/write_combine_cmds.py --card_dir datacards/tchihh_onedim -m CN
./scripts/write_combine_cmds.py --card_dir datacards/tchihh_twodim -m N1N2
./scripts/write_combine_cmds.py --card_dir datacards/t5hh_onedim -m T5HH
./scripts/write_combine_cmds.py --card_dir datacards/t5hh_twodim -m T5HH
~~~~

These can be run locally with `./scripts/run_commands.py` or submitted to the UCSB batch system with `auto_submit_jobs.py`. To combine outputs and make the limit plot

~~~~bash
cat datacards/tchihh_onedim/scan_point*/limit*txt | sort >> tchihh_onedim_resolved_limits.txt
cat datacards/tchihh_twodim/scan_point*/limit*txt | sort >> tchihh_twodim_resolved_limits.txt
cat datacards/t5hh_onedim/scan_point*/limit*txt | sort >> t5hh_onedim_resolved_limits.txt
cat datacards/t5hh_twodim/scan_point*/limit*txt | sort >> t5hh_twodim_resolved_limits.txt
./run/higgsino/plot_limit.exe -f tchihh_onedim_resolved_limits.txt --drawData -t tchihh_onedim_resolved
./run/higgsino/limit_scan.exe -f tchihh_twodim_resolved_limits.txt -m N1N2 -t tchihh_twodim_resolved --unblind
./run/higgsino/plot_limit.exe -f t5hh_onedim_resolved_limits.txt --drawData -m T5HH -t t5hh_onedim_resolved
./run/higgsino/limit_scan.exe -f t5hh_twodim_resolved_limits.txt -m T5HH -t t5hh_twodim_resolved --unblind
~~~~

### To generate AN tables/plots:

~~~~bash
./run/higgsino/an_plot_triggers.exe --unblind --year 2016 --string_options systematic,efficiency,cr
./run/higgsino/an_plot_triggers.exe --unblind --year 2017 --string_options systematic,efficiency,cr
./run/higgsino/an_plot_triggers.exe --unblind --year 2018 --string_options systematic,efficiency,cr
./run/higgsino/an_plot_triggers.exe --unblind --year run2
./run/higgsino/an_plot_selection.exe --unblind --year run2
./run/higgsino/an_plot_backgroundEstimation.exe --unblind --year run2
./run/higgsino/plot_kappas.exe -s search --scen mc -y run2 -o paper_style
./run/higgsino/an_plot_syst_ttbar.exe --unblind --year run2
./run/higgsino/an_plot_syst_zll.exe --unblind --year run2
./run/higgsino/an_plot_syst_qcd.exe --unblind --year run2
./scripts/datacard_syst_utils.py
./run/higgsino/plot_n_minus_1.exe --unblind --year run2 --sample search
./run/higgsino/plot_search_unblind.exe -u -a -o plot_baseline,paper_style,plot_in_btags,plot_in_btags_with_met_split
./run/higgsino/plot_kappas.exe -s search --scen data -y run2 -o paper_style,use_datacard_results,do_zbi --unblind
./scripts/plot_ra4style_results.py
~~~~

### To generate paper plots:

~~~~bash
./run/higgsino/plot_search_unblind.exe -u -a -o plot_baseline,paper_style,plot_in_btags,plot_in_btags_with_met_split
./run/higgsino/plot_kappas.exe -s search --scen mc -y run2 -o paper_style
./run/higgsino/plot_kappas.exe -s search --scen data -y run2 -o paper_style,use_datacard_results,do_zbi --unblind
./scripts/plot_ra4style_results.py
~~~~

### To generate supplementary plots:

~~~~bash
./run/higgsino/an_plot_triggers.exe -u -o plot,paper_style
./run/higgsino/an_plot_syst_zll.exe --year run2 --unblind -o paper_style
./run/higgsino/an_plot_syst_qcd.exe --year run2 --unblind -o paper_style
./run/higgsino/an_plot_syst_ttbar.exe --year run2 --unblind -o plot_data_vs_mc,paper_style
./run/higgsino/plot_kappas.exe --sample ttbar --scen data --unblind --year run2 -o paper_style
./run/higgsino/plot_kappas.exe --sample zll --scen data --unblind --year run2 -o paper_style
./run/higgsino/plot_kappas.exe --sample qcd --scen data --unblind --year run2 -o paper_style
./run/higgsino/plot_search_unblind.exe --year run2 -u --unblind_signal -o plot_in_btags,plot_in_btags_with_met_split,paper_style,supplementary
./scripts/plot_signal_efficiency.py
./run/higgsino/plot_phase_space.exe
./run/higgsino/plot_supplementary.exe -y run2 -o cutflow
./run/higgsino/plot_supplementary.exe -y run2 -o pie
./scripts/plot_ra4style_results.py --signal
~~~~

## Getting Higgsino cross-section (CN to N1N2) scale factors.

### Step 1. Download cross-section files

~~~~bash
scp -r lxplus:/afs/cern.ch/user/a/amete/public/EWKGauginoCrossSections_13TeV cross_section
~~~~

### Step 2. Fix bug in script

Set masses, xsecs, xsecUncs to 0 when initializing.

~~~~bash
cd cross_section
sed -i 's/std::vector<double>\* masses;/std::vector<double>\* masses=0;/' get_gaugino.C
sed -i 's/std::vector<double>\* xsecs;/std::vector<double>\* xsecs=0;/' get_gaugino.C
sed -i 's/std::vector<double>\* xsecUncs;/std::vector<double>\* xsecUncs=0;/' get_gaugino.C
~~~~

### Step 3. Run script for all mass points.

 model "CN" (mixing) or "N1N2" (no mixing).

~~~~bash
cd cross_section
../scripts/print_higgsino_cross_sections.py -i /net/cms25/cms25r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D -m CN
../scripts/print_higgsino_cross_sections.py -i /net/cms25/cms25r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D -m N1N2
../scripts/print_higgsino_cross_sections.py -i /net/cms25/cms25r0/pico/NanoAODv7/nano/2016/SMS-TChiHH_2D -c
~~~~

## Zgamma useful commands

