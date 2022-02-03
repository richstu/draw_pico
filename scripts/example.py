import drawpico
import ROOT

#set plot style
linear_plot = [drawpico.PlotOpt('txt/plot_styles.txt','CMSPaper2022')]

#assign data samples (processes)
print('assigning data samples')
procs = drawpico.SampleLoader().LoadNamedFunc("triggers_data",drawpico.met_trigger).LoadSamples('txt/samples.txt','HHMetPaper')

print('initializing PlotMaker')
pm = drawpico.PlotMaker()

print('making plot settings')
hig_cand_am_nf = drawpico.NamedFunc('hig_cand_am[0]')
print('making plot settings')
plot_axis = drawpico.Axis(20,0.0,200.0,hig_cand_am_nf,'<m_{bb}> [GeV]',[],[100.0,140.0])
print('making plot settings')
plot_selection = drawpico.final_pass_filters & 'met/mht<2 && met/met_calo<2&&weight<1.5&&ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&hig_cand_drmax[0]<=2.2&&hig_cand_dm[0]<=40&&nbt>=2&&nbm>=3'
print('making plot settings')
weight = drawpico.final_weight

print('pushing histogram')
pm.PushHist1D(plot_axis, plot_selection, procs, linear_plot).Weight(weight).Tag(
    'FixName:amjj_test').LuminosityTag('137.4')

pm.MakePlots(1.0)
