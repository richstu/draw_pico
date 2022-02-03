import drawpico

#set plot style
linear_plot = [drawpico.PlotOpt('txt/plot_styles.txt','CMSPaper2022')]

#assign data samples (processes)
procs = drawpico.SampleLoader().LoadNamedFunc("triggers_data",drawpico.met_trigger).LoadSamples('txt/samples.txt','HHMetPaper')

pm = drawpico.PlotMaker()

#define plots
pm.PushHist1D(
    drawpico.Axis(20,0.0,200.0,'hig_cand_am[0]','<m_{bb}> [GeV]',[],[100.0,140.0]),
    drawpico.final_pass_filters & 'met/mht<2 && met/met_calo<2&&weight<1.5&&ntk==0&&!low_dphi_met&&nvlep==0&&met>150&&njet>=4&&njet<=5&&hig_cand_drmax[0]<=2.2&&hig_cand_dm[0]<=40&&nbt>=2&&nbm>=3',
    procs, linear_plot).Weight(drawpico.final_weight).Tag(
    'FixName:amjj_test').LuminosityTag('137.4')

#make plots
pm.MakePlots(1.0)

