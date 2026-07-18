#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import matplotlib.lines as mlines
import mplhep as hep

import SigStrUtilities as SSU

print('log 0')

xlow = -15
xhigh = (23)//1 #27.5

#Set x-axis data and uncertainties, values ending in 0.50 are being rounded down, so I added 0.01 to fix this issue
xAxisOffset = 0
x = SSU.x 
x = [sigStr - xAxisOffset for sigStr in x]
xTotErr_ = [SSU.xTotErrLow,SSU.xTotErrHigh]
xSystErr_ = [SSU.xSystErrLow,SSU.xSystErrHigh]
xStatErr_ = [SSU.xStatErrLow,SSU.xStatErrHigh]

y = SSU.yAxisDistance(yAxisSpacing=3.3,offset=0)#, offset=0)
ylow = 0
yhigh = (max(y)*(7.15/6))//1 #(max(y)*(7.15/6))//1 

plot_width = 14.3 #12.8
plot_height = 8.4 #6.4
fig, ax = plt.subplots(figsize=(plot_width,plot_height))#9.64,6.4))
plt.tight_layout(pad=1.0, w_pad=1.0, h_pad=2.0)

print('log 1')

patch_common_settings = {'linewidth' : 1}
stat_patch_settings = {'facecolor':'#74bdfc', 'edgecolor':'#1389f0', 'label':'Stat.', **patch_common_settings}
syst_patch_settings = {'facecolor':'#ffda69', 'edgecolor':'#f2ba0f', 'label':'Syst.', **patch_common_settings}

syst_patch = mpatches.Patch(facecolor='#ffda69', edgecolor='#f2ba0f', label='Syst. unc.', **patch_common_settings)
stat_patch = mpatches.Patch(facecolor='#74bdfc', edgecolor='#1389f0', label='Stat. unc.', **patch_common_settings)
#total_unc  = mlines.Line2D([], [], color='black', linestyle='solid', label='Total unc.')
SMLine = mlines.Line2D([], [], color='black', linestyle='dashed', label='SM')
CombfitLine = mlines.Line2D([], [], color='magenta', linestyle='solid', label='Simultaneous fit') #, fontweight='bold')

plt.xlim(xlow, xhigh)
plt.ylim(ylow, yhigh)

plt.style.use([hep.style.ROOT, hep.style.firamath])
hep.cms.text("",fontsize=14)#Preliminary", fontsize=12)
#hep.cms.label(loc=0)

def place_rectangle(ax,x,xunc,y,ywidth,settings):
  yll = y-ywidth
  xll = x-xunc[0]
  xwidth = xunc[0]+xunc[1]
  rect = mpatches.Rectangle((xll,yll),xwidth,ywidth,**settings)
  ax.add_patch(rect)

def add_uncertainty_bars(ax, x, xunc, y, ywidth,settings):
  for idx in range(len(x)):
    place_rectangle(ax, x[idx], [xunc[0][idx], xunc[1][idx]], y[idx]+1, ywidth, settings)


linewidth = 2
edgewidth = 0.4
uppery = [yid + linewidth/2 for yid in y]
lowery = [yid - linewidth/2 for yid in y]

add_uncertainty_bars(ax, x, xStatErr_, uppery, linewidth, stat_patch_settings)
add_uncertainty_bars(ax, x, xSystErr_, lowery, linewidth, syst_patch_settings)
errB = plt.errorbar(x, y, xerr=xTotErr_,  fmt='sk', color='black', marker='o',markersize = 3, capsize=2*linewidth+edgewidth, label='Total unc.')
fig.subplots_adjust(left=0.16,right=0.98,top=0.94,bottom=0.1)#left=0.192,right=0.89,top=0.94,bottom=0.1)

common_settings = {'clip_on' : False, 'annotation_clip' : False}# {'bbox' : dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0)}
textcolor = 'black'

#Attempt to create scaling for horizontal text locations
x_scale = plot_width/9.64 

#horizontal offsets for annotations

label_x = 15.75 #15.75
prec_scale = 0.3/x_scale
pech = 0.6/x_scale
unc_offset = 2.0/x_scale 
colgap = 3.0/x_scale #3.5/x_scale
centering1 = 1.25/x_scale 
centering2 = -0.95/x_scale 
psign_offset = 0.42/x_scale
nsign_offset = 0.38/x_scale

#Vertical offsets for annotations
label_y_offset = 2.5 #0.5
label_y_unc_dist = 2.25#2.00
label_y_unc_offset = 0.8
yct_offset = 8.0

txt_size = 15 #12
unc_txt_size = 12 #8

#Adds uncertainty to the plot
def add_uncertainty(unc, xyplace, textsize=8, prec = 2, common_settings={}):
  plt.annotate('-',                         xy=(xyplace[0] - nsign_offset, xyplace[1]+label_y_unc_offset-label_y_unc_dist), size = textsize, **common_settings) 
  plt.annotate(("%." + str(prec) + "f") % round(unc[0],prec), xy=(xyplace[0],                xyplace[1]+label_y_unc_offset-label_y_unc_dist), size = textsize, **common_settings)
  
  plt.annotate('+',                     xy=(xyplace[0] - psign_offset, xyplace[1]+label_y_unc_offset+label_y_unc_dist), size = textsize, **common_settings)
  plt.annotate(("%." + str(prec) + "f") % round(unc[1],prec), xy=(xyplace[0],                xyplace[1]+label_y_unc_offset+label_y_unc_dist), size = textsize, **common_settings)

#Adds signal strength with uncertainty to the plot
def add_sigstr(xloc, yloc, sigstr, totunc, prec, common_settings={}):
  #if sigstr < 0:
  #  plt.annotate('-', xy=(xloc-psign_offset,yloc), size=12)
  
  pech = prec_scale*prec if prec > 1 else 0;

  sigstr_s = ("%." + str(prec) + "f") % round(sigstr,prec)  

  if sigstr < 0:
    plt.annotate(sigstr_s, xy=(xloc-psign_offset,yloc), size = txt_size, **common_settings)
  else: 
    plt.annotate(sigstr_s, xy=(xloc-pech/2,yloc), size = txt_size, **common_settings)
 
  add_uncertainty(unc=totunc, xyplace=(xloc+unc_offset+pech/2, yloc), textsize=unc_txt_size, prec=prec, common_settings=common_settings)  

#Adds signal strength, systematic uncertainty, and stat uncertainty to the plot.
def add_category_sigstr(xloc, yloc, sigstr, totunc, systunc, statunc, prec, common_settings={}):  
  #Add the signal strength to the plot
  add_sigstr(xloc, yloc,sigstr,totunc, prec, common_settings=common_settings)
  
  #Update the location to include padding for the sig str and the uncertainty
  xloc = xloc + unc_offset
  add_uncertainty(unc = statunc, xyplace = (xloc+colgap, yloc), textsize=unc_txt_size, prec=prec, common_settings=common_settings)
  add_uncertainty(unc = systunc, xyplace = (xloc+2*colgap,   yloc), textsize=unc_txt_size, prec=prec, common_settings=common_settings)

#Add column titles
plt.annotate('$\mu$',  xy=(label_x+centering1,                    y[0] + yct_offset), size = txt_size, fontweight='bold', fontstyle='italic', **common_settings) #weight='bold',  
plt.annotate('Stat.', xy=(label_x+centering2+unc_offset+colgap,   y[0] + yct_offset), size = txt_size, weight='bold', **common_settings) 
plt.annotate('Syst.', xy=(label_x+centering2+unc_offset+2*colgap, y[0] + yct_offset), size = txt_size, weight='bold', **common_settings) 
#Print each row
for i in range(len(x)):
  #Bold combined fit columns
  if i==13:
    common_settings['weight'] = 'bold'
    add_category_sigstr(label_x, y[i]-label_y_offset, x[i], [xTotErr_[0][i],xTotErr_[1][i]], [xSystErr_[0][i],xSystErr_[1][i]], [xStatErr_[0][i],xStatErr_[1][i]], 2, common_settings)
  else:
    add_category_sigstr(label_x, y[i]-label_y_offset, x[i], [xTotErr_[0][i],xTotErr_[1][i]], [xSystErr_[0][i],xSystErr_[1][i]], [xStatErr_[0][i],xStatErr_[1][i]], 1, common_settings)
  
  #Below line can help test alignment
  #plt.hlines(y[i], xlow, xhigh, colors='#9c9898').set_linewidth(1.0)

print('log 2')

#Set Y axis ticks
plt.yticks(y, SSU.categories, size=20.)
ticks = ax.get_yticklabels()
ticks[13].set_fontweight('bold')
plt.tick_params(axis="y",direction="in", length=3)

#Set X axis ticks
plt.tick_params(axis="x", direction="in", length=3)
plt.xticks(np.arange(xlow, xhigh, step=5),size=16)
ax.xaxis.set_minor_locator(MultipleLocator(1))

#Add lines for SM
yoffset = 1
ybounds = [ylow+yoffset, ylow-yoffset]
bound_vline = yhigh - 12 #24
plt.vlines(1, ylow, bound_vline, colors='k', linestyles ='dashed', dashes=(0, (5, 10))).set_linewidth(0.5)
plt.vlines(15, ylow, bound_vline, colors='#9c9898').set_linewidth(1.0)
plt.vlines(x[13], ylow, bound_vline, colors='magenta').set_linewidth(0.5)
plt.hlines(bound_vline, xlow, xhigh, colors='#9c9898').set_linewidth(1.0)



#fullFitLine, = ax.plot([x[13],x[13]], ybounds, color='magenta', label='Combined fit', linewidth=0)
#Create Legend
legend_handles = [stat_patch, syst_patch, errB, SMLine, CombfitLine]
#legend_settings = { 'loc' : 'upper right', 'prop' : {'size': 12}, 'ncol' : 3, 'columnspacing' : 0.7, 'handletextpad' : 0.4, 'borderaxespad' : 0.25}
legend_settings = { 'loc' : 'upper right', 'prop' : {'size': 14}, 'ncol' : 5, 'columnspacing' : 0.7, 'handletextpad' : 0.4, 'borderaxespad' : 0.25}
plt.legend(handles=legend_handles,**legend_settings)

lumi_settings = {'size' : 14, 'fontname' : "sans serif"} #old size 12
#plt.annotate('138 fb$^{-1}$ ($\sqrt{s}$ = 13 TeV),', xy=(xlow+0.75, yhigh-9), **lumi_settings)
#plt.annotate('62 fb$^{-1}$ ($\sqrt{s}$ = 13.6 TeV)', xy=(xlow+0.75, yhigh-19), **lumi_settings)

#plt.annotate('138 fb$^{-1}$ ($\sqrt{s}$ = 13 TeV), 62 fb$^{-1}$ ($\sqrt{s}$ = 13.6 TeV)', xy=(xlow+0.75, yhigh-10), **lumi_settings)
plt.annotate('138 fb$^{-1}$ ($\sqrt{s}$ = 13 TeV), 62 fb$^{-1}$ ($\sqrt{s}$ = 13.6 TeV)', xy=(xlow+0.2, yhigh-8.0), **lumi_settings)


#plt.rcParams["mathtext.default"] = "it"
#plt.rc('text', usetex = True)
plt.rcParams['text.usetex'] = True

plt.rcParams['mathtext.fontset'] = 'cm'
plt.xlabel(r'$\mathit{\mu}$', size=26) #12


print("at end")
#plt.savefig('./signal-strength-plot.png', dpi=600)
plt.savefig('./signal-strength-plot.pdf', dpi=1000)






#plt.annotate('SM', xy=(2, 3),size=12)


#for i in range(len(x)):
  #plt.annotate(str(x[i])+'$\pm$'+str(xerr1_[i]), xy=(x[i], y[i]), xytext=(8, y[i]), size = 12)
#  if i == 13:
#    plt.annotate('%.2f'%x[i], xy=(x[i], y[i]), xytext=(label_x, y[i]+label_y_offset), size = 12, color='magenta', facecolor='white') 
#  else:
#    plt.annotate('%.2f'%x[i], xy=(x[i], y[i]), xytext=(label_x, y[i]+label_y_offset), size = 12, facecolor='white')
#  plt.annotate('- ' + str(xerr1_[0][i]), xy=(x[i], y[i]), xytext=(label_x + label_x_offset, y[i]+0.4), size = 8)
#  plt.annotate('+' + str(xerr1_[1][i]), xy=(x[i], y[i]), xytext=(label_x + label_x_offset, y[i]+3.4), size = 8)


