import numpy as np
import matplotlib.pyplot as plt

input_file = "AN_tables_Sep2025/cutflow_resolved/Cutflow_resolved_4b_sigbkgmc_run2_cutflow_resolved_deepflavb_2dxsec_ptmiss_trig_ttjstitch_dyweight_wjstitch_lumi_1p0.tex"

categories = ['Others', 'QCD', 'Single Top', 'TT+X', 'Z+Jets', 'W+Jets']
colors = ["#717581", "#B9AC70", "#A96B59", "#832DB6", "#E76300", "#3F90DA"]

def generate_pie(file_path, cat, col, title, filename, line_num=18):
    with open(file_path, 'r') as f:
	for i, line in enumerate(f, 1):
	    if (i == line_num):
		print(line)
		total = float(line.split("&")[len(cat)+1])/100
		comp = []
		for j in range(len(cat)):
		    comp.append(float(line.split("&")[j+1])/total)
    plt.pie(comp, colors=col, autopct='%1.1f%%')
    plt.title(title)
    plt.legend(cat, fontsize='small')
    plt.axis('equal')
    plt.savefig(filename)

generate_pie(input_file, categories, colors, "Event composition in fitting SR", "AN_tables_Sep2025/resolved_fittingSR_eventcomp.pdf", line_num=30)

