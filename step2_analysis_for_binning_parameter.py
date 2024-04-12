import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import pandas as pd
import general_code as gc

input_json_file = "input_TAMCIS.json"
with open(input_json_file,"r") as f:
    dict_TAMCIS = json.load(f)

#--------------------------------------------------------------------------------------------------------------------
peptide = dict_TAMCIS["peptide"]
xyz_filename = dict_TAMCIS["xyz_raw_file"]

data = np.genfromtxt(xyz_filename,skip_header = 2)
with open(xyz_filename) as f:
    data_linewise = f.readlines()

parameter_list = dict_TAMCIS["parameters_list"]
print(parameter_list)

col_list = ['g','orange','b','gray']
i = 0

fig,ax = plt.subplots(1,4,figsize=(16,3))
for param_idx in range(1,len(parameter_list)+1):
    param_name = parameter_list[param_idx-1]
    x = data[:,param_idx]
    time = data[:,-1]
    #fig,ax = plt.subplots(figsize=(12,6.2),ncols=3,gridspec_kw={'width_ratios': [3, 3, 0.2]})
    
    #ax[i].hist(x,bins=20,color = col_list[i],alpha=0.5)
    ax[i].hist(x,bins=20,color = col_list[i], weights=np.ones(len(x)) / len(x))

    ax[i].yaxis.set_major_formatter(PercentFormatter(1))
    xlabel = param_name.replace("_"," ")
    xlabel = xlabel.replace("pep","")
    xlabel = xlabel.replace("H","HB :")
    xlabel = xlabel.replace("sidechain","sidechain-peptide")

    ax[i].set_xlabel(xlabel)
    if(i == 0):
        ax[i].set_ylabel("frequency")
    #formatter = FuncFormatter(to_percent)
    #plt.gca().yaxis.set_major_formatter(formatter)
    i = i + 1
    #h = ax[1].hist2d(time,x)
    #ax[1].set_ylabel(param_name)
    #ax[1].set_xlabel("time")
    #fig.colorbar(h[3],ax[2])

fig.suptitle(peptide.upper())
plt.savefig(peptide+"_analysis_4_bining_chbsc.png",dpi=300,bbox_inches='tight')
plt.show()
