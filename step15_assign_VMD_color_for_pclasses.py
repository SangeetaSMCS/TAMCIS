import numpy as np
import json

input_json_file = "input_TAMCIS.json"
with open(input_json_file,"r") as f:
    dict_TAMCIS = json.load(f)

peptide = dict_TAMCIS["peptide"]
num_frame = dict_TAMCIS["num_frame"]
tag = dict_TAMCIS["tag"]
input_file = dict_TAMCIS["cluster_alphanum_file"]
categories = dict_TAMCIS["prominent_class"]
len_promt_class = len(categories)


color_pclass_in_VMD = {}
for pkey,pvalue in categories.items():
    print("ID =",pvalue[2],"PClass name =",pkey,"Color assigned = ",pvalue[1])
    #print("Provide equivalent color ID in VMD for ",pvalue[1])
    VMD_cid = int(input("Provide equivalent color ID in VMD = "))
    color_pclass_in_VMD[pvalue[2]] = [VMD_cid,pkey]
json_file = "VMD_color_prominent_class.json"
with open(json_file,"w+") as f:
    json.dump(color_pclass_in_VMD,f,indent=4)

tcl_file = "number_prominent_class.tcl"
with open(tcl_file,"w+") as f:
    print("set num_pclass "+str(len(categories)),file=f)
