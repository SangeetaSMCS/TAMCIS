import numpy as np
import json
import matplotlib.pyplot  as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn import metrics
import time,os




#--------------------------------------------------------------------------------------------------------------------
#		UDER INPUT
#-------------------------------------------------------------------------------------------------------------------
input_json_file = "input_TAMCIS.json"
with open(input_json_file,"r") as f:
    dict_TAMCIS = json.load(f)

#--------------------------------------------------------------------------------------------------------------------
start = time.time()

peptide = dict_TAMCIS["peptide"]
xyz_filename = dict_TAMCIS["xyz_binned_file"]
time_slice = dict_TAMCIS["num_frame"]


# time gap for clustering
time_gap = 1

#----------------------------------------------------------------------------------------------------------------------
#		For the time being DO NOT TOUCH
#----------------------------------------------------------------------------------------------------------------------


data = np.genfromtxt(xyz_filename,skip_header = 3)
data_str = np.genfromtxt(xyz_filename,skip_header = 3,dtype=str)
with open(xyz_filename) as f:
    data_linewise = f.readlines()

#original file name
org_filename = data_linewise[0].replace("\n","")
# Get the parameter list
parameter_list = data_linewise[1].split(' ')
print(parameter_list)


labels = parameter_list[1:-1]
print(labels)

segname = data_str[:,0]

time_slice_list = [i for i in range(0,time_slice,time_gap)]

data_slice = []
data_slice_str = []
index_list = []
for idx in range(len(segname)):
    t = data[idx,-1]
    if(t in time_slice_list):
        data_slice.append(data[idx])
        data_slice_str.append(data_str[idx])
        index_list.append(idx)
        #print(time,data_slice_str[0][0])

data_slice = np.array(data_slice)
data_slice_str = np.array(data_slice_str)
  
segname_slice = data_slice_str[:,0]
print(segname_slice)

# number of hyper parameters
hyperparam_len = len(labels)

# loop over min_sample
for i in range(1,2):
    ID =segname
    save_path ="."
    
    #cluster=clt.meanshift_h(quan,data,ID,x,y,z,save_path,fiber_name,num_one_molecule_representation)
    
    #min_sample = i*50
    min_sample = i*1

    model = DBSCAN(eps=0.5, min_samples=min_sample)
    pred = model.fit_predict(data_slice[:,1:hyperparam_len+1])


    print("prediction complete")

    dict_cluster = {}
    dict_cluster["pred"] = pred.tolist()

    xy_name = ""
    for lab in labels:
        lab = lab.replace("pep_H_","")
        lab = lab.replace("_contact","")
        lab = lab.replace("end_end","ee")
        lab = lab.replace("end_mid","em")
        lab = lab.replace("mid_mid","mm")
        lab = lab.replace("sidechain","sc")
        lab = lab.replace("electrostatic","electro")
        xy_name = xy_name + lab +"_"
    
    file_store = peptide+"_cluster_data_"+xy_name+"min_sample_"+str(min_sample)+"_without_time_tgap_"+str(time_gap)+"_BINNED.dat"
    f = open(file_store,"w+")
    print(org_filename,file=f)
    print("idx","org_idx","c_id",file=f)
    for i in range(len(pred)):
        print(i,index_list[i],pred[i],file=f)

    f.close()

    end = time.time()-start
    print("PATH : ",file_store)
    print("time taken=",end)

dict_TAMCIS["cluster_file"] = file_store
with open(input_json_file,"w+") as f:
    dict_TAMCIS = json.dump(dict_TAMCIS,f,indent=4)
