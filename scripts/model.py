'''
https://faroit.com/keras-docs/1.2.0/backend/
'''
from collections import OrderedDict
from six.moves import range
import numpy as np
import os
import sys
import h5py

from keras.layers.core import Dense, Dropout, Flatten
from keras.layers import Embedding, BatchNormalization, CuDNNLSTM, LSTM, Bidirectional, Input, \
    Concatenate, Multiply, Dot, Reshape, Activation, Lambda, Masking,concatenate,Add,Conv2D,ConvLSTM2D
from keras.layers.convolutional import Convolution1D, MaxPooling1D,AveragePooling1D
from keras.models import Model
from keras.callbacks import EarlyStopping, ModelCheckpoint, Callback, TensorBoard
import tensorflow as tf
from keras import backend as K
from keras.initializers import random_normal
from keras import losses,layers,optimizers,regularizers
from keras.engine.topology import Layer
from tensorflow.python.framework import ops

from scipy.special import softmax
from multihead_attention import Attention
from attention_decoder import *

OUTPATH = None
# from sklearn.utils import shuffle

np.random.seed(1234)
encoding_baseout = OrderedDict([
    (1,'A'),
    (2,'C'),
    (3,'G'),
    (4,'T'),
    (0,'N'),

])

def generator_all(IVT_pos_file,IVT_neg_file,real_candidate_files,real_noncandidate_files,KO_candidate_files,KO_noncandidate_files,batch_size,bin=0.2,binsample=200,WT_label=1,
        binselecttype=1,sample_weights_flag=True,featuredim=3,windowsize=9,weight_exp=1,KO_noncandidate_weight=1,real_noncandidate_weight=1): #fixed start np.arnage(0,1,bin)
        print("generator_all")
        real_candidate_files=real_candidate_files.replace("\"",'')
        real_noncandidate_files=real_noncandidate_files.replace("\"",'')
        KO_candidate_files=KO_candidate_files.replace("\"",'')
        KO_noncandidate_files=KO_noncandidate_files.replace("\"",'')
        
        #IVT data
        IVT_pos_h5fs=h5py.File(IVT_pos_file,"r")
        IVT_neg_h5fs=h5py.File(IVT_neg_file,"r")
        samplesize_TN = len(IVT_neg_h5fs['X'])
        samplesize_TP = len(IVT_pos_h5fs['X'])
        index_IVT_TP = 0
        index_IVT_TN = 0
        
        #KO noncandidate files start
        KO_noncandidate_h5fs,KO_noncandidate_samplesize,index_noncandidate_KO=[],[],[]
        index=0
        for KO_noncandidate in KO_noncandidate_files.split(","):
            if os.path.isfile(KO_noncandidate):
                temp = h5py.File(KO_noncandidate,"r")
                if len(temp['X'])!=0:
                       KO_noncandidate_h5fs.append(temp)
                       temp_size=len(KO_noncandidate_h5fs[index]['X'])
                       KO_noncandidate_samplesize.append(temp_size)
                       select_index = np.random.choice(np.arange(0,temp_size,int(temp_size/10000)))
                       index_noncandidate_KO.append(select_index)
                       print("debug "+KO_noncandidate+" index_noncandidate_KO="+str(select_index)+"\n")
                       index+=1
                else:
                       print("debug error: KO_noncandidate_samplesize="+str(temp['X']))
                       print(KO_noncandidate)
        
        batch_size_KO_noncandidate=int(batch_size*KO_noncandidate_weight/len(KO_noncandidate_h5fs))
        print("debug batch_size_KO_noncandidate="+str(batch_size_KO_noncandidate))
        print("KO_noncandidate_file_num="+str(len(KO_noncandidate_h5fs)))
        #KO noncandidate files end
        #real noncandidate files start
        real_noncandidate_h5fs,real_noncandidate_samplesize,index_noncandidate_real={},{},{}
        batch_size_real_noncandidate={}
        for datatype in range(len(real_noncandidate_files.split(";"))): #use ; for different datatype
           real_noncandidate_h5fs[datatype],real_noncandidate_samplesize[datatype],index_noncandidate_real[datatype]=[],[],[]
           index=0
           for real_noncandidate in real_noncandidate_files.split(";")[datatype].split(","):
              if os.path.isfile(real_noncandidate):
                temp = h5py.File(real_noncandidate,"r")
                if len(temp['X'])==0:
                       print("Error: real_noncandidate_h5fs="+str(len(temp['X'])))
                       print(real_noncandidate)
                else:
                     real_noncandidate_h5fs[datatype].append(temp)
                     temp_size=len(real_noncandidate_h5fs[datatype][index]['X'])
                     real_noncandidate_samplesize[datatype].append(temp_size)
                     select_index = np.random.choice(np.arange(0,temp_size,int(temp_size/10000)))
                     index_noncandidate_real[datatype].append(select_index)
                     print("debug "+real_noncandidate+" index_noncandidate_real="+str(select_index)+"\n")
                     index+=1
           
           batch_size_real_noncandidate[datatype]=int(batch_size*real_noncandidate_weight/len(real_noncandidate_h5fs[datatype]))
           print("real_noncandidate_file_num="+str(len(real_noncandidate_h5fs[datatype])))
        
        #real noncandidate files end
        win_center_index=int(len(KO_noncandidate_h5fs[0]['X'][0])/2)
        nt = int(windowsize/2)
        
        AG_ratio_WT,AG_ratio_KO = [],[]
        real_candidate_h5fs=[]
        KO_candidate_h5fs=[]
        index=0
        real_candidate_file_size=[]
        for real_candidate in real_candidate_files.split(","):
          if os.path.isfile(real_candidate):
            print("openning "+real_candidate+"\n")
            temp = h5py.File(real_candidate,'r')
            file_length = len(temp['X'])
            real_candidate_file_size.append(file_length)
            if file_length==0:
                 print("Error: real candidate sample size is 0 ")
                 print(real_candidate)
                 real_candidate_file_size[index]=1 #avoid devide by 0
            else:
                AG_ratio_WT.append([])
                real_candidate_h5fs.append(temp)
                for line in real_candidate_h5fs[index]['info']:
                    try:
                         ratio=float(line.split("\t")[0].split(":")[1])
                    except:
                         print("Error:")
                         print(line)
                         exit()
                    
                    AG_ratio_WT[index].append(ratio)
                
                index+=1
        
        minsize=np.min(real_candidate_file_size)
        real_candidate_file_size=[int(np.ceil(x/minsize)) for x in real_candidate_file_size]
        print("debug real_candidate_file_size\n")
        print(real_candidate_file_size)
        real_candidate_file_num = index
        print("real_candidate_file_num="+str(real_candidate_file_num))
        index=0
        for KO_candidate in KO_candidate_files.split(","):
          if os.path.isfile(KO_candidate):
            print("openning "+KO_candidate+"\n")
            temp=h5py.File(KO_candidate,'r')
            if len(temp['X'])==0:
                 print("Error: KO candidate sample size is 0")
                 print(KO_candidate)
            else:
                KO_candidate_h5fs.append(temp)
                AG_ratio_KO.append([])
                for line in KO_candidate_h5fs[index]['info']:
                    AG_ratio_KO[index].append(float(line.split("\t")[0].split(":")[1]))
                
                index+=1
        
        KO_candidate_file_num = index
        print("KO_candidate_file_num="+str(KO_candidate_file_num))
        ratio_start = [round(x,4) for x in np.random.uniform(0.0, 1.0, binsample)]
        if binselecttype==4:
           ratio_start_index__=0
        
        ratio_bin_samplenum=np.zeros([real_candidate_file_num,len(ratio_start)])
        AG_ratio_WT_index,AG_ratio_KO_index = {},{}
        ratio_start_index=0
        for start in ratio_start:
           end = start+bin
           AG_ratio_WT_index[start]=[]
           flag_AG_ratio =0
           for _ in range(real_candidate_file_num):
                 index_start_end_WT = np.where((AG_ratio_WT[_]<=end) &(AG_ratio_WT[_]>start))[0]
                 AG_ratio_WT_index[start].append(index_start_end_WT) #AG_ratio_WT_index[start] element could be []
                 #print("len len(AG_ratio_WT_index) at AG ratio "+str(start)+"="+str(len(AG_ratio_WT_index[start][_]))+"\n")
                 ratio_bin_samplenum[_][ratio_start_index]=len(index_start_end_WT)
                 if len(index_start_end_WT) > 50:
                       flag_AG_ratio+=len(index_start_end_WT)
           
           if flag_AG_ratio==0:
                  del AG_ratio_WT_index[start]
           
           ratio_start_index+=1
        
        for start in ratio_start:
          end = start+bin
          AG_ratio_KO_index[start]=[]
          flag_AG_ratio = 0
          for _ in range(KO_candidate_file_num):
            index_start_end_KO = np.where((AG_ratio_KO[_]<=end) &(AG_ratio_KO[_]>start))[0]
            AG_ratio_KO_index[start].append(index_start_end_KO)
            #print("len len(AG_ratio_KO_index) at AG ratio "+str(start)+"="+str(len(AG_ratio_KO_index[start][_]))+"\n")
            if len(index_start_end_KO) > 50:
                  flag_AG_ratio+=len(index_start_end_KO)
          
          if flag_AG_ratio==0: #all files ==0 
                del AG_ratio_KO_index[start]
        
        for index,start in enumerate(ratio_start):
             if (start not in AG_ratio_KO_index) and (start not in AG_ratio_WT_index):
                    del ratio_start[index]
        
        #use all training samples in candidate WT and KO
        WT_select_sample_index={}
        max_WT_sample_perdatatype=np.zeros(real_candidate_file_num)
        for ratiokey in AG_ratio_WT_index:
           WT_select_sample_index[ratiokey]=[]
           for datatype_index in range(len(AG_ratio_WT_index[ratiokey])):
                WT_select_sample_index[ratiokey].append(0)
                if len(AG_ratio_WT_index[ratiokey][datatype_index])>max_WT_sample_perdatatype[datatype_index]:
                     max_WT_sample_perdatatype[datatype_index] = len(AG_ratio_WT_index[ratiokey][datatype_index])
        
        KO_select_sample_index={}
        max_KO_sample_perdatatype=np.zeros(KO_candidate_file_num)
        for ratiokey in AG_ratio_KO_index:
           KO_select_sample_index[ratiokey]=[]
           for datatype_index in range(len(AG_ratio_KO_index[ratiokey])):
                KO_select_sample_index[ratiokey].append(0)
                if len(AG_ratio_KO_index[ratiokey][datatype_index])>max_KO_sample_perdatatype[datatype_index]:
                     max_KO_sample_perdatatype[datatype_index] = len(AG_ratio_KO_index[ratiokey][datatype_index])
        
        ####use all training samples in candidate WT and KO end
        while 1:
           for datatype in range(real_candidate_file_num):
              for time_ in range(real_candidate_file_size[datatype]):
               print("time:"+str(time_)+"  datatypeindex="+str(datatype))
               sample_weights=[]
               #IVT data start
               x_TN_batch = IVT_neg_h5fs['X'][index_IVT_TN:(index_IVT_TN + batch_size),win_center_index-nt:win_center_index+nt+1,:featuredim]
               y_ref_TN_batch = IVT_neg_h5fs['y_ref'][index_IVT_TN:(index_IVT_TN + batch_size)]
               y_call_TN_batch = IVT_neg_h5fs['y_call'][index_IVT_TN:(index_IVT_TN + batch_size)]
               sample_weights.extend(np.ones(len(x_TN_batch)))
               x_TP_batch = IVT_pos_h5fs['X'][index_IVT_TP:(index_IVT_TP + batch_size),win_center_index-nt:win_center_index+nt+1,:featuredim]
               y_ref_TP_batch = IVT_pos_h5fs['y_ref'][index_IVT_TP:(index_IVT_TP+ batch_size)]
               y_call_TP_batch = IVT_pos_h5fs['y_call'][index_IVT_TP:(index_IVT_TP + batch_size)]
               sample_weights.extend(np.ones(len(x_TP_batch)))
               index_IVT_TP += batch_size
               index_IVT_TP %= samplesize_TP
               index_IVT_TN += batch_size
               index_IVT_TN %= samplesize_TN
               #IVT data end
               #KO noncandidate start
               x_KO_batch,y_ref_KO_batch,y_call_KO_batch=[],[],[]
               for _ in range(len(KO_noncandidate_h5fs)):
                   #print("index_noncandidate_KO[_]="+str(index_noncandidate_KO[_])+"\n")
                   x_KO_batch.append(KO_noncandidate_h5fs[_]['X'][index_noncandidate_KO[_]:(index_noncandidate_KO[_] + \
                   batch_size_KO_noncandidate),win_center_index-nt:win_center_index+nt+1,:featuredim])
                   y_ref_KO_batch.append(KO_noncandidate_h5fs[_]['y_ref'][index_noncandidate_KO[_]:(index_noncandidate_KO[_] + batch_size_KO_noncandidate)])
                   y_call_KO_batch.append(KO_noncandidate_h5fs[_]['y_call'][index_noncandidate_KO[_]:(index_noncandidate_KO[_] + batch_size_KO_noncandidate)])
                   
                   index_noncandidate_KO[_] +=batch_size_KO_noncandidate
                   if index_noncandidate_KO[_]>KO_noncandidate_samplesize[_]:
                       print("debug KO noncandidate to the end size. in file \n")
                   
                   index_noncandidate_KO[_] %= KO_noncandidate_samplesize[_]
               
               x_KO_batch = np.concatenate(x_KO_batch,0)
               #print("x_KO_batch="+str(len(x_KO_batch)))
               #print("current sample_weights 1="+str(len(sample_weights))+"\n")
               sample_weights.extend(np.ones(len(x_KO_batch)))
               y_ref_KO_batch = np.concatenate(y_ref_KO_batch,0)
               y_call_KO_batch = np.concatenate(y_call_KO_batch,0)
               #KO noncandidate end
               #real noncandidate start
               X_real_noncandidate,y_ref_real_noncandidate,y_call_real_noncandidate=[],[],[]
               for _ in range(len(real_noncandidate_h5fs[datatype])):
                   X_real_noncandidate.append(real_noncandidate_h5fs[datatype][_]['X'][index_noncandidate_real[datatype][_]:(index_noncandidate_real[datatype][_] + \
                   batch_size_real_noncandidate[datatype]),win_center_index-nt:win_center_index+nt+1,:featuredim])
                   y_ref_real_noncandidate.append(real_noncandidate_h5fs[datatype][_]['y_ref'][index_noncandidate_real[datatype][_]:(index_noncandidate_real[datatype][_] + batch_size_real_noncandidate[datatype])])
                   y_call_real_noncandidate.append(real_noncandidate_h5fs[datatype][_]['y_call'][index_noncandidate_real[datatype][_]:(index_noncandidate_real[datatype][_] + batch_size_real_noncandidate[datatype])])
                   
                   index_noncandidate_real[datatype][_] +=batch_size_real_noncandidate[datatype]
                   if index_noncandidate_real[datatype][_]>real_noncandidate_samplesize[datatype][_]:
                       print("debug real noncandidate to the end size. in file "+real_noncandidate_files.split(";")[datatype].split(",")[_]+"\n")
                   
                   index_noncandidate_real[datatype][_] %= real_noncandidate_samplesize[datatype][_]
               
               X_real_noncandidate = np.concatenate(X_real_noncandidate,0)
               #print("X_real_noncandidate="+str(len(X_real_noncandidate)))
               sample_weights.extend(np.ones(len(X_real_noncandidate)))
               #print("current sample_weights2="+str(len(sample_weights))+"\n")
               y_ref_real_noncandidate = np.concatenate(y_ref_real_noncandidate,0)
               y_call_real_noncandidate = np.concatenate(y_call_real_noncandidate,0)
               #KO noncandidate end
               #real ratio data
               if binselecttype==1:
                    print("binselecttype==1 select from uniform\n")
                    ratio_select = np.random.choice(ratio_start,1,replace=False)[0] #select one ratio from ratio list
               elif binselecttype==2: #use start ratio to select bin
                    print("binselecttype==2 select from 1-ratio_start\n")
                    ratio_select = np.random.choice(ratio_start,1,replace=False,p=softmax(1-np.asarray(ratio_start)))[0] #select one ratio from ratio list
               elif binselecttype==3: #use samplenumbers to select bin, bin with more samples, has high posibility to be selected
                    print("binselecttype==3 select from samplenum\n")
                    ratio_select = np.random.choice(ratio_start,1,replace=False,p=softmax([np.log10(x) for x in ratio_bin_samplenum[datatype]]))[0] #select one ratio from ratio list
               elif binselecttype==4:#uniform but select in an order
                    print("binselecttype==4 select from uniform but in increased order\n")
                    ratio_select = ratio_start[ratio_start_index__]
                    ratio_start_index__+=1
                    ratio_start_index__=ratio_start_index__%binsample
               
               print("select ratio"+str(ratio_select)+"\n")
               x_real,y_ref_real,y_call_real,mod_label_real,ratio_label_real=[],[],[],[],[]
               if ratio_select in AG_ratio_WT_index:
                  if len(AG_ratio_WT_index[ratio_select][datatype])==0:
                           print("skip file\n")
                           continue
                  
                  select_sample_size=np.min([batch_size,len(AG_ratio_WT_index[ratio_select][datatype])])
                  select_sampleindex_list = np.random.choice(AG_ratio_WT_index[ratio_select][datatype],select_sample_size,replace=False)#slow
                  for index in select_sampleindex_list:
                       info = real_candidate_h5fs[datatype]['info'][index]
                       #print(info)
                       label = info.split("\t")[0]
                       label1 = label.split(":")[0]
                       ratio = float(label.split(":")[1])
                       sample_weight = 1-np.abs((ratio-1)*np.log(1-ratio+1e-07)-ratio*np.log(ratio+1e-07))**weight_exp
                       sample_weights.append(sample_weight)
                       if WT_label == 1:
                            mod_label_real.append(1)
                       else:
                            mod_label_real.append(ratio)
                       
                       ratio_label_real.append(ratio)
                       x_real.append(real_candidate_h5fs[datatype]['X'][index,win_center_index-nt:win_center_index+nt+1,:featuredim])
                       y_ref_real.append(real_candidate_h5fs[datatype]['y_ref'][index]) 
                       y_call_real.append(real_candidate_h5fs[datatype]['y_call'][index])
               
               #for _ in range(KO_candidate_file_num): same as real_candidate_file_num
               if ratio_select in AG_ratio_KO_index:
                  if len(AG_ratio_KO_index[ratio_select][datatype])==0:
                           print("skip file\n")
                           continue
                  
                  select_sample_size=np.min([batch_size,len(AG_ratio_KO_index[ratio_select][datatype])])
                  select_sampleindex_list = np.random.choice(AG_ratio_KO_index[ratio_select][datatype],select_sample_size,replace=False) #not very sloww!
                  for index in select_sampleindex_list:
                      info = KO_candidate_h5fs[datatype]['info'][index]
                      label = info.split("\t")[0]
                      ratio = float(label.split(":")[1])
                      sample_weights.append(1)
                      mod_label_real.append(0)
                      ratio_label_real.append(0)
                      x_real.append(KO_candidate_h5fs[datatype]['X'][index,win_center_index-nt:win_center_index+nt+1,:featuredim])
                      y_ref_real.append(KO_candidate_h5fs[datatype]['y_ref'][index]) 
                      y_call_real.append(KO_candidate_h5fs[datatype]['y_call'][index])
               
               sample_weights = np.asarray(sample_weights)
               #print("current sample_weights3="+str(len(sample_weights))+"\n")
               x = np.concatenate([x_TN_batch,x_TP_batch,x_KO_batch,X_real_noncandidate,x_real],axis = 0)
               y_ref = np.concatenate([y_ref_TN_batch,y_ref_TP_batch,y_ref_KO_batch,y_ref_real_noncandidate,y_ref_real],axis = 0)
               y_call = np.concatenate([y_call_TN_batch,y_call_TP_batch,y_call_KO_batch,y_call_real_noncandidate,y_call_real],axis = 0)
               mod_label = np.concatenate([np.zeros(len(x_TN_batch)),np.ones(len(x_TP_batch)),np.zeros(len(x_KO_batch)),np.zeros(len(X_real_noncandidate)),mod_label_real],axis=0)
               ratio_label = np.concatenate([np.zeros(len(x_TN_batch)),np.ones(len(x_TP_batch))+0.1,np.zeros(len(x_KO_batch)),np.zeros(len(X_real_noncandidate)),ratio_label_real],axis=0)
               sample_weights_base = np.ones(len(sample_weights))
               print("pseudo_label\n")
               if sample_weights_flag:
                       yield([x,np.expand_dims(y_ref,axis=-1),np.expand_dims(y_call,axis=-1)],[np.expand_dims(y_ref,axis=-1), np.expand_dims(y_call, -1),np.expand_dims(ratio_label,-1),np.expand_dims(ratio_label,-1)],[sample_weights_base,sample_weights_base,sample_weights,sample_weights])
               else:
                       yield([x,np.expand_dims(y_ref,axis=-1),np.expand_dims(y_call,axis=-1)],[np.expand_dims(y_ref,axis=-1), np.expand_dims(y_call, -1),np.expand_dims(ratio_label,-1),np.expand_dims(ratio_label,-1)])

def generator_realdata(IVT_pos_file,IVT_neg_file,real_candidate_files,real_noncandidate_files,KO_candidate_files,KO_noncandidate_files,batch_size,bin=0.2,binsample=200,WT_label=1,
        binselecttype=1,sample_weights_flag=True,featuredim=3,windowsize=9,weight_exp=1,KO_noncandidate_weight=1,real_noncandidate_weight=1): #fixed start np.arnage(0,1,bin)
        print("generator_all")
        real_candidate_files=real_candidate_files.replace("\"",'')
        real_noncandidate_files=real_noncandidate_files.replace("\"",'')
        KO_candidate_files=KO_candidate_files.replace("\"",'')
        KO_noncandidate_files=KO_noncandidate_files.replace("\"",'')
        
        #KO noncandidate files start
        KO_noncandidate_h5fs,KO_noncandidate_samplesize,index_noncandidate_KO=[],[],[]
        index=0
        for KO_noncandidate in KO_noncandidate_files.split(","):
            if os.path.isfile(KO_noncandidate):
                temp = h5py.File(KO_noncandidate,"r")
                if len(temp['X'])!=0:
                       KO_noncandidate_h5fs.append(temp)
                       temp_size=len(KO_noncandidate_h5fs[index]['X'])
                       KO_noncandidate_samplesize.append(temp_size)
                       select_index = np.random.choice(np.arange(0,temp_size,int(temp_size/10000)))
                       index_noncandidate_KO.append(select_index)
                       print("debug "+KO_noncandidate+" index_noncandidate_KO="+str(select_index)+"\n")
                       index+=1
                else:
                       print("debug error: KO_noncandidate_samplesize="+str(temp['X']))
                       print(KO_noncandidate)
        
        batch_size_KO_noncandidate=int(batch_size*KO_noncandidate_weight/len(KO_noncandidate_h5fs))
        print("debug batch_size_KO_noncandidate="+str(batch_size_KO_noncandidate))
        print("KO_noncandidate_file_num="+str(len(KO_noncandidate_h5fs)))
        #KO noncandidate files end
        #real noncandidate files start
        real_noncandidate_h5fs,real_noncandidate_samplesize,index_noncandidate_real={},{},{}
        batch_size_real_noncandidate={}
        for datatype in range(len(real_noncandidate_files.split(";"))): #use ; for different datatype
           real_noncandidate_h5fs[datatype],real_noncandidate_samplesize[datatype],index_noncandidate_real[datatype]=[],[],[]
           index=0
           for real_noncandidate in real_noncandidate_files.split(";")[datatype].split(","):
              if os.path.isfile(real_noncandidate):
                temp = h5py.File(real_noncandidate,"r")
                if len(temp['X'])==0:
                       print("Error: real_noncandidate_h5fs="+str(len(temp['X'])))
                       print(real_noncandidate)
                else:
                     real_noncandidate_h5fs[datatype].append(temp)
                     temp_size=len(real_noncandidate_h5fs[datatype][index]['X'])
                     real_noncandidate_samplesize[datatype].append(temp_size)
                     select_index = np.random.choice(np.arange(0,temp_size,int(temp_size/10000)))
                     index_noncandidate_real[datatype].append(select_index)
                     print("debug "+real_noncandidate+" index_noncandidate_real="+str(select_index)+"\n")
                     index+=1
           
           batch_size_real_noncandidate[datatype]=int(batch_size*real_noncandidate_weight/len(real_noncandidate_h5fs[datatype]))
           print("real_noncandidate_file_num="+str(len(real_noncandidate_h5fs[datatype])))
        
        #real noncandidate files end
        win_center_index=int(len(KO_noncandidate_h5fs[0]['X'][0])/2)
        nt = int(windowsize/2)
        
        AG_ratio_WT,AG_ratio_KO = [],[]
        real_candidate_h5fs=[]
        KO_candidate_h5fs=[]
        index=0
        real_candidate_file_size=[]
        for real_candidate in real_candidate_files.split(","):
          if os.path.isfile(real_candidate):
            print("openning "+real_candidate+"\n")
            temp = h5py.File(real_candidate,'r')
            file_length = len(temp['X'])
            real_candidate_file_size.append(file_length)
            if file_length==0:
                 print("Error: real candidate sample size is 0 ")
                 print(real_candidate)
                 real_candidate_file_size[index]=1 #avoid devide by 0
            else:
                AG_ratio_WT.append([])
                real_candidate_h5fs.append(temp)
                for line in real_candidate_h5fs[index]['info']:
                    try:
                         ratio=float(line.split("\t")[0].split(":")[1])
                    except:
                         print("Error:")
                         print(line)
                         exit()
                    
                    AG_ratio_WT[index].append(ratio)
                
                index+=1
        
        minsize=np.min(real_candidate_file_size)
        real_candidate_file_size=[int(np.ceil(x/minsize)) for x in real_candidate_file_size]
        print("debug real_candidate_file_size\n")
        print(real_candidate_file_size)
        real_candidate_file_num = index
        print("real_candidate_file_num="+str(real_candidate_file_num))
        index=0
        for KO_candidate in KO_candidate_files.split(","):
          if os.path.isfile(KO_candidate):
            print("openning "+KO_candidate+"\n")
            temp=h5py.File(KO_candidate,'r')
            if len(temp['X'])==0:
                 print("Error: KO candidate sample size is 0")
                 print(KO_candidate)
            else:
                KO_candidate_h5fs.append(temp)
                AG_ratio_KO.append([])
                for line in KO_candidate_h5fs[index]['info']:
                    AG_ratio_KO[index].append(float(line.split("\t")[0].split(":")[1]))
                
                index+=1
        
        KO_candidate_file_num = index
        print("KO_candidate_file_num="+str(KO_candidate_file_num))
        ratio_start = [round(x,4) for x in np.random.uniform(0.0, 1.0, binsample)]
        if binselecttype==4:
           ratio_start_index__=0
        
        ratio_bin_samplenum=np.zeros([real_candidate_file_num,len(ratio_start)])
        AG_ratio_WT_index,AG_ratio_KO_index = {},{}
        ratio_start_index=0
        for start in ratio_start:
           end = start+bin
           AG_ratio_WT_index[start]=[]
           flag_AG_ratio =0
           for _ in range(real_candidate_file_num):
                 index_start_end_WT = np.where((AG_ratio_WT[_]<=end) &(AG_ratio_WT[_]>start))[0]
                 AG_ratio_WT_index[start].append(index_start_end_WT) #AG_ratio_WT_index[start] element could be []
                 ratio_bin_samplenum[_][ratio_start_index]=len(index_start_end_WT)
                 if len(index_start_end_WT) > 50:
                       flag_AG_ratio+=len(index_start_end_WT)
           
           if flag_AG_ratio==0:
                  del AG_ratio_WT_index[start]
           
           ratio_start_index+=1
        
        for start in ratio_start:
          end = start+bin
          AG_ratio_KO_index[start]=[]
          flag_AG_ratio = 0
          for _ in range(KO_candidate_file_num):
            index_start_end_KO = np.where((AG_ratio_KO[_]<=end) &(AG_ratio_KO[_]>start))[0]
            AG_ratio_KO_index[start].append(index_start_end_KO)
            #print("len len(AG_ratio_KO_index) at AG ratio "+str(start)+"="+str(len(AG_ratio_KO_index[start][_]))+"\n")
            if len(index_start_end_KO) > 50:
                  flag_AG_ratio+=len(index_start_end_KO)
          
          if flag_AG_ratio==0: #all files ==0 
                del AG_ratio_KO_index[start]
        
        for index,start in enumerate(ratio_start):
             if (start not in AG_ratio_KO_index) and (start not in AG_ratio_WT_index):
                    del ratio_start[index]
        
        #use all training samples in candidate WT and KO
        WT_select_sample_index={}
        max_WT_sample_perdatatype=np.zeros(real_candidate_file_num)
        for ratiokey in AG_ratio_WT_index:
           WT_select_sample_index[ratiokey]=[]
           for datatype_index in range(len(AG_ratio_WT_index[ratiokey])):
                WT_select_sample_index[ratiokey].append(0)
                if len(AG_ratio_WT_index[ratiokey][datatype_index])>max_WT_sample_perdatatype[datatype_index]:
                     max_WT_sample_perdatatype[datatype_index] = len(AG_ratio_WT_index[ratiokey][datatype_index])
        
        KO_select_sample_index={}
        max_KO_sample_perdatatype=np.zeros(KO_candidate_file_num)
        for ratiokey in AG_ratio_KO_index:
           KO_select_sample_index[ratiokey]=[]
           for datatype_index in range(len(AG_ratio_KO_index[ratiokey])):
                KO_select_sample_index[ratiokey].append(0)
                if len(AG_ratio_KO_index[ratiokey][datatype_index])>max_KO_sample_perdatatype[datatype_index]:
                     max_KO_sample_perdatatype[datatype_index] = len(AG_ratio_KO_index[ratiokey][datatype_index])
        
        ####use all training samples in candidate WT and KO end
        while 1:
           for datatype in range(real_candidate_file_num):
              for time_ in range(real_candidate_file_size[datatype]):
               print("time:"+str(time_)+"  datatypeindex="+str(datatype))
               sample_weights=[]
               #KO noncandidate start
               x_KO_batch,y_ref_KO_batch,y_call_KO_batch=[],[],[]
               for _ in range(len(KO_noncandidate_h5fs)):
                   #print("index_noncandidate_KO[_]="+str(index_noncandidate_KO[_])+"\n")
                   x_KO_batch.append(KO_noncandidate_h5fs[_]['X'][index_noncandidate_KO[_]:(index_noncandidate_KO[_] + \
                   batch_size_KO_noncandidate),win_center_index-nt:win_center_index+nt+1,:featuredim])
                   y_ref_KO_batch.append(KO_noncandidate_h5fs[_]['y_ref'][index_noncandidate_KO[_]:(index_noncandidate_KO[_] + batch_size_KO_noncandidate)])
                   y_call_KO_batch.append(KO_noncandidate_h5fs[_]['y_call'][index_noncandidate_KO[_]:(index_noncandidate_KO[_] + batch_size_KO_noncandidate)])
                   
                   index_noncandidate_KO[_] +=batch_size_KO_noncandidate
                   if index_noncandidate_KO[_]>KO_noncandidate_samplesize[_]:
                       print("debug KO noncandidate to the end size. in file \n")
                   
                   index_noncandidate_KO[_] %= KO_noncandidate_samplesize[_]
               
               x_KO_batch = np.concatenate(x_KO_batch,0)
               sample_weights.extend(np.ones(len(x_KO_batch)))
               y_ref_KO_batch = np.concatenate(y_ref_KO_batch,0)
               y_call_KO_batch = np.concatenate(y_call_KO_batch,0)
               #KO noncandidate end
               #real noncandidate start
               X_real_noncandidate,y_ref_real_noncandidate,y_call_real_noncandidate=[],[],[]
               for _ in range(len(real_noncandidate_h5fs[datatype])):
                   X_real_noncandidate.append(real_noncandidate_h5fs[datatype][_]['X'][index_noncandidate_real[datatype][_]:(index_noncandidate_real[datatype][_] + \
                   batch_size_real_noncandidate[datatype]),win_center_index-nt:win_center_index+nt+1,:featuredim])
                   y_ref_real_noncandidate.append(real_noncandidate_h5fs[datatype][_]['y_ref'][index_noncandidate_real[datatype][_]:(index_noncandidate_real[datatype][_] + batch_size_real_noncandidate[datatype])])
                   y_call_real_noncandidate.append(real_noncandidate_h5fs[datatype][_]['y_call'][index_noncandidate_real[datatype][_]:(index_noncandidate_real[datatype][_] + batch_size_real_noncandidate[datatype])])
                   
                   index_noncandidate_real[datatype][_] +=batch_size_real_noncandidate[datatype]
                   if index_noncandidate_real[datatype][_]>real_noncandidate_samplesize[datatype][_]:
                       print("debug real noncandidate to the end size. in file "+real_noncandidate_files.split(";")[datatype].split(",")[_]+"\n")
                   
                   index_noncandidate_real[datatype][_] %= real_noncandidate_samplesize[datatype][_]
               
               X_real_noncandidate = np.concatenate(X_real_noncandidate,0)
               #print("X_real_noncandidate="+str(len(X_real_noncandidate)))
               sample_weights.extend(np.ones(len(X_real_noncandidate)))
               #print("current sample_weights2="+str(len(sample_weights))+"\n")
               y_ref_real_noncandidate = np.concatenate(y_ref_real_noncandidate,0)
               y_call_real_noncandidate = np.concatenate(y_call_real_noncandidate,0)
               #KO noncandidate end
               #real ratio data
               if binselecttype==1:
                    print("binselecttype==1 select from uniform\n")
                    ratio_select = np.random.choice(ratio_start,1,replace=False)[0] #select one ratio from ratio list
               elif binselecttype==2: #use start ratio to select bin
                    print("binselecttype==2 select from 1-ratio_start\n")
                    ratio_select = np.random.choice(ratio_start,1,replace=False,p=softmax(1-np.asarray(ratio_start)))[0] #select one ratio from ratio list
               elif binselecttype==3: #use samplenumbers to select bin, bin with more samples, has high posibility to be selected
                    print("binselecttype==3 select from samplenum\n")
                    ratio_select = np.random.choice(ratio_start,1,replace=False,p=softmax([np.log10(x) for x in ratio_bin_samplenum[datatype]]))[0] #select one ratio from ratio list
               elif binselecttype==4:
                    print("binselecttype==4 select from uniform but in increased order\n")
                    ratio_select = ratio_start[ratio_start_index__]
                    ratio_start_index__+=1
                    ratio_start_index__=ratio_start_index__%binsample
               
               print("select ratio"+str(ratio_select)+"\n")
               x_real,y_ref_real,y_call_real,mod_label_real,ratio_label_real=[],[],[],[],[]
               if ratio_select in AG_ratio_WT_index:
                  if len(AG_ratio_WT_index[ratio_select][datatype])==0:
                           print("skip file\n")
                           continue
                  
                  select_sample_size=np.min([batch_size,len(AG_ratio_WT_index[ratio_select][datatype])])
                  select_sampleindex_list = np.random.choice(AG_ratio_WT_index[ratio_select][datatype],select_sample_size,replace=False)#slow
                  for index in select_sampleindex_list:
                       info = real_candidate_h5fs[datatype]['info'][index]
                       #print(info)
                       label = info.split("\t")[0]
                       label1 = label.split(":")[0]
                       ratio = float(label.split(":")[1])
                       sample_weight = 1-np.abs((ratio-1)*np.log(1-ratio+1e-07)-ratio*np.log(ratio+1e-07))**weight_exp
                       sample_weights.append(sample_weight)
                       if WT_label == 1:
                            mod_label_real.append(1)
                       else:
                            mod_label_real.append(ratio)
                       
                       ratio_label_real.append(ratio)
                       x_real.append(real_candidate_h5fs[datatype]['X'][index,win_center_index-nt:win_center_index+nt+1,:featuredim])
                       y_ref_real.append(real_candidate_h5fs[datatype]['y_ref'][index]) 
                       y_call_real.append(real_candidate_h5fs[datatype]['y_call'][index])
               
               #for _ in range(KO_candidate_file_num): same as real_candidate_file_num
               if ratio_select in AG_ratio_KO_index:
                  if len(AG_ratio_KO_index[ratio_select][datatype])==0:
                           print("skip file\n")
                           continue
                  
                  select_sample_size=np.min([batch_size,len(AG_ratio_KO_index[ratio_select][datatype])])
                  select_sampleindex_list = np.random.choice(AG_ratio_KO_index[ratio_select][datatype],select_sample_size,replace=False) #not very sloww!
                  for index in select_sampleindex_list:
                      info = KO_candidate_h5fs[datatype]['info'][index]
                      label = info.split("\t")[0]
                      ratio = float(label.split(":")[1])
                      sample_weights.append(1)
                      mod_label_real.append(0)
                      ratio_label_real.append(0)
                      x_real.append(KO_candidate_h5fs[datatype]['X'][index,win_center_index-nt:win_center_index+nt+1,:featuredim])
                      y_ref_real.append(KO_candidate_h5fs[datatype]['y_ref'][index]) 
                      y_call_real.append(KO_candidate_h5fs[datatype]['y_call'][index])
               
               sample_weights = np.asarray(sample_weights)
               x = np.concatenate([x_KO_batch,X_real_noncandidate,x_real],axis = 0)
               y_ref = np.concatenate([y_ref_KO_batch,y_ref_real_noncandidate,y_ref_real],axis = 0)
               y_call = np.concatenate([y_call_KO_batch,y_call_real_noncandidate,y_call_real],axis = 0)
               mod_label = np.concatenate([np.zeros(len(x_KO_batch)),np.zeros(len(X_real_noncandidate)),mod_label_real],axis=0)
               ratio_label = np.concatenate([np.zeros(len(x_KO_batch)),np.zeros(len(X_real_noncandidate)),ratio_label_real],axis=0)
               sample_weights_base = np.ones(len(sample_weights))
               print("pseudo_label\n")
               if sample_weights_flag:
                       yield([x,np.expand_dims(y_ref,axis=-1),np.expand_dims(y_call,axis=-1)],[np.expand_dims(y_ref,axis=-1), np.expand_dims(y_call, -1),np.expand_dims(ratio_label,-1),np.expand_dims(ratio_label,-1)],[sample_weights_base,sample_weights_base,sample_weights,sample_weights])
               else:
                       yield([x,np.expand_dims(y_ref,axis=-1),np.expand_dims(y_call,axis=-1)],[np.expand_dims(y_ref,axis=-1), np.expand_dims(y_call, -1),np.expand_dims(ratio_label,-1),np.expand_dims(ratio_label,-1)])

def generator_IVT_shuffle(IVT_pos_file,IVT_neg_file,real_candidate_files,real_noncandidate_files,KO_candidate_files,KO_noncandidate_files,batch_size,bin=0.2,binsample=200,WT_label=1,
        binselecttype=1,sample_weights_flag=True,featuredim=3,windowsize=9): #fixed start np.arnage(0,1,bin)
        print("generator_IVT_shuffle")
        #IVT data
        IVT_pos_h5fs=h5py.File(IVT_pos_file,"r")
        IVT_neg_h5fs=h5py.File(IVT_neg_file,"r")
        win_center_index=int(len(IVT_neg_h5fs['X'][0])/2)
        nt = int(windowsize/2)
        samplesize_TN = len(IVT_neg_h5fs['X'])
        samplesize_TP = len(IVT_pos_h5fs['X'])
        TN_sample_list = np.arange(samplesize_TN)
        TP_sample_list = np.arange(samplesize_TP)
        index_IVT_TP = 0
        index_IVT_TN = 0
        print("samplesize_TN="+str(samplesize_TN)+"\n")
        print("samplesize_TP="+str(samplesize_TP)+"\n")
        while 1:
            if index_IVT_TP < batch_size:
                 print("shuffle TP dataset\n")
                 np.random.shuffle(TP_sample_list)
            
            if index_IVT_TN < batch_size:
                 print("shuffle TN dataset\n")
                 np.random.shuffle(TN_sample_list)
            
            x_TN_batch=[]
            x_TP_batch=[]
            y_ref_TN_batch=[]
            y_call_TN_batch=[]
            y_ref_TP_batch=[]
            y_call_TP_batch=[]
            
            for _ in range(batch_size):
                index=TN_sample_list[(index_IVT_TN+_)%samplesize_TN]
                #IVT data start
                x_TN_batch.append(IVT_neg_h5fs['X'][index,win_center_index-nt:win_center_index+nt+1,:featuredim])
                y_ref_TN_batch.append(IVT_neg_h5fs['y_ref'][index])
                y_call_TN_batch.append(IVT_neg_h5fs['y_call'][index])
                
                index=TP_sample_list[(index_IVT_TP+_)%samplesize_TP]
                x_TP_batch.append(IVT_pos_h5fs['X'][index,win_center_index-nt:win_center_index+nt+1,:featuredim])
                y_ref_TP_batch.append(IVT_pos_h5fs['y_ref'][index])
                y_call_TP_batch.append(IVT_pos_h5fs['y_call'][index])
                
            
            index_IVT_TP += batch_size
            index_IVT_TP %= samplesize_TP
            index_IVT_TN += batch_size
            index_IVT_TN %= samplesize_TN
            #IVT data end
            
            x = np.concatenate([x_TN_batch,x_TP_batch],axis = 0)
            y_ref = np.concatenate([y_ref_TN_batch,y_ref_TP_batch],axis = 0)
            y_call = np.concatenate([y_call_TN_batch,y_call_TP_batch],axis = 0)
            mod_label = np.concatenate([np.zeros(len(x_TN_batch)),np.ones(len(x_TP_batch))+0.1],axis=0)
            ratio_label = np.concatenate([np.zeros(len(x_TN_batch)),np.ones(len(x_TP_batch))+0.1],axis=0)
            yield([x,np.expand_dims(y_ref,axis=-1),np.expand_dims(y_call,axis=-1)],[np.expand_dims(y_ref,axis=-1), np.expand_dims(y_call, -1),mod_label,np.expand_dims(ratio_label,-1)])


def generator_IVT(IVT_pos_file,IVT_neg_file,real_candidate_files,real_noncandidate_files,KO_candidate_files,KO_noncandidate_files,batch_size,bin=0.2,binsample=200,WT_label=1,
        binselecttype=1,sample_weights_flag=True,featuredim=3,windowsize=9): #fixed start np.arnage(0,1,bin)
        print("generator_IVT")
        #IVT data
        IVT_pos_h5fs=h5py.File(IVT_pos_file,"r")
        IVT_neg_h5fs=h5py.File(IVT_neg_file,"r")
        win_center_index=int(len(IVT_neg_h5fs['X'][0])/2)
        nt = int(windowsize/2)
        samplesize_TN = len(IVT_neg_h5fs['X'])
        samplesize_TP = len(IVT_pos_h5fs['X'])
        index_IVT_TP = 0
        index_IVT_TN = 0
        while 1:
            #IVT data start
            x_TN_batch = IVT_neg_h5fs['X'][index_IVT_TN:(index_IVT_TN + batch_size),win_center_index-nt:win_center_index+nt+1,:featuredim]
            y_ref_TN_batch = IVT_neg_h5fs['y_ref'][index_IVT_TN:(index_IVT_TN + batch_size)]
            y_call_TN_batch = IVT_neg_h5fs['y_call'][index_IVT_TN:(index_IVT_TN + batch_size)]
            x_TP_batch = IVT_pos_h5fs['X'][index_IVT_TP:(index_IVT_TP + batch_size),win_center_index-nt:win_center_index+nt+1,:featuredim]
            y_ref_TP_batch = IVT_pos_h5fs['y_ref'][index_IVT_TP:(index_IVT_TP+ batch_size)]
            y_call_TP_batch = IVT_pos_h5fs['y_call'][index_IVT_TP:(index_IVT_TP + batch_size)]
            index_IVT_TP += batch_size
            index_IVT_TP %= samplesize_TP
            index_IVT_TN += batch_size
            index_IVT_TN %= samplesize_TN
            #IVT data end
            
            x = np.concatenate([x_TN_batch,x_TP_batch],axis = 0)
            y_ref = np.concatenate([y_ref_TN_batch,y_ref_TP_batch],axis = 0)
            y_call = np.concatenate([y_call_TN_batch,y_call_TP_batch],axis = 0)
            mod_label = np.concatenate([np.zeros(len(x_TN_batch)),np.ones(len(x_TP_batch))+0.1],axis=0)
            ratio_label = np.concatenate([np.zeros(len(x_TN_batch)),np.ones(len(x_TP_batch))+0.1],axis=0)
            yield([x,np.expand_dims(y_ref,axis=-1),np.expand_dims(y_call,axis=-1)],[np.expand_dims(y_ref,axis=-1), np.expand_dims(y_call, -1),mod_label,np.expand_dims(ratio_label,-1)])


def batch_predict(hdf5filelist,modellist,outputfile,maxoutputlen=13,batch_size=1000,verbose=0,featuredim=3,windowsize=9,printref=False):
    output = open(outputfile,'w')
    hdf5filelist=hdf5filelist.replace("\"",'')
    hdf5files = hdf5filelist.split(',')
    for hdf5file in hdf5files:
         try:
           h5f = h5py.File(hdf5file,"r")
         except:
             print('read '+hdf5file+" error.")
             continue
         win_center_index=int(len(h5f['X'][0])/2)
         nt = int(windowsize/2)
         total_size = len(h5f['X'])
         for ndx in range(0,total_size,batch_size):
           prossratio = float(ndx)/total_size
           if verbose==1:
              print(str(prossratio))
           
           batch_data = h5f['X'][ndx:(ndx + batch_size),win_center_index-nt:win_center_index+nt+1,:featuredim]#batch_generator.next()
           batch_info = h5f['info'][ndx:(ndx+batch_size)]
           if printref:
                batch_y_ref = h5f['y_ref'][ndx:(ndx+batch_size)]
           
           fake_input = np.zeros([len(batch_data),maxoutputlen,1])
           predict_mod = np.zeros(len(batch_data))
           for model in modellist:
                 predicty = model.model.predict([batch_data,fake_input,fake_input],batch_size=batch_size)            
                 predict_mod += predicty[2][:,1]
           
           predict_mod=predict_mod/len(modellist)
           
           if printref:
               for pindex in range(len(predict_mod)):
                  output.write(batch_info[pindex]+"\t"+str("".join([encoding_baseout[x] for x in batch_y_ref[pindex][:-2]]))+"\t"+str(predict_mod[pindex])+"\n")
           else:
              for pindex in range(len(predict_mod)):
                  if len(batch_info[pindex])==0:
                       print("debug batch_info error in file "+hdf5file+" ndx="+str(ndx)+" pindex="+str(pindex)+"\n")
                  output.write(batch_info[pindex]+"\t"+str(predict_mod[pindex])+"\n")
         
         h5f.close()

    
    output.close()

def ratio_loss(y_true,y_pred):
    y_true_rank = ops.convert_to_tensor(y_true).get_shape().ndims
    #If y_true is (num_samples,)
    if y_true_rank==1:
          y_true = tf.expand_dims(y_true,-1)
    
    mask_value_1 = K.variable(1.1)
    mask_value_0 = K.variable(0.0)
    mask = K.any(tf.concat([K.equal(y_true, mask_value_0),K.equal(y_true,mask_value_1)],-1), axis=-1) #mask samples that have 0 or 1 labels
    mask = 1 - K.cast(mask, K.floatx())
    def ratio_loss_temp():
        #count y_pred 1 num and 0 num use cutoff 0.5
        ratio_true_sum = K.sum(tf.squeeze(y_true) * mask)
        y_pred_label = K.all(K.greater(y_pred,K.variable(0.5)),axis = -1)
        y_pred_label = K.cast(y_pred_label, K.floatx())
        y_pred_WT_sum =  K.sum(tf.squeeze(y_pred_label) * mask)
        return K.abs(y_pred_WT_sum-ratio_true_sum)/K.sum(mask)
    def no_ratio_loss():
        return 0.0
    
    return tf.cond(K.equal(K.sum(mask),0), no_ratio_loss,ratio_loss_temp)


def neg_ratio_loss(y_true,y_pred):
    y_true_rank = ops.convert_to_tensor(y_true).get_shape().ndims
    #If y_true is (num_samples,)
    if y_true_rank==1:
          y_true = tf.expand_dims(y_true,-1)
    
    mask_value_0 = K.variable(0.0)
    mask_keep0=K.cast(K.all(K.equal(y_true, mask_value_0),axis=-1),K.floatx())
    def neg_ratio_loss_temp():
        y_pred_label = K.all(K.greater(y_pred,K.variable(0.5)),axis = -1)
        y_pred_label = K.cast(y_pred_label, K.floatx())
        y_pred_neg_sum=K.sum(tf.squeeze(y_pred_label) * mask_keep0)
        return K.abs(y_pred_neg_sum)/K.sum(mask_keep0)
    def no_ratio_loss():
        return 0.0
    
    return tf.cond(K.equal(K.sum(mask_keep0),0), no_ratio_loss,neg_ratio_loss_temp)

def pos_ratio_loss(y_true,y_pred):
    y_true_rank = ops.convert_to_tensor(y_true).get_shape().ndims
    #If y_true is (num_samples,)
    if y_true_rank==1:
          y_true = tf.expand_dims(y_true,-1)
    
    mask_value_1 = K.variable(1.1)
    mask_keep1=K.cast(K.all(K.equal(y_true, mask_value_1),axis=-1),K.floatx())
    def pos_ratio_loss_temp():
        y_pred_label = K.all(K.greater(y_pred,K.variable(0.5)),axis = -1)
        y_pred_label = K.cast(y_pred_label, K.floatx())
        y_pred_pos_sum=K.sum(tf.squeeze(y_pred_label) * mask_keep1)
        return K.abs(y_pred_pos_sum-K.sum(mask_keep1)*mask_value_1)/K.sum(mask_keep1)
    def no_ratio_loss():
        return 0.0
    
    return tf.cond(K.equal(K.sum(mask_keep1),0), no_ratio_loss,pos_ratio_loss_temp)

def self_rank_loss(y_true,y_pred): #y_true is ratio
    y_true_rank = ops.convert_to_tensor(y_true).get_shape().ndims
    #If y_true is (num_samples,)
    if y_true_rank==1:
         y_true = tf.expand_dims(y_true,-1)
    
    mask_value_1 = K.variable(1.1)
    mask_value_0 = K.variable(0.0)
    #mask = K.all(K.equal(y_true, mask_value_0),axis=-1) #y_true (100 Bitwise reduction (logical AND
    mask = K.any(tf.concat([K.equal(y_true, mask_value_0),K.equal(y_true,mask_value_1)],-1), axis=-1) #mask samples that have 0 or 1 labels
    mask = 1 - K.cast(mask, K.floatx())
    def false_fn():
        ratio_true_sum = K.sum(tf.squeeze(y_true) * mask)
        WT_numbers=K.cast(tf.floor(ratio_true_sum),dtype='int32')
        y_pred_WTmasked =  tf.squeeze(y_pred[:,1]) * mask #make all KO labels to 0 
        def allzero():
            return tf.squeeze(y_true) * mask * 0 #K.zeros(K.shape(y_pred_WTmasked))
        
        def pseudolabel():
           y_pseudo_label = K.greater_equal(y_pred_WTmasked,tf.math.top_k(y_pred_WTmasked,WT_numbers).values[-1])
           y_pseudo_label = K.cast(y_pseudo_label,K.floatx())+tf.squeeze(y_true) * (1-mask)
           return y_pseudo_label
        
        y_pseudo_label = tf.cond(K.equal(WT_numbers,0),allzero,pseudolabel)
        
        return y_pseudo_label
    
    def true_fn():
        return y_true
    
    y_pseudo_label=tf.cond(K.equal(K.sum(mask),0), true_fn,false_fn)
    y_pseudo_label = K.clip(y_pseudo_label,0,1) #for IVT data y_true = 1.1 clip it to 1
    loss=K.mean(losses.sparse_categorical_crossentropy(y_pseudo_label,y_pred))
    return loss

class multihead_attention:
    
    def __init__(self, max_len, maxoutputlen, featuredim,nb_classes, save_path, kfold_index):
        self.max_len = max_len
        self.maxoutputlen=maxoutputlen
        self.featuredim = featuredim
        self.nb_classes = nb_classes
        global OUTPATH
        OUTPATH = save_path
        self.kfold_index = kfold_index
        self.att_weight_var = K.variable(1.0)
    
    def build_model_seq2seq(self,
                            dim_attention,
                            load_weights=False, weight_dir=None,
                            drop_input=0,
                            drop_cnn=0.2,
                            dim_lstm=50,
                            embedding_output_dim=30,
                            BatchNorm=False,
                            is_monotonic=False,
                            teacherforcing=True,
                            teacher_forcing_ratio=0.5,
                            Att_regularizer_weight=0.0005,
                            headnum=10,
                            W_regularizer=0.0005,
                            drop_flat=0,
                            fc_dim=100,
                            fcnum=0,
                            loss_weights=[1., 1.,1.,1.],
                            ):
        
        center_index = int(self.max_len/2) #max_len is windowsize
        alphabet_size=4+2 #ACGTI+ end
        #print("max_len:"+str(self.max_len))
        #print("featuredim:"+str(self.featuredim))
        input = Input(shape=(self.max_len,self.featuredim), dtype='float32')
        outp_true_ref = Input(shape=(None,1), dtype='int64')
        outp_true_call = Input(shape=(None,1), dtype='int64')
        input_drop = Dropout(drop_input)(input)
        LSTM1 = Dropout(drop_cnn)(BatchNormalization()(Bidirectional(LSTM(dim_lstm, kernel_initializer="orthogonal", recurrent_initializer="orthogonal",return_sequences=True), merge_mode='sum')(input_drop)))
        LSTM2 = Dropout(drop_cnn,name="embeddinglayer")(BatchNormalization()(Bidirectional(LSTM(dim_lstm, kernel_initializer="orthogonal", recurrent_initializer="orthogonal",return_sequences=True), merge_mode='sum')(LSTM1)))
        
        if is_monotonic:
           preds1 = AttentionDecoder(name="ref_base",units=dim_attention,alphabet_size=alphabet_size,embedding_dim=embedding_output_dim,is_monotonic=True,normalize_energy=True)(LSTM2)
        else:
           preds1 = AttentionDecoder(name="ref_base",units=dim_attention,alphabet_size=alphabet_size,embedding_dim=embedding_output_dim)([LSTM2,outp_true_ref],use_teacher_forcing=teacherforcing,teacher_forcing_ratio=teacher_forcing_ratio)
        
        if is_monotonic:
           preds2 = AttentionDecoder(name="called_base",units=dim_attention,alphabet_size=alphabet_size,embedding_dim=embedding_output_dim,is_monotonic=True,normalize_energy=True)(LSTM2)
        else:
           preds2 = AttentionDecoder(name="called_base",units=dim_attention,alphabet_size=alphabet_size,embedding_dim=embedding_output_dim)([LSTM2,outp_true_call],use_teacher_forcing=teacherforcing,teacher_forcing_ratio=teacher_forcing_ratio)
        
        att = Attention(hidden=LSTM2.get_shape()[-1].value, da=dim_attention, r=headnum, init='glorot_uniform', activation='tanh',
                    W1_regularizer=regularizers.l2(W_regularizer),
                    W2_regularizer=regularizers.l2(W_regularizer),
                    W1_constraint=None, W2_constraint=None, return_attention=False,
                    attention_regularizer_weight=Att_regularizer_weight)(LSTM2)
        
        attbn = BatchNormalization()(att)
        output = Dropout(drop_flat)(Flatten()(attbn))
        fc = output
        for _ in range(fcnum):
             fc = Dense(fc_dim,activation='relu')(fc)
             fc = Dropout(drop_flat)(fc)
        
        preds3 = Dense(self.nb_classes,activation='softmax',name='mod_pred')(fc)
        
        Gslice = Lambda(lambda x, i, j: x[:,i:j], output_shape=[1], name="ratio_pred")  # Define your lambda layer
        Gslice.arguments = {'i': 1, 'j':2}  # Update extra arguments to Gslice
        # Call F
        preds3_ratio = Gslice(preds3)
        self.model = Model(inputs=[input,outp_true_ref,outp_true_call], outputs=[preds1,preds2,preds3,preds3_ratio])
        
        from keras import optimizers
        optim = optimizers.nadam()
        self.model.compile(loss=['sparse_categorical_crossentropy','sparse_categorical_crossentropy',
        self_rank_loss,ratio_loss], loss_weights=loss_weights,optimizer='adam',
        metrics={'ref_base':new_sparse_categorical_accuracy,'called_base':new_sparse_categorical_accuracy,'mod_pred':[pos_accuracy,neg_accuracy,pseudo_accuracy2,true_accuracy2]})
        if load_weights and weight_dir is not None and os.path.exists(weight_dir):
              print("loading weights: "+str(weight_dir))
              self.model.load_weights(weight_dir)
        elif load_weights:
             print("Error in loading weights, check weights "+str(weight_dir))
        
        self.bn = False
        #self.model.summary()
    
    def train(self,samplesize,real_candidate_files,real_noncandidate_files,KO_candidate_files,KO_noncandidate_files,IVT_pos_file=None,IVT_neg_file=None,batch_size=1000,epochs=100,x_valid=None,y_valid_ref=None,y_valid_call=None,val_mod_label=None,val_ratio_label=None,early_stop=0,bin=0.2,binsample=100,traintype="basecall",WT_label=1,binselecttype=1,steps_per_epoch=1000,sample_weights_flag=True,shuffleIVT=False,weight_exp=1,KO_noncandidate_weight=1,real_noncandidate_weight=1):
        print("debug training\n")
        print("steps_per_epoch="+str(steps_per_epoch)+"\n")
        if real_candidate_files is not None:
            real_candidate_files=real_candidate_files.replace("\"",'')
            real_noncandidate_files=real_noncandidate_files.replace("\"",'')
            KO_candidate_files=KO_candidate_files.replace("\"",'')
            KO_noncandidate_files=KO_noncandidate_files.replace("\"",'')
        
        print("samplesize:"+str(samplesize)) #use the smallest dataset (real data) as sample size!
        if x_valid is not None:
           #print("Monitor on validation data\n")
           val_prefix = "val_"
        else:
           val_prefix=""
        
        if early_stop!=0:
           if traintype=="basecall":
              early_stopping = EarlyStopping(monitor=val_prefix+'loss', patience=early_stop)
           elif traintype == "totalloss":
                early_stopping = EarlyStopping(monitor=val_prefix+'loss', patience=early_stop)
           elif traintype =="ratioloss":
              early_stopping = EarlyStopping(monitor=val_prefix+'mod_pred_loss', patience=early_stop)
           elif traintype == "modloss":
               early_stopping = EarlyStopping(monitor=val_prefix+'mod_pred_loss', patience=early_stop)
        
        if traintype=="basecall":
               best_model_path = os.path.join(OUTPATH,'model_basecall_best.ckpt')
        elif traintype == "totalloss":
                best_model_path = os.path.join(OUTPATH,'model_best.ckpt')
        else:   
               best_model_path = os.path.join(OUTPATH,'model_best.ckpt')
        
        if traintype=="basecall":
                model_checkpoint = ModelCheckpoint(best_model_path, monitor=val_prefix+'loss', save_best_only=True, verbose=1)
        elif traintype=="totalloss":
                model_checkpoint = ModelCheckpoint(best_model_path, monitor=val_prefix+'loss', save_best_only=True, verbose=1)
        elif traintype=="ratioloss":
                model_checkpoint = ModelCheckpoint(best_model_path, monitor=val_prefix+'mod_pred_loss', save_best_only=True, verbose=1)
        elif traintype == "modloss":
                model_checkpoint = ModelCheckpoint(best_model_path, monitor=val_prefix+'mod_pred_loss', save_best_only=True, verbose=1)
        
        model_checkpoint_everystep = ModelCheckpoint(os.path.join(OUTPATH,'model_current.ckpt'), save_best_only=False, verbose=1)
        
        if IVT_pos_file is None or not os.path.exists(IVT_pos_file):
            generator_fun = generator_realdata
        else:
            if real_candidate_files is None or not os.path.exists(real_candidate_files.split(',')[0]):
                 if shuffleIVT:
                       generator_fun=generator_IVT_shuffle #generator_IVT
                 else:
                       generator_fun=generator_IVT
            else:
                 generator_fun = generator_all #generator_IVT_data#generator_IVT_data_nonshuffle
        
        if x_valid is None:
           #print("no val")
           if early_stop!=0:
               history=self.model.fit_generator(generator_fun(IVT_pos_file,IVT_neg_file,real_candidate_files,real_noncandidate_files,KO_candidate_files,KO_noncandidate_files,batch_size=batch_size,bin=bin,binsample=binsample,WT_label=WT_label,binselecttype=binselecttype,sample_weights_flag=sample_weights_flag,featuredim=self.featuredim,windowsize=self.max_len,weight_exp=weight_exp,KO_noncandidate_weight=KO_noncandidate_weight,real_noncandidate_weight=real_noncandidate_weight), steps_per_epoch=steps_per_epoch, epochs=epochs, verbose=1,
                   callbacks=[model_checkpoint,early_stopping,model_checkpoint_everystep], shuffle=True)
           else:
               history=self.model.fit_generator(generator_fun(IVT_pos_file,IVT_neg_file,real_candidate_files,real_noncandidate_files,KO_candidate_files,KO_noncandidate_files,batch_size=batch_size,bin=bin,binsample=binsample,WT_label=WT_label,binselecttype=binselecttype,sample_weights_flag=sample_weights_flag,featuredim=self.featuredim,windowsize=self.max_len,weight_exp=weight_exp,KO_noncandidate_weight=KO_noncandidate_weight,real_noncandidate_weight=real_noncandidate_weight), steps_per_epoch=steps_per_epoch, epochs=epochs, verbose=1,
                   callbacks=[model_checkpoint,model_checkpoint_everystep], shuffle=True)
        elif early_stop!=0:
            #print("val and early_stop")
            #print("x_valid shape="+str(x_valid.shape))
            #print("y_valid_ref shape="+str(y_valid_ref.shape))
            val_fake_input = np.zeros([len(x_valid),self.maxoutputlen,1])
            history=self.model.fit_generator(generator_fun(IVT_pos_file,IVT_neg_file,real_candidate_files,real_noncandidate_files,KO_candidate_files,KO_noncandidate_files,batch_size=batch_size,bin=bin,binsample=binsample,WT_label=WT_label,binselecttype=binselecttype,sample_weights_flag=sample_weights_flag,featuredim=self.featuredim,windowsize=self.max_len,weight_exp=weight_exp,KO_noncandidate_weight=KO_noncandidate_weight,real_noncandidate_weight=real_noncandidate_weight), steps_per_epoch=steps_per_epoch, epochs=epochs, verbose=1,
                        callbacks=[model_checkpoint,early_stopping,model_checkpoint_everystep], shuffle=True,\
                        validation_data=[[x_valid,val_fake_input,val_fake_input],[np.expand_dims(y_valid_ref,axis=-1),np.expand_dims(y_valid_call,axis=-1),val_mod_label,np.expand_dims(val_ratio_label,-1)]])
        else:
            print("val")
            val_fake_input = np.zeros([len(x_valid),self.maxoutputlen,1])
            history=self.model.fit_generator(generator_fun(IVT_pos_file,IVT_neg_file,real_candidate_files,real_noncandidate_files,KO_candidate_files,KO_noncandidate_files,batch_size=batch_size,bin=bin,binsample=binsample,WT_label=WT_label,binselecttype=binselecttype,sample_weights_flag=sample_weights_flag,featuredim=self.featuredim,windowsize=self.max_len,weight_exp=weight_exp,KO_noncandidate_weight=KO_noncandidate_weight,real_noncandidate_weight=real_noncandidate_weight), steps_per_epoch=steps_per_epoch, epochs=epochs, verbose=1,
                        callbacks=[model_checkpoint,model_checkpoint_everystep], shuffle=True,\
                        validation_data=[[x_valid,val_fake_input,val_fake_input],[np.expand_dims(y_valid_ref,axis=-1),np.expand_dims(y_valid_call,axis=-1),val_mod_label,np.expand_dims(val_ratio_label,-1)]])
        
        hist = history.history
        if traintype=="basecall":
             self.model.save_weights(os.path.join(OUTPATH,"basecall_final.ckpt"))
             print("Model %s saved." % (os.path.join(OUTPATH,"basecall_final.ckpt")))
        else:
             self.model.save_weights(os.path.join(OUTPATH,"final.ckpt"))
             print("Model %s saved." % (os.path.join(OUTPATH,"final.ckpt")))
        
        self.model.load_weights(best_model_path)
        print("Model training finished!\n")
