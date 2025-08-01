import numpy as np
from sklearn.cluster import SpectralClustering,KMeans
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import matthews_corrcoef,accuracy_score,precision_score,recall_score,confusion_matrix
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
hash_candidat={}
hash_candidat['AFG-H1_directRNA']='H1-AFG.candidate_sites.tab'
hash_candidat['AFG-H9_directRNA']='H9-AFG.candidate_sites.tab'
hash_candidat["PGC-H1_directRNA"]='H1-PGC.candidate_sites.tab'
hash_candidat["DE-H1_directRNA"]='H1-DE.candidate_sites.tab'
hash_candidat["DE-H9_directRNA"]='H9-DE.candidate_sites.tab'
hash_candidat["GM12878_directRNA"]='GM12878.candidate_sites.tab'
hash_candidat["H1-hESC_directRNA"]='H1-hESC.candidate_sites.tab'
hash_candidat["H9-hESC_directRNA"]='H9-hESC.candidate_sites.tab'
hash_candidat['HEK293T_DKO_directRNA']='HEK293T_WT.candidate_sites.tab'
hash_candidat["HEK293T_WT_directRNA"]='HEK293T_WT.candidate_sites.tab'
hash_candidat["HEK_WT_pass"]='HEK293T_WT.candidate_sites.tab'

rediportal = "/users/PCON0009/duolin/REDIportal/hg38/REDIportal_hg38.txt"
REDIportal_sites={}
input = open(rediportal,'r')
for line in input:
        chr_ = line.split()[0]
        pos_ = line.split()[1]
        REDIportal_sites[chr_+"-"+pos_]=1



cutoff = 0.5#0.5
featuredim=5
windowsize=9
shorread_coverage_cutoff=30#30#30 #100
longread_coverage_cutoff=0#6
use_coverage_cutoff=1

for datatype in ['HEK293T_WT_directRNA']:
    modelname="hg38retrain_HEK293T_KO1_chrsplitHEK293T_WT_directRNA_epochs1_100_epochs2_88"
    filename="Eventfeature_"+str(featuredim)+"feature_win"+str(windowsize)+"/"+modelname+"/test_"+datatype+"_onlycandidate.txt"
    if datatype=="HCT116_WT":
        candidatefile="/fs/ess/scratch/PCON0009/duolin/real_data/Dinopore_data/HCT116_hg38.bed"
    else:
         candidatefile="/fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/7-reditools2_final/"+hash_candidat[datatype]
    
    print(datatype)
    input = open(filename)
    pos_coverage={}
    coverage = {}
    cutoff = 0.5#0.5
    weighted_predict_score={}
    for line in input:
        score=float(line.split("\t")[-1])
        if score>=cutoff:
           predict_label=1
        else:
           predict_label=0
        
        transid= line.split("\t")[2]
        transpos = line.split("\t")[3]
        chrpos = transid+"-"+transpos
        if chrpos not in weighted_predict_score:
            weighted_predict_score[chrpos]=score
        else:
            weighted_predict_score[chrpos]+=score
        
        if chrpos not in pos_coverage.keys():
            pos_coverage[chrpos]=predict_label
            coverage[chrpos]=1
        else:
            pos_coverage[chrpos]+=predict_label
            coverage[chrpos]+=1
    
    AG_ratio_per_site={}
    shortreadcoverage={}
    input = open(candidatefile,'r')
    for line in input:
        chr_ = line.split()[0]
        pos_ = line.split()[1]
        chrpos=chr_+"-"+pos_
        AG_ratio_per_site[chrpos]=float(line.split("\t")[3])
        if datatype !="HCT116_WT":
           shortreadcoverage[chrpos]=float(line.split("\t")[4])
    
    pos_coverage_list=[]
    coverage_list=[]
    select_site=[]
    for site in AG_ratio_per_site:
      if AG_ratio_per_site[site]>=0.95:
         if use_coverage_cutoff==1:
           if datatype !="HCT116_WT":
              if shortreadcoverage[site]<shorread_coverage_cutoff:
                 continue
         
         if site in pos_coverage and coverage[site]>=longread_coverage_cutoff:
              print(site)
              pos_coverage_list.append(pos_coverage[site])
              coverage_list.append(coverage[site])
              select_site.append(site)
    
    print(pos_coverage_list)
    print(coverage_list)
    pos_read_num=np.sum(pos_coverage_list)
    coverage_read=np.sum(coverage_list)
    recall =  pos_read_num/coverage_read
    print("recall:"+str(recall)+"\n")
    filename="Eventfeature_"+str(featuredim)+"feature_win"+str(windowsize)+"/"+modelname+"/test_"+datatype+"_noncandidate.txt"
    input = open(filename)
    neg_coverage={}
    KO_coverage = {}
    num=0
    max_read=coverage_read #100000
    KO_coverage_list=[]
    neg_read_list=[]
    for line in input:
        if num>=coverage_read*100:
            break
        
        num+=1
        score=float(line.split("\t")[-1])
        if score>=cutoff:
           predict_label=1
        else:
           predict_label=0
        
        transid= line.split("\t")[2]
        transpos = line.split("\t")[3]
        chrpos = transid+"-"+transpos
        if chrpos in REDIportal_sites:
            continue
        
        KO_coverage_list.append(1)
        neg_read_list.append(predict_label)
    
    KO_coverage_read=0
    neg_read_num=0
    neg_read_num=np.sum(neg_read_list[:coverage_read])
    KO_coverage_read=np.sum(KO_coverage_list[:coverage_read])
    print("sitenum="+str(len(coverage_list)))
    FPR=neg_read_num/KO_coverage_read
    print("FPR="+str(FPR)+"\n")
    print("confusion_matrix:"+str(pos_read_num)+","+str(coverage_read-pos_read_num)+","+str(neg_read_num)+","+str(KO_coverage_read-neg_read_num)+"\n")
    print("confusion_matrix_acc:"+str(pos_read_num/coverage_read)+","+str((coverage_read-pos_read_num)/coverage_read)+","+str(neg_read_num/KO_coverage_read)+","+str((KO_coverage_read-neg_read_num)/KO_coverage_read)+"\n")
    accuracy=(pos_read_num+KO_coverage_read-neg_read_num)/(coverage_read+KO_coverage_read)
    print("accuracy="+str(accuracy)+"\n")