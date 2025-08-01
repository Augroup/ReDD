
import numpy as np
from sklearn.cluster import SpectralClustering,KMeans
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import matthews_corrcoef,accuracy_score,precision_score,recall_score,confusion_matrix
from scipy.stats import pearsonr

import matplotlib.pyplot as plt

def evaluate(predict_proba,true_label,predict_label = None,cutoff=0.5):
    from sklearn.metrics import roc_auc_score,average_precision_score,matthews_corrcoef
    roc_score=roc_auc_score(true_label,predict_proba)
    pr_score=average_precision_score(true_label,predict_proba)
    if predict_label is None:
          predict_label = [np.argmax([cutoff,x]) for x in predict_proba]
    
    mcc = matthews_corrcoef(true_label,predict_label)
    recall = recall_score(true_label,predict_label)
    precision = precision_score(true_label,predict_label)
    acc = accuracy_score(true_label,predict_label)
    tn, fp, fn, tp = confusion_matrix(true_label,predict_label).ravel()
    FPR = fp/(tn+fp)
    sensitivity = tp/(fn+tp)
    return roc_score,pr_score,mcc,recall,precision,acc,FPR,sensitivity


name_mapping={}
name_mapping['AFG-H1_directRNA']='H1-AFG'
name_mapping['AFG-H9_directRNA']='H9-AFG'
name_mapping["PGC-H1_directRNA"]='H1-PGC'
name_mapping["DE-H1_directRNA"]='H1-DE'
name_mapping["DE-H9_directRNA"]='H9-DE'
name_mapping["GM12878_directRNA"]='GM12878'
name_mapping["H1-hESC_directRNA"]='H1-hESC'
name_mapping["H9-hESC_directRNA"]='H9-hESC'
name_mapping['HEK293T_DKO_directRNA']='HEK293T'
name_mapping["HEK293T_WT_directRNA"]='HEK293T'

# datatype ='GM12878_directRNA' 'H1-hESC_directRNA' 'AFG-H1_directRNA' 'DE-H1_directRNA' 'H9-hESC_directRNA' 'AFG-H9_directRNA' 'DE-H9_directRNA' 'HEK293T_WT_directRNA' 'PGC-H1_directRNA'

#datatype='AFG-H1_directRNA'
#datatype='H1-hESC_directRNA'
#datatype='PGC-H1_directRNA'
#datatype='DE-H1_directRNA'
#datatype='DE-H9_directRNA'
#datatype='H9-hESC_directRNA'
#datatype='AFG-H9_directRNA'
#datatype='GM12878_directRNA'
#datatype='HEK293T_DKO_directRNA'
#datatype='HEK293T_WT_directRNA'
#datatype="HEK_WT_pass"
#filename="Eventfeature_5feature_win9/hg38_merge9alldata5_noearlystop_ep40/"+datatype+"_onlycandidate.txt"
#filename="Eventfeature_5feature_win9/hg38_merge9alldata5_noearlystop_ep40_finalpredict/"+datatype+"_onlycandidate.txt"
#filename="Eventfeature_5feature_win9/hg38_merge9alldata5_del_H1_negsite/"+datatype+"_onlycandidate.txt"
#filename="Eventfeature_5feature_win9/hg38_HEK293T_finalpredict/"+datatype+"_onlycandidate.txt"
#incremental
#datatype='H1-hESC_directRNA'
#datatype='AFG-H9_directRNA'
#datatype="DE-H9_directRNA"

#filename="Eventfeature_5feature_win9/hg38_incre5train__HEK293T_WT_directRNA/"+datatype+"_onlycandidate.txt"
#filename="Eventfeature_5feature_win9/hg38_incre5train__HEK293T_WT_directRNA_AFG-H1_directRNA/"+datatype+"_onlycandidate.txt"
#filename="Eventfeature_5feature_win9/hg38_incre5train__HEK293T_WT_directRNA_AFG-H1_directRNA_DE-H1_directRNA/"+datatype+"_onlycandidate.txt"

#filename="Eventfeature_5feature_win9/hg38_incre5train__HEK293T_WT_directRNA_H9-hESC_directRNA/"+datatype+"_onlycandidate.txt"
#filename="Eventfeature_5feature_win9/hg38_incre5train__HEK293T_WT_directRNA_AFG-H1_directRNA_DE-H1_directRNA_H9-hESC_directRNA/"+datatype+"_onlycandidate.txt"
#modelname="hg38_merge9alldata5_del_H9_negsite"

#modelname="hg38_HEK293T_entropyratio_weight1000_finalpredict"
#modelname="hg38_HEK293T_finalpredict"
#modelname="hg38_HEK293T_entropyratio_weight1000_sampleweight" #haven't finished
#modelname="hg38_HEK293T_entropyratio_mergeposneg_weight1000"
#modelname="hg38_HEK293T_entropyratio_weight0"
#modelname="hg38_HEK293T_entropyratio_weight1"
#modelname="hg38_HEK293T_entropyratio_weight0_sampleweight"
#modelname="hg38_HEK293T_entropyratio_weight1000_sampleweight"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_ep49"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_GM12878_directRNA_epochs5_lossflag1_ratiolossweight1000"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_GM12878_directRNA_epochs5"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_GM12878_directRNA_epochs5_lossflag1_ratiolossweight1_sampleweight_wexp2"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_GM12878_directRNA_epochs5_lossflag1_ratiolossweight1_sampleweight"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_"+datatype+"_epochs5"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_"+datatype+"_epochs5_lossflag2_ratiolossweight1"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_GM12878_directRNA_epochs5_lossflag1_ratiolossweight1_sampleweight_2wexp2"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_"+datatype+"_epochs5_lossflag3_ratiolossweight1_sampleweight_2wexp2"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_"+datatype+"_epochs5_lossflag4_ratiolossweight1_sampleweight_2wexp2"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_"+datatype+"_epochs5_lossflag2_ratiolossweight1_sampleweight_2wexp2"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_"+datatype+"_epochs5_lossflag4_ratiolossweight1_sampleweight_2wexp2"
#modelname="hg38_merge9alldata5_del_H9_negsite_epochs30_lossflag4_ratiolossweight1_sampleweight_2wexp2"
#modelname="hg38_merge9alldata5_del_H9_epochs30_lossflag4_ratiolossweight1_sampleweight_2wexp2"
#modelname='hg38_merge9alldata5_del_H1_negsite_epochs30_lossflag4_ratiolossweight1_sampleweight_2wexp2'
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_H1-hESC_directRNA_epochs200_persitemodel"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_H1-hESC_directRNA_epochs20_persitemodel_sampleratio"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_HEK293T_WT_directRNA_epochs5_selectallcandidate"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_HEK293T_WT_directRNA_epochs5_selectallcandidate_shufeachratio"
#modelname="hg38_merge9alldata5_del_H1_negsite"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_HEK293T_WT_directRNA_epochs5_KOnonw5"
#modelname="hg38_merge9alldata5_del_H9_negsite_KOnonw5_withLSTM"
#modelname="hg38_merge9alldata5_del_H1_negsite_KOnonw5_withLSTM"
modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_mergedH_epochs35_KOnonw5withLSTM"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_HEK293T_WT_directRNA_epochs5_KOnonw5_withLSTM"
#modelname="hg38_merge9alldata5_noearlystop_run1ep40_run2ep60_negsite_GM12878_directRNA_epochs5_KOnonw5_withLSTM"
#modelname="hg38_HEK293T_KOnonw5_withLSTM"
#for datatype in ["DE-H1_directRNA",'H1-hESC_directRNA','AFG-H1_directRNA','PGC-H1_directRNA']:
#for datatype in ["DE-H9_directRNA",'H9-hESC_directRNA','AFG-H9_directRNA']:
#for datatype in ["PGC-H1_directRNA","GM12878_directRNA","HEK293T_WT_directRNA"]:

modelname="hg38retrain_HEK293T_chrsplitHEK293T_WT_directRNA_epochs1_100_epochs2_fist35_sec25" #in release_code
modelname="hg38retrain_HEK293T_chrsplitHEK293T_WT_directRNA_epochs1_100_epochs2_fist50_sec13"
#for datatype in ['H9-hESC_directRNA','AFG-H9_directRNA','DE-H9_directRNA','H1-hESC_directRNA','AFG-H1_directRNA','DE-H1_directRNA','PGC-H1_directRNA','GM12878_directRNA','HEK293T_WT_directRNA']:

for datatype in ['HEK293T_WT_directRNA','HEK293T_DKO_directRNA']:
    outputfolder="Eventfeature_5feature_win9/"+modelname+"/"
    #filename="Eventfeature_5feature_win9/"+modelname+"/"+datatype+"_onlycandidate.txt"
    filename="Eventfeature_5feature_win9/"+modelname+"/test_"+datatype+"_onlycandidate.txt"
    try:
       input = open(filename)
    except:
        print("skipping "+datatype)
        continue
    
    pos_coverage={}
    coverage = {}
    cutoff = 0.5
    for line in input:
        if float(line.split("\t")[-1])>cutoff:
           predict_label=1
        else:
           predict_label=0
        
        transid= line.split("\t")[2]
        transpos = line.split("\t")[3]
        chrpos = transid+"-"+transpos
        
        if chrpos not in pos_coverage.keys():
            pos_coverage[chrpos]=predict_label
            coverage[chrpos]=1
        else:
            pos_coverage[chrpos]+=predict_label
            coverage[chrpos]+=1
    
    coverage_cutoff=5
    use_coverage_cutoff=1
    pos_coverage_ratio = 0 #0.2 #0.2#0.2
    label=[]
    predict =[]
    AG_ratio = []
    predict_label=[]
    WT_index=[]
    KO_index=[]
    input = open(filename)
    index=0
    for line in input:
        transid= line.split("\t")[2]
        transpos = line.split("\t")[3]
        chrpos = transid+"-"+transpos
        if use_coverage_cutoff==1:
           if coverage[chrpos]<=coverage_cutoff:
                continue
        
        if float(line.split("\t")[-1])>cutoff and (pos_coverage[chrpos]/coverage[chrpos])>pos_coverage_ratio:
           predict_label.append(1)
        else:
           predict_label.append(0)
        
        labeltype = line.split("\t")[0].split(":")[0]
        if ":" in line:
             ratio = float(line.split("\t")[0].split(":")[1])
        else:
             ratio = 0
        
        if "KO" in labeltype:
           KO_index.append(index)
           #ratio = 0
        else:
           WT_index.append(index)
        
        if "KO" not in labeltype: #0.8: #0 only KO 0, 0.8 ratio>0.8 1
            label.append(1)
        else:
            label.append(0)
        
        
        AG_ratio.append(ratio)
        predict.append(float(line.split("\t")[-1]))
        index+=1
    
    WT_index = np.asarray(WT_index)
    KO_index = np.asarray(KO_index)
    
    label = np.asarray(label)
    predict = np.asarray(predict)
    predict_label = np.asarray(predict_label)
    AG_ratio =np.asarray(AG_ratio)
    
    
    #break
    
    #####random start from 0-1
    predict_ratio = []
    weighted_AG_ratios = []
    np.random.seed(seed=1234)
    ratio_list = np.random.uniform(0.0, 1.0, 100)
    KO_FPR = []
    weighted_AG_ratios_FDR = []
    sample_size=[]
    weighted_AG_ratio_eval = []
    MCClist = []
    FPRlist = []
    precisionlist = []
    recalllist = []
    sensitivitylist=[]
    for index in range(len(ratio_list)):
        start = ratio_list[index]
        end = start+0.1
        index_start_end = np.intersect1d(np.where((AG_ratio<=end) &(AG_ratio>start))[0],WT_index)
        if len(index_start_end)>=100:
           #index_start_end=index_start_end[:100]
           weighted_AG_ratio = np.mean(AG_ratio[index_start_end])
           weighted_AG_ratios.append(weighted_AG_ratio)
           ratio = len(np.where(predict_label[index_start_end]==1)[0])/len(predict_label[index_start_end]) 
           predict_ratio.append(ratio)
           sample_size.append(len(predict_label[index_start_end]))
           print("predict modification ratio in WT between AG ratio:%2.1f< <=%2.1f:"%(start,end))
           print("%5.2f=%5.2f/%5.2f\t%5.2f"%(ratio,len(np.where(predict_label[index_start_end]==1)[0]),len(predict_label[index_start_end]),weighted_AG_ratio))
    
    
    plt.clf()
    f, (ax1,ax2) = plt.subplots(2, 1,figsize=(6, 6),gridspec_kw={'height_ratios': [0.8, 2]})
    #ax2.bar(weighted_AG_ratios,sample_size,width=0.1,align='edge')
    ax1.set(xlim=(0, 1))
    ax1.set(ylim=(0, 20))
    ax1.set_yticks([5,10,15,20])
    
    ax1.set_title(name_mapping[datatype])
    #ax1.vlines(weighted_AG_ratios, ymin=np.zeros(len(weighted_AG_ratios)),ymax=[np.log2(x) for x in sample_size])
    #ax1.set_ylabel('log2(sample size)')
    ax1.vlines(weighted_AG_ratios, ymin=np.zeros(len(weighted_AG_ratios)), ymax = [np.log2(x) for x in sample_size],color="black")
    ax1.set_ylabel('Sample size \n(log2)')
    
    ax1.get_xaxis().set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['bottom'].set_visible(False)
    
    ax2.scatter(weighted_AG_ratios, predict_ratio,5,color="black")
    ax2.plot([0, 1], [0, 1], ls="--",color="black")
    hi_coverage_site_num= len(np.where(np.asarray(list(coverage.values()))>coverage_cutoff)[0])
    #ax2.set_xlabel('weighted_AG_ratio\nTotal number of candidate sites='+str(len(coverage))+"(>"+str(coverage_cutoff)+":"+str(hi_coverage_site_num)+")")
    pcc=pearsonr(predict_ratio,weighted_AG_ratios)[0]
    mse=(np.square(np.asarray(predict_ratio)-np.asarray(weighted_AG_ratios))).mean()
    ax2.set_xlabel('Expected A-to-I ratio\n(PCC=%.4f mse=%.4f)'%(pcc,mse))
    ax2.set_ylabel('Estimated A-to-I ratio')
    ax2.set_xticks(np.arange(0,1.1,0.2))
    ax2.set_yticks(np.arange(0,1.1,0.2))
    ax2.set(xlim=(0, 1), ylim=(0, 1))
    plt.subplots_adjust(hspace=0.03,wspace=0.2,top=0.805,bottom=0.395,left=0.410,right=0.790)
    #plt.show()
    #plt.savefig(outputfolder+datatype+"_binlevel_ratio.png")
    plt.savefig(outputfolder+"test_"+datatype+"_binlevel_ratio.png")

            