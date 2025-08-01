
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


coverage_cutoff=6 #5 10
covfilename=5
balance_neg=True
def evaluate(predict_proba,true_label,predict_label = None,cutoff=0.5):
    from sklearn.metrics import roc_auc_score,average_precision_score,matthews_corrcoef
    roc_score=-1
    #roc_score=roc_auc_score(true_label,predict_proba)
    pr_score=-1
    #pr_score=average_precision_score(true_label,predict_proba)
    if predict_label is None:
          predict_label = [np.argmax([cutoff,x]) for x in predict_proba]
    
    
    #mcc = matthews_corrcoef(true_label,predict_label)
    mcc=-1
    recall = recall_score(true_label,predict_label)
    precision=-1
    #precision = precision_score(true_label,predict_label)
    acc=-1
    #acc = accuracy_score(true_label,predict_label)
    tn, fp, fn, tp = confusion_matrix(true_label,predict_label).ravel()
    FPR = fp/(tn+fp)
    #sensitivity = tp/(fn+tp)
    sensitivity=-1
    return roc_score,pr_score,mcc,recall,precision,acc,FPR,sensitivity

featuredim=5
windowsize=9
includesnp=True
#for datatype in ['H1-hESC_directRNA','AFG-H1_directRNA','DE-H1_directRNA','PGC-H1_directRNA','H9-hESC_directRNA','AFG-H9_directRNA','DE-H9_directRNA']: 
for datatype in ['HEK293T_WT_directRNA']:#,'
    modelname="hg38retrain_HEK293T_chrsplitHEK293T_WT_directRNA_epochs1_100_epochs2_fist50_sec13"
    modelname="hg38retrain_HEK293T_chrsplitHEK293T_WT_directRNA_epochs1_100_epochs2_fist50_sec50"
    modelname="hg38retrain_HEK293T_chrsplitHEK293T_WT_directRNA_epochs1_100_epochs2_fist100_sec50"
    #modelname="hg38retrain_HEK293T_chrsplitHEK293T_WT_directRNA_epochs1_100_epochs2_50"
    modelname="hg38retrain_HEK293T_KO1_chrsplitHEK293T_WT_directRNA_epochs1_100_epochs2_88"
    filename1="Eventfeature_"+str(featuredim)+"feature_win"+str(windowsize)+"/"+modelname+"/test_"+datatype+"_onlycandidate_cov"+str(covfilename)+"_ratio0_modcov0.sitelev.bed" #background must use no filter except coverage
    filename2="Eventfeature_"+str(featuredim)+"feature_win"+str(windowsize)+"/"+modelname+"/test_"+datatype+"_noncandidate_cov"+str(covfilename)+"_ratio0_modcov0.sitelev.bed"
    
    outputfolder="Eventfeature_"+str(featuredim)+"feature_win"+str(windowsize)+"/"+modelname+"/"
    
    candidatefile="/fs/project/PCON0009/Au-scratch2/ying/StemCell/RNAdirect/reditools2/7-reditools2_final/"+hash_candidat[datatype]
    negative_site = "/fs/ess/scratch/PCON0009/duolin/real_data/Dinopore/test_HEK293T_WT_negative_sites_rthred0.1.txt"
    snp_site = "/fs/ess/scratch/PCON0009/duolin/real_data/Dinopore/test_HEK293T_WT_snp_sites_rthred0.1.txt"
    positive_site = "/fs/ess/scratch/PCON0009/duolin/real_data/Dinopore/test_HEK293T_WT_positive_sites_rthred0.1.txt"
    using_mean_prediction_score=True
    if using_mean_prediction_score:
       score_index=-1
    else:
       score_index=-2
    
    true_candidate_ratios={}
    input = open(candidatefile,'r')
    for line in input:
        chr_ = line.split()[0]
        pos_ = line.split()[1].strip()
        true_candidate_ratios[chr_+"-"+pos_]=float(line.split("\t")[3])
    
    nonmodification_sites={}
    input = open(negative_site,'r')
    for line in input:
        chr_ = line.split("-")[0]
        pos_ = line.split("-")[1].strip()
        nonmodification_sites[chr_+"-"+pos_]=1
    
    if includesnp:
        input = open(snp_site,'r')
        for line in input:
            chr_ = line.split("-")[0]
            pos_ = line.split("-")[1].strip()
            nonmodification_sites[chr_+"-"+pos_]=1
    
    positive_sites={}
    input = open(positive_site,'r')
    for line in input:
        chr_ = line.split("-")[0]
        pos_ = line.split("-")[1].strip()
        positive_sites[chr_+"-"+pos_]=1
    
    
    true_value_all=[]
    predict_value_all=[]
    num_candidate=0
    index_neg=[]
    inputs=[open(filename1),open(filename2)]
    for input in inputs:
      for line in input:
        readcoverage=int(line.split()[3])
        if readcoverage<coverage_cutoff:
             continue
        
        chr_ = line.split()[0]
        pos_ = line.split()[1]
        if len(line.split("\t")) == 6:
            predictscore=float(line.split("\t")[score_index]) #last score is mean prediciton score
        else:
            predictscore=float(line.split("\t")[-1]) #for old result there is no mean prediction scores
        
        if chr_+"-"+pos_ in positive_sites:
             true_value_all.append(true_candidate_ratios[chr_+"-"+pos_])
        elif chr_+"-"+pos_ in nonmodification_sites:
             true_value_all.append(0)
             index_neg.append(len(true_value_all)-1)
        else:
             continue
        
        predict_value_all.append(predictscore)
    
    true_value_all=np.asarray(true_value_all)
    predict_value_all=np.asarray(predict_value_all)
    index_neg=np.asarray(index_neg)
    import matplotlib.pyplot as plt
    plt.clf()
    plt.cla()
    plt.figure(figsize=(6, 6))
    params = {
         'figure.figsize': (6, 6),
         'axes.labelsize': '12',
         'axes.titlesize':'12',
         'xtick.labelsize':'12',
         'ytick.labelsize':'12'}
    
    plt.rcParams.update(params)
    from sklearn.metrics import roc_curve, auc
    index=0
    colorlist = ['b','g','r','c','m','y','k','w']
    for pos_coverage_ratio in [0.1,0.2,0.3,0.5]:
        index_pos=np.where(true_value_all >= pos_coverage_ratio)[0]
        index_pos_len=len(index_pos)
        if includesnp:
            select_neg=index_neg[:2*index_pos_len]
        else:
            select_neg=index_neg[:index_pos_len]
        
        predict_value_selected=predict_value_all[np.concatenate([index_pos,select_neg])]
        true_value_selected=true_value_all[np.concatenate([index_pos,select_neg])]
        true_value_selected = (true_value_selected >= pos_coverage_ratio).astype(int) 
        fpr, tpr,_ = roc_curve(true_value_selected, predict_value_selected)
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, color=colorlist[index], label='cutoff=%0.1f,AUC=%0.4f' % (pos_coverage_ratio,roc_auc))
        index+=1
    
    plt.plot([0, 1], [0, 1], color='navy')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    #plt.title(str(name_mapping)+' ROC curve at different AGratio cutoffs\n'+str(num_noncandidate)+' noncandidate sites and '+str(num_candidate)+" candidate sites (read coverage>10) \n")
    plt.title(str(name_mapping[datatype]))
    plt.legend(loc="lower right")
    plt.subplots_adjust(hspace=0.2,wspace=0.2,top=0.605,bottom=0.210,left=0.435,right=0.900)
    #plt.show()
    plt.savefig(outputfolder+datatype+"_cov"+str(coverage_cutoff)+"_ROC_readcoverage"+str(coverage_cutoff)+"_includesnp"+str(includesnp)+".png")
    
    #####precision recall curve
    import matplotlib.pyplot as plt
    plt.clf()
    plt.cla()
    plt.figure(figsize=(6, 6))
    params = {
         'figure.figsize': (6, 6),
         'axes.labelsize': '12',
         'axes.titlesize':'12',
         'xtick.labelsize':'12',
         'ytick.labelsize':'12'}
    
    plt.rcParams.update(params)
    from sklearn.metrics import roc_curve, auc,average_precision_score,precision_recall_curve
    index=0
    colorlist = ['b','g','r','c','m','y','k','w']
    for pos_coverage_ratio in [0.1,0.2,0.3,0.5]:
        index_pos=np.where(true_value_all >= pos_coverage_ratio)[0]
        index_pos_len=len(index_pos)
        if includesnp:
            select_neg=index_neg[:2*index_pos_len]
        else:
            select_neg=index_neg[:index_pos_len]
        
        predict_value_selected=predict_value_all[np.concatenate([index_pos,select_neg])]
        true_value_selected=true_value_all[np.concatenate([index_pos,select_neg])]
        true_value_selected = (true_value_selected >= pos_coverage_ratio).astype(int) 
        precision, recall, thresholds = precision_recall_curve(true_value_selected, predict_value_selected)
        roc_auc = average_precision_score(true_value_selected, predict_value_selected, average="micro") #auc(precision, recall)
        plt.plot(recall, precision, color=colorlist[index], label='cutoff=%0.1f,AUC=%0.4f' % (pos_coverage_ratio,roc_auc))
        index+=1
    
    #plt.plot([0, 1], [0, 1], color='navy')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall (sensitivity)')
    plt.ylabel('Precision')
    #plt.title(str(name_mapping[datatype])+' ROC curve at different AGratio cutoffs\n'+str(num_noncandidate)+' noncandidate sites and '+str(num_candidate)+" candidate sites (read coverage>10) \n")
    plt.title(str(name_mapping[datatype]))
    plt.legend(loc="lower left")
    #plt.show()
    plt.subplots_adjust(hspace=0.2,wspace=0.2,top=0.605,bottom=0.210,left=0.435,right=0.900)
    plt.savefig(outputfolder+datatype+"_cov"+str(coverage_cutoff)+"_precision_recall_readcoverage"+str(coverage_cutoff)+"_includesnp"+str(includesnp)+".png")
    