#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --account=kinfai1 --mem=90GB
#

source activate /nfs/turbo/umms-kinfai/haorli/mambaforge/envs/ReDD_CPU
echo $message
echo $weights_dir_step2
echo "steps_per_epoch"
echo $steps_per_epoch
echo "auto_steps_per_epoch?"
echo $auto_steps_per_epoch
echo $batch_size
echo $bin_size
echo $binselecttype
echo $windowsize
echo $featuredim
echo $outputdir

if [ "$outputdir" = '' ];then
outputdir="./"
fi 

if [ "$windowsize" = '' ];then
windowsize=9
fi

if [ "$KO_noncandidate_weight" = '' ];then
KO_noncandidate_weight=1
fi

if [ "$real_noncandidate_weight" = '' ];then
real_noncandidate_weight=1
fi

if [ "$featuredim" = '' ];then
featuredim=5
fi

if [ "$early_stop1" = '' ];then
early_stop1=3
fi

if [ "$early_stop2" = '' ];then
early_stop2=3
fi

if [ "$epochs1" = '' ];then
epochs1=100
fi

if [ "$epochs2" = '' ];then
epochs2=100
fi

if [ "$IVT_pos_file" = '' ];then
IVT_pos_file="None"
fi

if [ "$IVT_neg_file" = '' ];then
IVT_neg_file="None"
fi


if [ "$auto_steps_per_epoch" = '' ];then
   if [ "$weights_dir_step2" = '' ];then
   python /scratch/kinfai_root/kinfai0/haorli/software/ReDD/scripts/Train.py --early_stop1 ${early_stop1} --epochs1 ${epochs1} --epochs2 ${epochs2} --early_stop2 ${early_stop2} --dim_attention 150 --featuredim $featuredim --dim_lstm 150 --embedding_output_dim 30 --windowsize $windowsize --drop_input 0 --drop_cnn 0.25 --drop_flat 0.2 \
   --message ${message} \
   -outputdir ${outputdir} \
   --bin_size ${bin_size} \
   -steps_per_epoch ${steps_per_epoch} \
   --IVT_pos_file ${IVT_pos_file} --IVT_neg_file ${IVT_neg_file} \
   --real_candidate_files  ${real_candidate_files} \
   --KO_candidate_files ${KO_candidate_files} \
   --real_noncandidate_files ${real_noncandidate_files} \
   --KO_noncandidate_files ${KO_noncandidate_files} \
   --batch_size ${batch_size} \
   -real_noncandidate_weight ${real_noncandidate_weight} -KO_noncandidate_weight ${KO_noncandidate_weight}
   
   else
   python /scratch/kinfai_root/kinfai0/haorli/software/ReDD/scripts/Train.py --load_pretrain --early_stop1 ${early_stop1} --epochs1 ${epochs1} --epochs2 ${epochs2} --early_stop2 ${early_stop2} --dim_attention 150 --featuredim $featuredim --dim_lstm 150 --embedding_output_dim 30 --windowsize $windowsize --drop_input 0 --drop_cnn 0.25 --drop_flat 0.2 \
   --message ${message} \
   -outputdir ${outputdir} \
   --bin_size ${bin_size} \
   -steps_per_epoch ${steps_per_epoch} \
   --IVT_pos_file ${IVT_pos_file} --IVT_neg_file ${IVT_neg_file} \
   --real_candidate_files  ${real_candidate_files} \
   --KO_candidate_files ${KO_candidate_files} \
   --real_noncandidate_files ${real_noncandidate_files} \
   --KO_noncandidate_files ${KO_noncandidate_files} \
   --weights_dir_step2 ${weights_dir_step2} \
   --batch_size ${batch_size} \
   -real_noncandidate_weight ${real_noncandidate_weight} -KO_noncandidate_weight ${KO_noncandidate_weight}
   fi
else
     if [ "$weights_dir_step2" = '' ];then
      python /scratch/kinfai_root/kinfai0/haorli/software/ReDD/scripts/Train.py -auto_steps_per_epoch --binselecttype $binselecttype --sample_weights_flag --early_stop1 3 --epochs1 100 --epochs2 100 --early_stop2 3 --dim_attention 150 --featuredim $featuredim --dim_lstm 150 --embedding_output_dim 30 --windowsize $windowsize --drop_input 0 --drop_cnn 0.25 --drop_flat 0.2 \
      --message ${message} \
      -outputdir ${outputdir} \
      --bin_size ${bin_size} \
      -steps_per_epoch ${steps_per_epoch} \
      --IVT_pos_file ${IVT_pos_file} --IVT_neg_file ${IVT_neg_file} \
      --real_candidate_files  ${real_candidate_files} \
      --KO_candidate_files ${KO_candidate_files} \
      --real_noncandidate_files ${real_noncandidate_files} \
      --KO_noncandidate_files ${KO_noncandidate_files} \
      --batch_size ${batch_size} \
      -real_noncandidate_weight ${real_noncandidate_weight} -KO_noncandidate_weight ${KO_noncandidate_weight}
      
      else
      python /scratch/kinfai_root/kinfai0/haorli/software/ReDD/scripts/Train.py -auto_steps_per_epoch --binselecttype $binselecttype --load_pretrain --sample_weights_flag --early_stop1 3 --epochs1 100 --epochs2 100 --early_stop2 3 --dim_attention 150 --featuredim $featuredim --dim_lstm 150 --embedding_output_dim 30 --windowsize $windowsize --drop_input 0 --drop_cnn 0.25 --drop_flat 0.2 \
      --message ${message} \
      -outputdir ${outputdir} \
      --bin_size ${bin_size} \
      -steps_per_epoch ${steps_per_epoch} \
      --IVT_pos_file ${IVT_pos_file} --IVT_neg_file ${IVT_neg_file} \
      --real_candidate_files  ${real_candidate_files} \
      --KO_candidate_files ${KO_candidate_files} \
      --real_noncandidate_files ${real_noncandidate_files} \
      --KO_noncandidate_files ${KO_noncandidate_files} \
      --weights_dir_step2 ${weights_dir_step2} \
      --batch_size ${batch_size} \
      -real_noncandidate_weight ${real_noncandidate_weight} -KO_noncandidate_weight ${KO_noncandidate_weight}
      fi

fi
