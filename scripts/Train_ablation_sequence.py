'''
change pn alldata4:
datatype with more samples train more steps
dont use auto_steps_per_epoch

For HDF5 lock error, solved by :
export HDF5_USE_FILE_LOCKING='FALSE'
keras version 2.2.4 
tensorflow '1.12.0'
python 3.5.2

taiyaki https://github.com/nanoporetech/taiyaki#installing-system-prerequisites
'''
import datetime
import itertools
from collections import OrderedDict
import argparse
import os
import sys
from keras.utils.io_utils import HDF5Matrix
import h5py
from scipy.signal import gaussian
from scipy.ndimage.filters import maximum_filter1d

#os.environ["CUDA_VISIBLE_DEVICES"] = "0"
#sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import tensorflow as tf
from keras.backend.tensorflow_backend import set_session
import keras.backend as K
from model import *  #add noncanddiate site 

#batch_size = 1000 #batch size must < 148 because smallest ratio is 148 for 0.2 bin
nb_classes = 2 #P1 P2 N

# starts training in CNN model
def run_model(dataset, **kwargs):
    '''config CPU parallelization'''
    if kwargs['device'] == 'CPU':
       config = tf.ConfigProto(intra_op_parallelism_threads=kwargs['threads'],
                     inter_op_parallelism_threads=kwargs['threads'],
                     allow_soft_placement=True,
                     device_count={'CPU': kwargs['threads']})
       
       sess = tf.Session(config=config)
       set_session(session=sess)
    elif kwargs['device'] == 'GPU':
         gpu_options = tf.GPUOptions()
         gpu_options.allow_growth = True
         config = tf.ConfigProto(gpu_options=gpu_options)
         #config.gpu_options.per_process_gpu_memory_fraction=0.4
         sess = tf.Session(config=config)
         set_session(session=sess)
    
    featuredim = kwargs['featuredim']
    '''
    if not kwargs['predict']:
        val_x=None
        val_y_ref=None
        val_y_call=None
        val_mod_label=None
        val_ratio_label=None
        if kwargs['val_hdf5'] is not None:
             h5f = h5py.File(kwargs['val_hdf5'],"r+")
             val_x = h5f['X'][:]
             val_y_ref = h5f['y_ref'][:]
             val_y_call = h5f['y_call'][:]
             h5f.close()
        
        if kwargs['val_IVT_neg_file'] is not None:
             h5f_IVT_neg = h5py.File(kwargs['val_IVT_neg_file'],"r+")
             h5f_IVT_pos = h5py.File(kwargs['val_IVT_pos_file'],"r+")
             win_center_index=int(len(h5f_IVT_pos['X'][0])/2)
             nt = int(kwargs['windowsize']/2)
             featuredim=int(kwargs['featuredim'])
             val_x = np.concatenate([h5f_IVT_neg['X'][:,win_center_index-nt:win_center_index+nt+1,:featuredim],h5f_IVT_pos['X'][:,win_center_index-nt:win_center_index+nt+1,:featuredim]],axis = 0)
             val_y_ref = np.concatenate([h5f_IVT_neg['y_ref'][:],h5f_IVT_pos['y_ref'][:]],axis = 0)
             val_y_call = np.concatenate([h5f_IVT_neg['y_call'][:],h5f_IVT_pos['y_call'][:]],axis = 0)
             val_mod_label = np.concatenate([np.zeros(len(h5f_IVT_neg['X'])),np.ones(len(h5f_IVT_pos['X']))],axis=0)
             val_ratio_label = np.concatenate([np.zeros(len(h5f_IVT_neg['X'])),np.ones(len(h5f_IVT_pos['X']))],axis=0)
             h5f_IVT_neg.close()
             h5f_IVT_pos.close()
        
        #h5f = h5py.File(kwargs['real_file'],"r") ##use the smallest dataset (real data) as sample size!
        #for merge data
        #h5f.close()
    '''
    samplesize=-1 #4006966
    if kwargs['predict']:
         load_pretrain=True
    else:
         if kwargs['load_pretrain']:
            load_pretrain = True
         else:
            load_pretrain = False
    
    model = multihead_attention(kwargs['windowsize'], kwargs['maxoutputlen'],featuredim,nb_classes, OUTPATH, kfold_index=0)  # initialize
    if kwargs['predict']:
        model.build_model_seq2seq(dim_attention=kwargs['dim_attention'], #10 vs 150
                              load_weights=load_pretrain, weight_dir=kwargs['weights_dir'],
                              drop_input=kwargs['drop_input'],
                              drop_cnn=kwargs['drop_cnn'],
                              drop_flat=kwargs['drop_flat'],
                              BatchNorm=kwargs['BatchNorm'],
                              dim_lstm=kwargs['dim_lstm'],
                              embedding_output_dim=kwargs['embedding_output_dim'],
                              is_monotonic=kwargs['is_monotonic'],
                              Att_regularizer_weight=kwargs['Att_regularizer_weight'],
                              headnum=kwargs['headnum'],
                              W_regularizer=kwargs['W_regularizer'],
                              fc_dim = kwargs['fc_dim'],
                              fcnum = kwargs['fcnum']
                              )
    else: #training
       steps_per_epoch=kwargs['steps_per_epoch']
       if kwargs['auto_steps_per_epoch']:
          if  os.path.exists(kwargs['IVT_pos_file']):
             IVT_pos_h5fs=h5py.File(kwargs['IVT_pos_file'],"r")
             steps_per_epoch=int(np.ceil(float(len(IVT_pos_h5fs['info']))/float(kwargs['batch_size']))) 
          else:
            steps_per_epoch=kwargs['bin_sample']
            print("auto steps_per_epoch="+str(steps_per_epoch)+"\n")
       
       if kwargs['train_type']!='total_loss' and kwargs['weights_dir_step2'] is None:
           print("Training basecalling procedure....\n")
           
           model.build_model_seq2seq(dim_attention=kwargs['dim_attention'], #10 vs 150
                                  load_weights=load_pretrain, weight_dir=kwargs['weights_dir'],
                                  drop_input=kwargs['drop_input'],
                                  drop_cnn=kwargs['drop_cnn'],
                                  drop_flat=kwargs['drop_flat'],
                                  BatchNorm=kwargs['BatchNorm'],
                                  dim_lstm=kwargs['dim_lstm'],
                                  embedding_output_dim=kwargs['embedding_output_dim'],
                                  is_monotonic=kwargs['is_monotonic'],
                                  Att_regularizer_weight=kwargs['Att_regularizer_weight'],
                                  headnum=kwargs['headnum'],
                                  W_regularizer=kwargs['W_regularizer'],
                                  fc_dim = kwargs['fc_dim'],
                                  fcnum = kwargs['fcnum'],
                                  loss_weights=[1., 1.,0.,0.], #train model focus on base calling prediction
                                  )
           
           print("before basecalling training\n")
           
           model.train(samplesize,kwargs['real_candidate_files'],kwargs['real_noncandidate_files'],
                       kwargs['KO_candidate_files'],
                       kwargs['KO_noncandidate_files'],
                       IVT_pos_file=kwargs['IVT_pos_file'],
                       IVT_neg_file=kwargs['IVT_neg_file'],
                       batch_size=kwargs['batch_size'], epochs=kwargs['epochs1'],
                       #no validation for basecalling training x_valid=val_x,y_valid_ref=val_y_ref,y_valid_call=val_y_call,val_mod_label=val_mod_label,val_ratio_label=val_ratio_label,
                       early_stop=kwargs['early_stop1'],bin = kwargs['bin_size'],binsample = kwargs['bin_sample'],
                       traintype="basecall",WT_label=kwargs['WT_label'],
                       sample_weights_flag=kwargs['sample_weights_flag'],
                       steps_per_epoch=steps_per_epoch,
                       shuffleIVT=kwargs['shuffleIVT']
                       )
           
       
       #step 2 training for modification
       if kwargs['weights_dir_step2'] is not None:
          weight_dir_step2 = kwargs['weights_dir_step2']
       else:
          weight_dir_step2 = os.path.join(OUTPATH,'model_basecall_best.ckpt')
       
       if kwargs['train_type']!='total_loss':
          traintype="ratioloss"
       
       if kwargs['real_candidate_files'] is None:
           loss_weights=[1., 1.,1.,0]
           traintype="modloss"
       else:
            loss_weights=[float(x) for x in kwargs['trainloss'].split("-")]
       
       print("Training for modification procedure...\n")
       print(loss_weights)
       print(weight_dir_step2)
       model.build_model_seq2seq(dim_attention=kwargs['dim_attention'], #10 vs 150
                              load_weights=True, weight_dir=weight_dir_step2, #load best weight
                              drop_input=kwargs['drop_input'],
                              drop_cnn=kwargs['drop_cnn'],
                              drop_flat=kwargs['drop_flat'],
                              BatchNorm=kwargs['BatchNorm'],
                              dim_lstm=kwargs['dim_lstm'],
                              embedding_output_dim=kwargs['embedding_output_dim'],
                              is_monotonic=kwargs['is_monotonic'],
                              Att_regularizer_weight=kwargs['Att_regularizer_weight'],
                              headnum=kwargs['headnum'],
                              W_regularizer=kwargs['W_regularizer'],
                              fc_dim = kwargs['fc_dim'],
                              fcnum = kwargs['fcnum'],
                              loss_weights=loss_weights,
                              )
       
       print("before training\n")
       model.train(samplesize,kwargs['real_candidate_files'],kwargs['real_noncandidate_files'],\
          kwargs['KO_candidate_files'],kwargs['KO_noncandidate_files'],IVT_pos_file=kwargs['IVT_pos_file'],\
          IVT_neg_file=kwargs['IVT_neg_file'],batch_size=kwargs['batch_size'], epochs=kwargs['epochs2'],\
          #x_valid=val_x,y_valid_ref=val_y_ref,y_valid_call=val_y_call,val_mod_label=val_mod_label,val_ratio_label=val_ratio_label,\
          val_real_candidate_files=kwargs['val_real_candidate_files'],val_real_noncandidate_files=kwargs['val_real_noncandidate_files'],\
          val_KO_candidate_files=kwargs['val_KO_candidate_files'],val_KO_noncandidate_files=kwargs['val_KO_noncandidate_files'],
          val_IVT_pos_file=kwargs['val_IVT_pos_file'],val_IVT_neg_file=kwargs['val_IVT_neg_file'],
          early_stop=kwargs['early_stop2'],\
          bin = kwargs['bin_size'],binsample = kwargs['bin_sample'],\
          traintype=traintype,\
          WT_label=kwargs['WT_label'],binselecttype=kwargs['binselecttype'],\
          sample_weights_flag=kwargs['sample_weights_flag'],
          steps_per_epoch=steps_per_epoch,
          shuffleIVT=kwargs['shuffleIVT'],
          weight_exp=kwargs['weight_exp'],
          KO_noncandidate_weight=kwargs['KO_noncandidate_weight'],
          real_noncandidate_weight=kwargs['real_noncandidate_weight']
          )
    
    if kwargs['predict']:
       batch_predict(kwargs['test_cachedata'],model,kwargs['test_outputfile'],kwargs['maxoutputlen'],kwargs['batch_size'],kwargs['verbose'],featuredim,kwargs['windowsize'],kwargs['printref'])
    
    K.clear_session()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    '''Model parameters'''
    parser.add_argument('--windowsize', type=int, default=11,
                        help="windowsize")
    parser.add_argument('--maxoutputlen', type=int, default=None,help="maxoutputlen")
    parser.add_argument('--dim_lstm', type=int, default=150, help='hidden size of lstm')
    parser.add_argument('--dim_attention', type=int, default=10, help='dim_attention')
    parser.add_argument('--headnum', type=int, default=10, help='number of multiheads') #select one from 3
    parser.add_argument('--dim_capsule', type=int, default=4, help='capsule dimention')
    parser.add_argument('--drop_rate', type=float, default=0.1, help='dropout ratio')
    parser.add_argument('--drop_input', type=float, default=0.2, help='dropout ratio')
    parser.add_argument('--drop_cnn', type=float, default=0.5, help='dropout ratio')
    parser.add_argument('--drop_flat', type=float, default=0.2, help='dropout ratio')
    parser.add_argument('--W_regularizer', type=float, default=0.0005, help='W_regularizer')
    parser.add_argument('--Att_regularizer_weight', type=float, default=0.0005, help='Att_regularizer_weight')
    parser.add_argument('--dataset', type=str, default='../preData/Gnid10__', help='input sequence data for training')
    parser.add_argument('--IVT_pos_file', type=str, default=None, help='directRNA features from IVT positive data')
    parser.add_argument('--IVT_neg_file', type=str, default=None, help='directRNA features from IVT negatvie data')
    parser.add_argument('--real_candidate_files', type=str, default=None, help='directRNA features (WT) on candidate A-to-I editing sites detected by bulk data')
    parser.add_argument('--real_noncandidate_files', type=str, default=None, help='directRNA features (WT) on sites that are not detected by bulk data')
    parser.add_argument('--KO_candidate_files', type=str, default=None, help='directRNA features (KO) on candidate A-to-I editing sites detected by bulk data')
    parser.add_argument('--KO_noncandidate_files', type=str, default=None, help='directRNA features (KO) on sites that are not detected by bulk data')
    #parser.add_argument('--val_hdf5', type=str, default=None, help='val_hdf5!')
    parser.add_argument('--selectrate',type=float,default=1,help='input data rate')
    parser.add_argument('--epochs1', type=int, default=100, help='')
    parser.add_argument('--epochs2', type=int, default=100, help='')
    parser.add_argument('--nb_filters', type=int, default=32, help='number of CNN filters') #select one from 3
    parser.add_argument('--filters_length', type=int, default=10, help='size of CNN filters') #select one from 3
    parser.add_argument('--filters_length1', type=int, default=10, help='size of CNN filters') #select one from 3
    parser.add_argument('--filters_length2', type=int, default=10, help='size of CNN filters') #select one from 3
    parser.add_argument('--filters_length3', type=int, default=10, help='size of CNN filters') #select one from 3
    parser.add_argument('--pooling_size', type=int, default=3, help='pooling_size') #select one from 3
    parser.add_argument('--featuredim', type=int, default=3, help='featuredim') #select one from 3
    parser.add_argument('--num_grid', type=int, default=300, help='num_grid')
    parser.add_argument('--smoothdim', type=int, default=-1, help='smoothdim if -1 no smoothing')
    parser.add_argument('--att_weight', type=float, default=1, help='number of att_weight') #select one from 3
    parser.add_argument("--BatchNorm", action="store_true",help="use BatchNorm")
    parser.add_argument("--is_monotonic", action="store_true",help="is_monotonic")
    parser.add_argument('--Resblocknum', type=int, default=3, help='Resblocknum')
    parser.add_argument('--fcnum', type=int, default=0, help='fcnum')
    parser.add_argument('--fc_dim', type=int, default=50, help='fc_dim')
    parser.add_argument('--embedding_output_dim', type=int, default=30, help='embedding_output_dim')
    parser.add_argument('--maxneg', type=int, default=30, help='maxneg')
    parser.add_argument('--predict', action="store_true", help='predict')
    parser.add_argument('--test_cachedata', type=str, default='./', help='test_cachedata!')
    parser.add_argument('--test_inputdata', type=str, default='../preData/Gnid10__', help='input sequence data for testing')
    parser.add_argument('--test_outputfile', type=str, default='../preData/Gnid10__', help='test_outputfile')
    parser.add_argument('--dataset_TN', type=str, default=None, help='input TN data for training')
    parser.add_argument('--h5py_path_TN', type=str, default=None, help='output hdf5 TN data for training')
    parser.add_argument('--get_embed', type=int, default=0, help='get_embed 1: get embed. 0: no embed')
    parser.add_argument('--bin_size',type=float,default=0.1,help='bin_size for training')
    parser.add_argument('--bin_sample',type=float,default=200,help='bin_sample size for training')
    parser.add_argument('--early_stop1',type=int,default=3, help='earlystop patience, 0 for no early_stop')
    parser.add_argument('--early_stop2',type=int,default=3, help='earlystop patience, 0 for no early_stop')
    '''Experiment parameters'''
    parser.add_argument("--message", type=str, default="", help="append to the dir name")
    parser.add_argument("--load_pretrain", action="store_true",
                        help="load pretrained CNN weights to the first convolutional layers")
    
    parser.add_argument('--batch_size',type=int,default=1000, help='batch_size < 1000 for training ')
    parser.add_argument("--weights_dir", type=str, default="",
                        help="Must specificy pretrained weights dir, if load_pretrain is set to true.") 
    
    parser.add_argument("--weights_dir_step2", type=str, default=None,help="Must specificy pretrained weights dir, if load_pretrain is set to true.") 
    
    parser.add_argument('--verbose',type=int,default=1, help='verbosity mode, 0 or 1. If 0 print nothing')
    parser.add_argument("--train_type", type=str, default="", help="train_type == total_loss?")
    parser.add_argument("--WT_label", type=int, default=1, help="WT_label:1 all 1 0: ratio")
    parser.add_argument("--binselecttype", type=int, default=1, help="binselecttype:1 uniform 2 ratio ")
    parser.add_argument("--sample_weights_flag", action="store_true", help='sample_weights_flag if added set it to true')
    parser.add_argument('--threads',type=int,default=1, help='Number of threads used for running by CPU [default: 1]')
    parser.add_argument('--device',type=str,default='GPU', help='Whether CPU or GPU is used for running ["CPU", "GPU"] [default: GPU]')
    parser.add_argument('--printref', action="store_true", help='printref')
    parser.add_argument("-steps_per_epoch", type=int, default=1000, help="steps_per_epoch")
    
    parser.add_argument('--val_IVT_pos_file', type=str, default=None, help='val_IVT_pos_file!')
    parser.add_argument('--val_IVT_neg_file', type=str, default=None, help='val_IVT_neg_file!')
    parser.add_argument('--val_real_candidate_files', type=str, default=None, help='directRNA features (WT) on candidate A-to-I editing sites detected by bulk data')
    parser.add_argument('--val_real_noncandidate_files', type=str, default=None, help='directRNA features (WT) on sites that are not detected by bulk data')
    parser.add_argument('--val_KO_candidate_files', type=str, default=None, help='directRNA features (KO) on candidate A-to-I editing sites detected by bulk data')
    parser.add_argument('--val_KO_noncandidate_files', type=str, default=None, help='directRNA features (KO) on sites that are not detected by bulk data')
    
    parser.add_argument('-shuffleIVT', action="store_true", help='shuffleIVT will be very slow for IVT data!! not recomanded to use')
    parser.add_argument('-auto_steps_per_epoch', action="store_true", help='auto_steps_per_epoch: if setted, will not use steps_per_epoch')
    parser.add_argument('-weighted_ratio_loss_flag',type=int,default=0,help="0: no ratio loss weight;1: 1-((p-1)*K.log(1-p+K.epsilon())-p*K.log(p+K.epsilon()))")
    
    parser.add_argument('-ratio_loss_weight',type=float,default=1,help="ratio_loss_weight")
    parser.add_argument('-weight_exp',type=float,default=1,help="weight_exp")
    parser.add_argument('-KO_noncandidate_weight',type=float,default=1,help="KO_noncandidate_weight")
    parser.add_argument('-real_noncandidate_weight',type=float,default=1,help="real_noncandidate_weight")
    parser.add_argument('-outputdir',type=str,default='./',help="output directory")
    parser.add_argument('-trainloss', type=str, default="1,1,1,1", help='-trainloss')
    
    args = parser.parse_args()
    
    OUTPATH = os.path.join(args.outputdir,args.message + '/')
    if not os.path.exists(OUTPATH):
        os.makedirs(OUTPATH)
    print('OUTPATH:', OUTPATH)
    
    if args.predict:
       if not os.path.exists(os.path.dirname(args.test_outputfile)):
          os.makedirs(os.path.dirname(args.test_outputfile))
    
    if args.maxoutputlen is None:
       args.maxoutputlen=args.windowsize + 2
    
    configout = open(os.path.join(OUTPATH,"config.txt"),'w')
    for k, v in vars(args).items():
        print(k, ':', v)
        configout.write(str(k)+":"+str(v)+"\n")
    
    configout.close()
    
    import datetime
    starttime = datetime.datetime.now()
    run_model(**vars(args))
    endtime = datetime.datetime.now()
    print("running time in seconds:\n")
    print((endtime - starttime).seconds)

