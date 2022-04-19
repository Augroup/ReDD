'''
keras version 2.2.4 
tensorflow '1.12.0'
python 3.5.2
'''

import datetime
import itertools
from collections import OrderedDict
import argparse
import os
import sys
from keras.utils.io_utils import HDF5Matrix
import h5py
#os.environ["CUDA_VISIBLE_DEVICES"] = "0"
#sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import tensorflow as tf
from keras.backend.tensorflow_backend import set_session
import keras.backend as K
from model import *  #add noncanddiate site 

#batch_size = 1000 #batch size must < 148 because smallest ratio is 148 for 0.2 bin
nb_classes = 2 #P1 P2 N

# starts training in CNN model
def run_model(**kwargs):
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
    weights_dir_list = kwargs['weights_dir'].replace("\"","").split(",")
    modellist=[]
    for weight_dir in weights_dir_list:
        print("built model by weights from "+weight_dir+"\n")
        model = multihead_attention(kwargs['windowsize'], kwargs['maxoutputlen'],featuredim,nb_classes, "./", kfold_index=0)  # initialize
        model.build_model_seq2seq(dim_attention=kwargs['dim_attention'], #10 vs 150
                              load_weights=True, weight_dir=weight_dir,
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
        
        modellist.append(model)
    
    batch_predict(kwargs['test_cachedata'],modellist,kwargs['test_outputfile'],kwargs['maxoutputlen'],kwargs['batch_size'],kwargs['verbose'],featuredim,kwargs['windowsize'],kwargs['printref'])
    
    K.clear_session()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    '''Model parameters'''
    parser.add_argument('--windowsize', type=int, default=9,
                        help="windowsize")
    parser.add_argument('--maxoutputlen', type=int, default=None,help="maxoutputlen")
    parser.add_argument('--dim_lstm', type=int, default=150, help='hidden size of lstm')
    parser.add_argument('--dim_attention', type=int, default=150, help='dim_attention')
    parser.add_argument('--headnum', type=int, default=10, help='number of multiheads') #select one from 3
    parser.add_argument('--dim_capsule', type=int, default=4, help='capsule dimention')
    parser.add_argument('--drop_rate', type=float, default=0.1, help='dropout ratio')
    parser.add_argument('--drop_input', type=float, default=0, help='dropout ratio')
    parser.add_argument('--drop_cnn', type=float, default=0.25, help='dropout ratio')
    parser.add_argument('--drop_flat', type=float, default=0.2, help='dropout ratio')
    parser.add_argument('--W_regularizer', type=float, default=0.0005, help='W_regularizer')
    parser.add_argument('--Att_regularizer_weight', type=float, default=0.0005, help='Att_regularizer_weight')
    parser.add_argument('--nb_filters', type=int, default=32, help='number of CNN filters') 
    parser.add_argument('--filters_length', type=int, default=10, help='size of CNN filters') 
    parser.add_argument('--filters_length1', type=int, default=10, help='size of CNN filters') 
    parser.add_argument('--filters_length2', type=int, default=10, help='size of CNN filters') 
    parser.add_argument('--filters_length3', type=int, default=10, help='size of CNN filters') 
    parser.add_argument('--pooling_size', type=int, default=3, help='pooling_size') 
    parser.add_argument('--featuredim', type=int, default=5, help='featuredim') 
    parser.add_argument('--num_grid', type=int, default=300, help='num_grid')
    parser.add_argument('--smoothdim', type=int, default=-1, help='smoothdim if -1 no smoothing')
    parser.add_argument('--att_weight', type=float, default=1, help='number of att_weight')
    parser.add_argument("--BatchNorm", action="store_true",help="use BatchNorm")
    parser.add_argument("--is_monotonic", action="store_true",help="is_monotonic")
    parser.add_argument('--Resblocknum', type=int, default=3, help='Resblocknum')
    parser.add_argument('--fcnum', type=int, default=0, help='fcnum')
    parser.add_argument('--fc_dim', type=int, default=50, help='fc_dim')
    parser.add_argument('--embedding_output_dim', type=int, default=30, help='embedding_output_dim')
    parser.add_argument('--test_cachedata', type=str, help='test_cachedata')
    parser.add_argument('--test_outputfile', type=str, help='test_outputfile')
    parser.add_argument('--get_embed', type=int, default=0, help='get_embed 1: get embed. 0: no embed')
    parser.add_argument('--batch_size',type=int,default=10000, help='batch_size in prediction')
    parser.add_argument("--weights_dir", type=str, default="",
                        help="Must specificy pretrained weights dir, if load_pretrain is set to true.") 
    
    parser.add_argument('--verbose',type=int,default=1, help='verbosity mode, 0 or 1. If 0 print nothing')
    parser.add_argument("--train_type", type=str, default="", help="train_type == total_loss?")
    parser.add_argument("--WT_label", type=int, default=1, help="WT_label:1 all 1 0: ratio")
    parser.add_argument('--threads',type=int,default=1, help='Number of threads used for running by CPU [default: 1]')
    parser.add_argument('--device',type=str,default='GPU', help='Whether CPU or GPU is used for running ["CPU", "GPU"] [default: GPU]')
    parser.add_argument('--printref', action="store_true", help='printref')
    
    args = parser.parse_args()
    
    if not os.path.exists(os.path.dirname(args.test_outputfile)):
          os.makedirs(os.path.dirname(args.test_outputfile))
    
    
    if args.maxoutputlen is None:
       args.maxoutputlen=args.windowsize+2 
    
    import datetime
    starttime = datetime.datetime.now()
    run_model(**vars(args))
    endtime = datetime.datetime.now()
    print("running time in seconds:\n")
    print((endtime - starttime).seconds)




