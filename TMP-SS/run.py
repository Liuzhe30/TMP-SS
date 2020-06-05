#! /usr/bin/env python
# -*- coding:utf-8 -*-

import tensorflow as tf
from keras import layers, models, optimizers
from keras.utils import to_categorical
from keras.layers import *
from keras.models import *
from keras.callbacks import Callback
from keras.initializers import Constant
from keras import backend as K
K.clear_session()
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from keras import callbacks
from keras import backend as K 
from keras.utils.np_utils import to_categorical
import time
K.set_image_data_format('channels_last')
np.random.seed(0)

# topology structure map: 0 for Non-transmembrane region, 1 for Transmembrane region
topo_dict = {0: 'N', 1: 'T'}

# secondary structure map: 0 for Helix, 1 for Strand, 2 for Coil
ss_dict = {0: 'H', 1: 'E', 2: 'C'}

# one-hot map
dict = {'C':0, 'D':1, 'S':2, 'Q':3, 'K':4,
        'I':5, 'P':6, 'T':7, 'F':8, 'N':9,
        'G':10, 'H':11, 'L':12, 'R':13, 'W':14,
        'A':15, 'V':16, 'E':17, 'Y':18, 'M':19}

# metric
def cc(y_true, y_pred):
    x = y_true
    y = y_pred
    mx = K.mean(x)
    my = K.mean(y)
    xm, ym = x-mx, y-my
    r_num = K.sum(tf.multiply(xm,ym))
    r_den = K.sqrt(tf.multiply(K.sum(K.square(xm)), K.sum(K.square(ym))))
    r = r_num / r_den

    r = K.maximum(K.minimum(r, 1.0), -1.0)
    return K.square(r)

# attention mechanism
def attention_3d_block(inputs):
    a = Permute((2, 1))(inputs)
    a = Dense(nb_time_steps, activation='relu', name = 'attention_dense')(a)
    a_probs = Permute((2, 1), name='attention_vec')(a)
    #output_attention_mul = merge([inputs, a_probs], name='attention_mul', mode='mul')
    output_attention_mul = multiply([inputs, a_probs], name='attention_mul')
    return output_attention_mul

#parameters for LSTM
window_length = 19
nb_lstm_outputs = 700  
rows, cols = window_length, 52
nb_time_steps = window_length

if __name__ == "__main__":
        
    '''
    cmd = python run.py --fasta sample/sample.fasta --hhblits_path sample/hhblits/ --output_path results/
    cmd = python run.py -f sample/sample.fasta -p sample/hhblits/ -o results/
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', default='sample/sample.fasta')
    parser.add_argument('-p', '--hhblits_path', default='sample/hhblits/')
    parser.add_argument('-o', '--output_path', default='results/')
    args = parser.parse_args()
    #print(args)    
    
    # calculate running time
    time_start=time.time()
    
    # generate data
    from pre_processing import Processor
    processor = Processor()
    fasta = args.fasta
    hhblits_path = args.hhblits_path   
    output_path = args.output_path   
    x_test = processor.data_pre_processing(fasta, hhblits_path, window_length) 
    
    # build model
    input1 = Input(shape=(window_length, 21),name = 'input1')
    re = Reshape((window_length, 21, 1))(input1)
    conv1 = Conv2D(filters=16, kernel_size=3, strides=1, padding='same', activation='relu')(re)  
    conv1 = BatchNormalization(axis=-1)(conv1)
    conv1 = Conv2D(filters=32, kernel_size=5, strides=1, padding='same', activation='relu')(conv1)  
    conv1 = MaxPooling2D((3, 3), strides=(1, 1), padding='same')(conv1)
       
    #print('out1.get_shape()', out1.get_shape())
    out1 = Reshape((window_length, -1))(conv1)
    
    input2 = Input(shape=(window_length, 31),name = 'input2') 
    re = Reshape((window_length, 31, 1))(input2)
    conv1 = Conv2D(filters=15, kernel_size=3, strides=1, padding='same', activation='relu')(re)
    conv1 = BatchNormalization(axis=-1)(conv1)
    conv1 = Conv2D(filters=32, kernel_size=5, strides=1, padding='same', activation='relu')(conv1)  
    conv1 = MaxPooling2D((3, 3), strides=(1, 1), padding='same')(conv1)  
    conv1 = Conv2D(filters=64, kernel_size=7, strides=1, padding='same', activation='relu')(conv1)  
    conv1 = MaxPooling2D((3, 3), strides=(1, 1), padding='same')(conv1)  
    out2 = Reshape((window_length, -1))(conv1)
    #print('lstm_out2.get_shape()', lstm_out2.get_shape())  
    
    merged = concatenate([out1,out2], axis=-1)
    lstm_out = Bidirectional(LSTM(nb_lstm_outputs, return_sequences=True), name='bilstm1')(merged)
    lstm_out2 = Bidirectional(LSTM(nb_lstm_outputs, return_sequences=True), name='bilstm2')(lstm_out)
    attention_mul = attention_3d_block(lstm_out2)
    attention_flatten = Flatten()(attention_mul)
    drop2 = Dropout(0.3)(attention_flatten)
    fc1 = Dense(700, activation='relu',kernel_initializer='random_uniform',
                bias_initializer='zeros')(drop2)       
    fc2 = Dense(100, activation='relu')(fc1)
    output1 = Dense(3, activation='softmax', name = 'output_1')(fc2)
    output2 = Dense(2, activation='softmax', name = 'output_2')(fc2)
    
    model = Model(inputs=[input1, input2], outputs=[output1, output2])  
    model.summary()
    
    # load weights
    model.load_weights('models/ss_topo_weights.h5')
    
    # predict
    [y_pred1,y_pred2] = model.predict([x_test[:,:,0:21], x_test[:,:,21:52]], batch_size=128)
    time_end=time.time()
    
    # normalize results
    ss_pred = np.argmax(y_pred1, axis=-1)  
    topo_pred = np.argmax(y_pred2, axis=-1)
    print('finished!')
    
    with open(output_path + '/result.txt','w+') as w:
        w.write('Thank you for using TMPSS!\n')
        w.write('Running time: ' + str(time_end-time_start) + ' s\n')
        w.write('Note:\n')
        w.write('1. secondary structure map: H for Helix, E for Strand, C for Coil\n')
        w.write('2. topology structure map: N for Non-transmembrane region, T for Transmembrane region\n')
        with open(fasta) as get_fasta:
            temp = get_fasta.readline()
            pdb_id = ""
            index = 0
            while temp:
                if(temp[0] == ">"):
                    pdb_id = temp[1:].strip()
                    temp = get_fasta.readline()
                    w.write('\nPDB_ID: ' + str(pdb_id) + '\n')
                    w.write('Residue\tSS\tTOPO\t\n')
                    continue
                for i in temp:
                    if(i != '\n'):
                        ss = ss_dict[ss_pred[index]]
                        topo = topo_dict[topo_pred[index]]
                        index += 1
                        w.write(i + '\t' + ss + '\t' + topo + '\t\n')
                temp = get_fasta.readline()        
        
     
