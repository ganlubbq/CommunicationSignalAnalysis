# -*- coding: utf-8 -*-
"""
Created on Thu May  9 23:21:42 2019

@author: glingi
"""

# Source from
# https://github.com/radioML/examples/blob/master/modulation_recognition/RML2016.10a_VTCNN2_example.ipynb

# Import all the things we need ---
#   by setting env variables before Keras import 
#  tensorflow backend and which GPU it uses
#%matplotlib inline
import random
import numpy as np
#import theano as th
#import theano.tensor as T
from keras.utils import np_utils
from keras.layers.core import Reshape,Dense,Dropout,Activation,Flatten
from keras.layers.noise import GaussianNoise
from keras.layers.convolutional import Convolution2D, MaxPooling2D, ZeroPadding2D
from keras.regularizers import *
from keras.optimizers import adam
import matplotlib.pyplot as plt
import seaborn as sns
#import h5py
import pickle
#import  keras1

# Load the dataset ...
#  You will need to seperately download or generate this file


file_path = 'D:\\python_script\\datasets\\RML2016.10a\\RML2016.10a_dict.pkl'
with open(file_path, 'rb') as file_pointer:
#    u = pickle._Unpickler(file_pointer)
#    u.encoding = 'latin1'
#    p = u.load()
    Xd = pickle.load(file_pointer, encoding='latin1')

"""

# for hdf5 format
Data are stored in hdf5 format as complex floating point values, 
with 2 million examples, each 1024 samples long.

"""
#file_path = 'D:\\python_script\\datasets\\2018.01.OSC.0001_1024x2M.h5\\2018.01\\'
#Xd = h5py.File(file_path + "GOLD_XYZ_OSC.0001_1024.hdf5",'r')
#
## Display all groups
#print('Keys: {}'.format(Xd.keys()))
#group_key = list(Xd.keys())[0]
#
## Get the data
#num_samples = 10
#sample_length = 1024
#data=list(Xd[group_key][0:num_samples])
#
#plt.plot(data[8])

# pair of modulation scheme and signal-to-noise ratio
group_keys = sorted(list(Xd.keys()))

signal_iq = list()
labels_modulation = list()
labels_snr = list()
for modulation, snr in group_keys:
    signal_iq.append(Xd[(modulation, snr)])
    for i in range(Xd[(modulation, snr)].shape[0]):
        labels_modulation.append(modulation)
        labels_snr.append(snr)

# make one column 
X = np.vstack(signal_iq)       

# Partition the data
#  into training and test sets of the form we can train/test on 
#  while keeping SNR and Mod labels handy for each
np.random.seed(2016)
num_examples = X.shape[0] # 220,000
num_train = int(num_examples * 0.5)    # prepare 50% of data for training set 
# prepare indices for data
data_idx = np.random.choice(range(0, num_examples), size=num_train, replace=False)
train_idx = data_idx[:int(num_train / 2)]
test_idx = data_idx[int(num_train / 2) + 1:]

# training data, test data
train_data = X[train_idx]
test_data = X[test_idx]

# what is the function of map function?
classes = list(set(labels_modulation))
one_hot_labels = np.zeros(shape=(len(labels_modulation), len(classes)))
for idx, modulation in enumerate(labels_modulation):
    m_idx = classes.index(modulation)
    one_hot_labels[idx, m_idx] = 1
    
# training label, test label
train_labels = one_hot_labels[train_idx, :]
test_labels = one_hot_labels[test_idx, :]

input_shape = list(train_data.shape[1:])
print(train_data.shape, input_shape)


# Build VT-CNN2 Neural Net model using Keras primitives -- 
#  - Reshape [N,2,128] to [N,1,2,128] on input
#  - Pass through 2 2DConv/ReLu layers
#  - Pass through 2 Dense layers (ReLu and Softmax)
#  - Perform categorical cross entropy optimization

import keras
from keras import layers
from keras import models

#Model construction
dr = 0.5 # dropout rate (%)
model = models.Sequential()
model.add(Reshape(input_shape+[1], input_shape=input_shape))
model.add(layers.ZeroPadding2D(padding=(0, 2))) # height, width
model.add(layers.Conv2D(256, (1, 3), 
                        #border_mode='valid',
                        activation='relu',
                        name='conv1'))
                        #init='glorot_uniform'))
model.add(layers.Dropout(dr))
#model.add(layers.MaxPooling2D(1, 1))
model.add(layers.ZeroPadding2D(padding=(0, 2)))
model.add(layers.Conv2D(80, (2, 3), 
                        #border_mode='valid',
                        activation='relu',
                        name='conv2'))
                        #init='glorot_uniform'))
model.add(layers.Dropout(dr))
#model.add(layers.MaxPooling2D(1, 2))
model.add(layers.Flatten())
model.add(layers.Dense(256,
                       activation='relu',
                       #init='he_normal',
                       name="dense1"))
model.add(layers.Dropout(dr))
#model.add(layers.MaxPooling2D(1, 2))
model.add(layers.Dense(len(classes),
                       #init='he_normal',
                       activation='softmax',
                       name="dense2" ))
model.add(Reshape([len(classes)]))
# 모델 컴파일
model.compile(loss='categorical_crossentropy',
              #optimizer='adam',
              optimizer='rmsprop',
              metrics=['accuracy'])
model.summary()

# Set up some params 
nb_epoch = 100     # number of epochs to train on
batch_size = 1024  # training batch size

# perform training ...
#   - call the main training loop in keras for our network+dataset
weight_file_path = 'D:\\python_script\\datasets\\'
weight_file_path = weight_file_path + 'conv_modulation_recognition_nets_CNN2_0.5.wts.h5'

train_data = X[train_idx]
test_data = X[test_idx]
history = model.fit(train_data,
                    train_labels,
                    batch_size = batch_size,
                    epochs = nb_epoch,
                    verbose = 2,
                    validation_data=(test_data, test_labels),
                    callbacks = [
                        keras.callbacks.ModelCheckpoint(
                                weight_file_path,
                                monitor='val_loss', # val_loss가 좋아지지 않으면 모델파일을 덮어쓰지 않는다.
                                verbose=0,
                                save_best_only=True, # 훈련하는 동아 가장 좋은 모델이 저장됨
                                mode='auto'),
                        keras.callbacks.EarlyStopping(
                                monitor='val_loss', # 모델의 검증 소실을 모니터링
                                patience=5, 
                                verbose=0, 
                                mode='auto')
        ])
        
# we re-load the best weights once training is finished
model.load_weights(weight_file_path)

model.compile(loss='categorical_crossentropy',
              #optimizer='adam',
              optimizer='rmsprop',
              metrics=['accuracy'])
# Evaludate and Plot Model Performance
# Show simple version and performance
import matplotlib.pyplot as plt

plt.figure(1);
loss = history.history['loss']
val_loss = history.history['val_loss']

epochs = range(1, len(loss) + 1)

plt.plot(epochs, loss, 'bo', label='Training loss')
plt.plot(epochs, val_loss, 'b', label='Validation loss')
plt.title('Training and Validation loss')
plt.xlabel('Epochs')
plt.ylabel('loss')
plt.legend()

plt.figure(2);
loss = history.history['acc']
val_loss = history.history['val_acc']

epochs = range(1, len(loss) + 1)

plt.plot(epochs, loss, 'bo', label='Training acc')
plt.plot(epochs, val_loss, 'b', label='Validation acc')
plt.title('Training and Validation accuracy')
plt.xlabel('Epochs')
plt.ylabel('Accuracy')
plt.legend()

plt.show()

def plot_confusion_matrix(cm, title='Confusion matrix', cmap=plt.cm.Blues, labels=[]):
    
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title)
    plt.colorbar()
    tick_marks = np.arange(len(labels))
    plt.xticks(tick_marks, labels, rotation=45)
    plt.yticks(tick_marks, labels)
    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    
# Plot confusion matrix
test_labels_hat = model.predict(test_data, batch_size=batch_size)
conf = np.zeros([len(classes), len(classes)])
confnorm = np.zeros([len(classes), len(classes)])
for idx in range(0, test_data.shape[0]):
    j = list(test_labels[idx, :]).index(1)
    k = int(np.argmax(test_labels_hat[idx, :]))
    conf[j, k] = conf[j, k] + 1
for idx in range(0, len(classes)):
    confnorm[idx, :] = conf[idx, :] / np.sum(conf[idx, :])

plot_confusion_matrix(confnorm, labels=classes)

# Plot confusion matrix
acc = {}
test_SNRs = np.array(np.array(labels_snr)[test_idx])
types_of_test_SNRs = np.array(sorted(set(np.array(labels_snr)[test_idx])))
for snr in types_of_test_SNRs:
    # extract classes at SNR    
    test_data_i = test_data[np.where(test_SNRs == snr)]
    test_labels_i = test_labels[np.where(test_SNRs == snr)]
    
    # estimate classes
    test_labels_i_hat = model.predict(test_data_i)
    conf = np.zeros([len(classes), len(classes)])
    confnorm = np.zeros([len(classes), len(classes)])
    for i in range(0, test_data_i.shape[0]):
        j = list(test_labels_i[i, :]).index(1)
        k = int(np.argmax(test_labels_i_hat[i, :]))
        conf[j, k] = conf[j, k] + 1
        
    for i in range(0, len(classes)):
        confnorm[i, :] = conf[i, :] / np.sum(conf[i, :])
    plt.figure()
    plot_confusion_matrix(confnorm, labels=classes, 
                          title='ConvNet Confusion Matrix (SNR={})'.format(snr))
    correct = np.sum(np.diag(conf))
    uncorrect = np.sum(conf) - correct
    print('Overall Accuracy {}'.format(correct / (correct + uncorrect)))
    acc[snr] = 1.0 * correct / (correct + uncorrect)

# Save results to a pickle file for plotting later
print(acc)
fd = open('results_cnn2_d0.5.dat','wb')
pickle.dump( ("CNN2", 0.5, acc) , fd )
plt.figure()
# Plot accuracy curve
plt.plot(types_of_test_SNRs, list(acc))
plt.xlabel("Signal to Noise Ratio")
plt.ylabel("Classification Accuracy")
plt.title("CNN2 Classification Accuracy on RadioML 2016.10 Alpha")
