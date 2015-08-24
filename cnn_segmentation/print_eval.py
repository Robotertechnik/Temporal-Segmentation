#!/usr/bin/env python
"""
Generate an image out of a fully convolutional neural network (CNN).
"""
import caffe
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
# from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

# load image, switch to BGR, subtract mean, and make dims C x H x W for Caffe
im = np.array(Image.open('caffe-master/data/camvid/701_StillsRaw_full/0006R0_f01590.png'))
#im = np.array(Image.open('../../../data/camvid/701_StillsRaw_full/0006R0_f01590.png').resize((240, 180), Image.ANTIALIAS))
in_ = np.array(im, dtype=np.float32)
in_ = in_[:,:,::-1]
in_ -= np.array((104.00698793,116.66876762,122.67891434))
in_ = in_.transpose((2,0,1))

# load net
net = caffe.Net('caffe-master/models/fcn_bvlc_alexnet/Alexnet_32s/deploy.prototxt', 'caffe-master/models/fcn_bvlc_alexnet/Alexnet_32s/snapshot/train_iter_500.caffemodel', caffe.TEST)
# shape for input (data blob is N x C x H x W), set data
net.blobs['data'].reshape(1, *in_.shape)
net.blobs['data'].data[...] = in_
# run net and take argmax for prediction
net.forward()
out = net.blobs['score'].data[0].argmax(axis=0)

#print out
plt.imshow(out)
plt.show()