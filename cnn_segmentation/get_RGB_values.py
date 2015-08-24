from PIL import Image
import numpy as np

#opens txt file which contains path to camvid testing dataset
f = open('camvid_test.txt','r')
inputs = f.read().splitlines()
f.close()

N[:,:] #init N vector
n_labels= 11 #number of labels for semantic segmentation

# start for loop
for in_idx, in_ in enumerate(inputs):
im = np.array(Image.open(in_))
in_ = np.array(im, dtype=np.float32)
in_ = in_[:,:,::-1]
in_ -= np.array((104.00698793,116.66876762,122.67891434))
in_ = in_.transpose((2,0,1))

# load net
net = caffe.Net('deploy.prototxt', 'snapshot/train_iter_500.caffemodel', caffe.TEST)
# shape for input (data blob is N x C x H x W), set data
net.blobs['data'].reshape(1, *in_.shape)
net.blobs['data'].data[...] = in_
# run net and take argmax for prediction
net.forward()
out = net.blobs['score'].data[0].argmax(axis=0)
out.size=[X,Y]

if
    out[x,y]
else
    x=x+1
if x>X
    x=0
    y=y+1

