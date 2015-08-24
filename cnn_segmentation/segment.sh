#!/usr/bin/env sh
mkdir raw_images
mkdir segmented_images
ffmpeg -i Google\ Drive/A\&R\ TU\ Dortmund/4th\ Semester\ 2015/Master\ thesis/Computer\ Vision\ Thesis/CamVid\ Database/Videos/0005VD.MXF -s 960x720 raw_images/image%05d.png
./caffe-master/segmentation/semantic_segmentation.py raw_images  \
caffe-master/models/fcn_bvlc_alexnet/Alexnet_32s/deploy.prototxt  \
caffe-master/models/fcn_bvlc_alexnet/Alexnet_32s/snapshot/train_iter_500.caffemodel  \
segmented_images

