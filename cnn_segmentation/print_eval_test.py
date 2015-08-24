#!/usr/bin/env python
"""
Generate an image out of a fully convolutional neural network (CNN).
"""
import caffe
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def parse_args():
    """Parse input arguments
    """

    parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_image_file',
                        help='Input image file')
    parser.add_argument('input_deploy_file',
                        help='Input deploy file')
    parser.add_argument('input_net_file',
                        help='Input deploy file')
    parser.add_argument('output_image_file',
                        help='Output image file')


    args = parser.parse_args()
    return args

def main():
	# load image, switch to BGR, subtract mean, and make dims C x H x W for Caffe
    args = parse_args()
    im = np.array(Image.open(args.input_image_file))
    in_ = np.array(im, dtype=np.float32)
    in_ = in_[:,:,::-1]
    in_ -= np.array((104.00698793,116.66876762,122.67891434))
    in_ = in_.transpose((2,0,1))

    # load net
    net = caffe.Net(args.input_deploy_file, args.input_net_file, caffe.TEST)
    # shape for input (data blob is N x C x H x W), set data
    net.blobs['data'].reshape(1, *in_.shape)
    net.blobs['data'].data[...] = in_
    # run net and take argmax for prediction
    net.forward()
    out = net.blobs['score'].data[0].argmax(axis=0)
    # save image
    plt.imshow(out)
    plt.savefig(args.output_image_file)
    

if __name__ == '__main__':
    main()