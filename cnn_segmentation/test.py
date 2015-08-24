#!/usr/bin/env python
"""
Generate an image out of a fully convolutional neural network (CNN).
"""
import caffe
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import os

def parse_args():
    """Parse input arguments
    """

    parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('input_image_folder',
                        help='Input image folder')
    parser.add_argument('input_deploy_file',
                        help='Input deploy file')
    parser.add_argument('input_net_file',
                        help='Input deploy file')
    parser.add_argument('output_image_file',
                        help='Output image file')
    parser.add_argument('output_folder_file',
                        help='Output image file')


    args = parser.parse_args()
    return args
