#! /usr/bin/python3
import basic_code_elements as basic
from scipy import ndimage
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt

img = misc.imread('Images/Example_curve_n0.png').astype(np.float)

a = basic.atd(img)

basic.save_pic(a[:,:,0], 'test')

