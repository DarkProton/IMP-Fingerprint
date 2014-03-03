#! /usr/bin/python3
import basic_code_elements as basic
from scipy import ndimage
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt

#fingerprint5
#Images/Example_curve_n0
img = misc.imread('arch.png').astype(float)

#a = basic.pnd(img)
b = basic.atd(img,window=3)

basic.save_pic(b[:,:,0], 'ATD_u1', colourmap='jet')
basic.save_pic(b[:,:,1], 'ATD_u2', colourmap='jet')
#basic.save_pic(a[:,:,0], 'PND_u1', colourmap='jet')
#basic.save_pic(a[:,:,1], 'PND_u2', colourmap='jet')
plt.show()
