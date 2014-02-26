#! /usr/bin/python3
import basic_code_elements as basic
from scipy import ndimage
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt

img = misc.imread('Images/Example_curve_n0.png').astype(float)
a = basic.pnd(img)


#a = basic.atd(img)

#basic.save_pic(a[:,:,0], 'ATD_u1', colourmap='jet')
#basic.save_pic(a[:,:,1], 'ATD_u2', colourmap='jet')
#plt.show()
