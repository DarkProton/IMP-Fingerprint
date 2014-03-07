#! /usr/bin/python3
import basic_code_elements as basic
from scipy import ndimage
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt
import time

def main():#'fingerprint5_small', 'Example_curve_n0', 'Example_curve_n50'
    img = misc.imread('Images/arch.png').astype(float)
    xSize , ySize = np.shape(img)
    #Load in the image and get its size.

    tangents = basic.atd(img, window=5)

##    plt.figure()
##    plt.imshow(img[0:25,50:100], interpolation='none', \
##               cmap=plt.get_cmap('gray'))
##    plt.quiver(tangents[0:25,50:100,0], tangents[0:25,50:100,1],
##               pivot='mid', color='r', units='inches', scale=5)
##    plt.show()
    
    tx = 60
    ty = 5
    #tx and ty are the locations where the ridge finding starts

    plt.imshow(img[0:70,50:90], interpolation='none', cmap=plt.get_cmap('gray'))
    #Show the image
    plt.quiver(tangents[0:70,50:90,0], tangents[0:70,50:90,1],
               pivot='mid', color='r', units='inches', scale=5)
    
    oC= basic.followRidge(tangents, tx, ty, img, mu=2, rad=2)
    rx = [x[0]-50 for x in oC]
    ry = [x[1] for x in oC]
    plt.plot(rx,ry,'b-')

    plt.plot(tx-50, ty,'ro')
    #plt.show()
    #Plot where the ridge finding starts

    plt.savefig('RF_ATD_arch.svg')
    #print('Figure', 'RF_ATD_arch.svg', 'saved.')
    #Show image

if __name__ == '__main__':
        main()
