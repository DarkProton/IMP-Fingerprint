#! /usr/bin/python3
import basic_code_elements as basic
from scipy import ndimage
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt
import time

def main():#'fingerprint5_small''Example_curve_n0', 'Example_curve_n50', 
    for n in range(4):
        files = ('arch','x')
        img = misc.imread('Images/'+files[n]+'.png').astype(float)
        xSize , ySize = np.shape(img)
        #Load in the image and get its size.

        tangents = basic.atd(img)

        plt.figure()
        plt.imshow(img[0:25,50:100], interpolation='none', \
                   cmap=plt.get_cmap('gray'))
        plt.quiver(tangents[0:25,50:100,0], tangents[0:25,50:100,1],
                   pivot='mid', color='r', units='inches', scale=5)
        #plt.show()
        plt.savefig('PND_'+files[n]+'.svg')
        print('Figure', 'PND_'+files[n]+'.svg', 'saved.')
##        
##        tx = 140
##        ty = 110
##        #tx and ty are the locations where the ridge finding starts
##
##        
##        oC= basic.followRidge(tangents, tx, ty, img, mu=1, rad=50)
##        rx = [x[0] for x in oC]
##        ry = [x[1] for x in oC]
##        plt.plot(rx,ry,'b-')
##
##        plt.imshow(img,interpolation='nearest',cmap=plt.get_cmap('gray'))
##        #Show the image
##
##        plt.plot([tx], [ty],'ro')
##        #Plot where the ridge finding starts
##
##
##        plt.show()
##        #Show image

if __name__ == '__main__':
        main()
