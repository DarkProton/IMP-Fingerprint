#! /usr/bin/python3
import basic_code_elements as basic
from scipy import ndimage
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt
import time

def main():
        CAL_TANGENTS = False
        #If this is true, it will calculate the tangents. Otherwise, it will load them from tan.npy

        img = misc.imread('Images/fingerprint5.png').astype(float)
        xSize , ySize = np.shape(img)
        #Load in the image and get its size.

        if CAL_TANGENTS:
                tangents = basic.pnd(img)
                np.save('tan.npy',tangents)
                print ('Tangents calculated and saved')
        else:
                tangents = np.load('tan.npy')
                print('Tangents loaded')

##        plt.figure()
##        plt.imshow(img,interpolation='none',cmap=plt.get_cmap('gray'))
##        plt.quiver(tangents[:,:,0], tangents[:,:,1],
##                   pivot='mid', color='r', units='inches', scale=1.2)
##        plt.show()
        
        tx = 850
        ty = 1190
        #tx and ty are the locations where the ridge finding starts

        
        oC= basic.followRidge(tangents, tx, ty, img, mu=1, rad=30)
        rx = [x[0] for x in oC]
        ry = [x[1] for x in oC]
        plt.plot(rx,ry,'b-')

        plt.imshow(img,interpolation='nearest',cmap=plt.get_cmap('gray'))
        #Show the image

        plt.plot([tx], [ty],'ro')
        #Plot where the ridge finding starts


        plt.show()
        #Show image
if __name__ == '__main__':
        main()
