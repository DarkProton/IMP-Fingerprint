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

        img = misc.imread('Images/arch.png').astype(float)
        xSize , ySize = np.shape(img)
        #Load in the image and get its size.

        if CAL_TANGENTS:
                tangents = basic.atd(img,3)
                np.save('tan.npy',tangents)
                print ('Tangents calculated and save')
        else:
                tangents = np.load('tan.npy')
                print('Tangents loaded')

        tx = 209
        ty = 53
        #tx and ty are the locations where the ridge finding starts


        for angle_offset in [-np.pi/2,0,np.pi/2]:
                #For all the angle offsets, calculate and plot the ridge
                rx, ry = basic.followRidge(tangents,tx,ty,1,angle_offset)
                plt.plot(rx,ry,'bo')

        plt.imshow(img,interpolation='nearest',cmap=plt.get_cmap('gray'))
        #Show the image

        plt.plot([tx], [ty],'ro')
        #Plot where the ridge finding starts

        plt.show()
        #Show image
if __name__ == '__main__':
        main()
