#! /usr/bin/python3
import basic_code_elements as basic
from scipy import ndimage
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt
import time

def main():
        img = misc.imread('Images/arch.png').astype(float)

        xSize , ySize = np.shape(img)

        tangents = basic.atd(img,9)
        np.save('tan.npy',tangents)
        print ('Tangents done')

        tangents = np.load('tan.npy')
        tx = 66
        ty = 128
        rx, ry = basic.followRidge(tangents,tx,ty,2)
        plt.imshow(img,interpolation='nearest',cmap=plt.get_cmap('gray'))

        plt.plot([tx], [ty],'ro')
        plt.plot(rx,ry,'r--')
        #plt.quiver(tangents[::4,::4,0], tangents[::4,::4,1])

        plt.show()
if __name__ == '__main__':
        main()
