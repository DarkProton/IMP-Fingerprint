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

        tangents = basic.atd(img)
        placesWhereTheTangentIsNotZero = np.where( tangents[:,:,1] + tangents[:,:,0] != 0) 

        tx = 98
        ty = 152
        rx, ry = basic.followRidge(tangents,tx,ty)
        plt.imshow(img,interpolation='nearest',cmap=plt.get_cmap('gray'))

        plt.plot([tx], [ty],'ro')
        plt.plot(rx,ry,'bo')
        plt.show()
if __name__ == '__main__':
        main()
