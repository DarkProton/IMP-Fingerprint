#! /usr/bin/python3
import basic_code_elements as basic
from scipy import ndimage
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt
import time

def main():
        img = misc.imread('Images/Example_curve_n0.png').astype(float)
        basic.travers(img,2)
        return

        a = basic.atd(img)

        plt.quiver(a[:,:,0],a[:,:,1])
        plt.show()

if __name__ == '__main__':
        main()
