#! /usr/bin/python3
import basic_code_elements as basic
from scipy import ndimage
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt
import time

def main():
        img = misc.imread('Images/Example_curve_n0.png').astype(float)

        Tangents  = basic.atd(img)
        print ('Got tangents')

        a = basic.travers(img,Tangents)
        plt.plot( [b[0] for b in a],[b[1] for b in a],'ro')
        plt.show()
if __name__ == '__main__':
        main()
