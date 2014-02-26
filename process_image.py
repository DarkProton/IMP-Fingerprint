#! /usr/bin/python3
import basic_code_elements as basic
from scipy import ndimage
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt
import time

def main():
        img = misc.imread('Images/Example_curve_n0.png').astype(float)

        directionTangents  = basic.atd(img)
        a = basic.travers(img,basic.atd(img))
        print (a)
if __name__ == '__main__':
        main()
