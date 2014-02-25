#! /usr/bin/python3
import basic_code_elements as basic
from scipy import ndimage
from scipy import misc
import numpy as np
import matplotlib.pyplot as plt

def main():
	img = misc.imread('Images/Example_curve_n0.png').astype(np.float)

	a = basic.atd(img)
	print(a)

	basic.show_pic(a[:][:][0])

if __name__ == "__main__":
	main()

