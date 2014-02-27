#! /usr/bin/python3

## Basic fragments of code that will probably be useful

## Uses axes, imshow, colorbar, title, savefig, close, figure
## from matplotlib.pyplot
from matplotlib.pyplot import axes, imshow, colorbar, title, \
     savefig, close, figure
def show_pic(image,plot_name='Test Image',colourmap=None):
    """Method to simply shwo the image. Useful for tests"""
    figure(facecolor='white', figsize=(5,4))

    # Remove axes from the plot
    ax = axes(frameon=False)
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)

    # Give the image a title
    title(plot_name)
    
    if colourmap != None:
        # Make image with a colourmap
        imshow(image, interpolation='none',cmap=colourmap)
        colorbar()
    else:
        # Make an image without a colourmap
        imshow(image, interpolation='none')


def save_pic(image, plot_name, colourmap=None):
    """A function that nicley abstracts producing an image object"""
    figure(facecolor='white', figsize=(5,4))

    # Remove axes from the plot
    ax = axes(frameon=False)
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)

    # Give the image a title
    title(plot_name)
    
    if colourmap != None:
        # Make image with a colourmap
        imshow(image, interpolation='none',cmap=colourmap)
        colorbar()
    else:
        # Make an image without a colourmap
        imshow(image, interpolation='none')
    savefig("{:s}.pdf".format(plot_name).replace(' ', '_'))
    close()

## Uses sqrt, sum
def normalise(array):
    """Normalises a 3D array"""
    from numpy import sqrt
    yl = len(array)
    xl = len(array[0])
    for n in range(yl):
        for m in range(xl):
            y = array[n,m,0]
            x = array[n,m,1]

            if x == 0 and y == 0:
                array[n,m,:] = 0 , 0
            else:
                array[n,m,:] = array[n,m,:]*(x**2 + y**2)**(-0.5)

    return array

## Uses zeros and array from numpy
def pnd(image):
    """Returns the point normal determination"""
    from numpy import zeros,array
    # Get image dimensions
    y_pix = len(image)
    x_pix = len(image[0])

    # Make array to hold normals
    normals = zeros((y_pix-1, x_pix-1,2))

    # Loop through and calculate the normal vector
    for n in range(1, y_pix):
        for m in range(1, x_pix):
            n1 = (-image[n, m] + image[n-1, m] + image[n-1, m-1] \
                  - image[n, m-1])/4
            n2 = (-image[n, m] - image[n-1, m] + image[n-1, m-1] \
                  + image[n, m-1])/4

            normals[n-1, m-1, :] = n1, n2

    normals = normalise(normals)

    return normals

## Uses zeros, sum, sqrt from numpy
def atd(image, window=9):
    """Returns the averaged tangent diraction of a normal array"""
    from numpy import zeros, sum, sqrt, array
    # Get image dimensions
    y_pix = len(image)
    x_pix = len(image[0])

    # Calculate window half-size
    k1 = int(window/2)
    if k1 == window/2:
        k2 = int(window/2) - 1
    else:
        k2 = k1

    # Make array to hold tangents
    tangents = zeros((y_pix-1, x_pix-1,2))

    # Calcualte normals
    normals = pnd(image)

    # Test to see if a pixel is in the image
    inImageTest = lambda ky,kx : not ((ky < 0) and (kx < 0) and \
                                      (ky > y_pix - k1 - 1) and \
                                      (kx > x_pix - k1 - 1))

    # Loop through and calculate ATD
    for n in range(1,y_pix-(k1+1)):
        for m in range(1,x_pix-(k1+1)):
            # Extract window from image
            window = [normals[ky,kx] for ky in range(n-k1,n+k2) \
                      for kx in range(m-k1,m+k2) if inImageTest(ky,kx)]

            # Make window a numpy array for easier indexing
            window = array(window)

            # Calculate normalised ATD
            A = sum(window[:,0]**2)
            B = sum(window[:,1]**2)
            C = sum(window[:,0]*window[:,1])

            # Diagonal matrix case
            if C == 0 and A < B:
                u1 = 1
                u2 = 0

            elif C == 0 and A > B:
                u1 = 0
                u2 = 1

            # Plateau
            elif A == B:
                u1 = 0
                u2 = 0

            # Non-diagonal case
            else:
                D = (B - A)/(2*C) - sqrt(1 + ((B - A)/(2*C))**2) 
                u1 = 1/sqrt(1 + D**2)
                u2 = u1*D

            tangents[n,m,:] = u1, u2

    return tangents

def OC1(A, B, C, D, E, xl, yl):
    """Calculates the curvature point via method 1 - needs ADT results"""
    # Straight line case
    if A*B == C**2:
        x = xl/2
        y = yl/2

    # Curved line case
    else:
        x = (B*D - C*E)/(A*B - C**2)
        y = (A*E  - C*D)/(A*B - C**2)

    pc = x, y
    
    return pc

def OC2(A, B, C, D, E, M, xl, yl):
    """Calculates  the curvature point via method 2 - needs ADT results"""
    # Two points case
    if not M == 0:
        a = D**2/E + (A - B)*D - E
        b = (B - A)*M - D**2 - E**2 + 2*M*D/E
        c = (C*M/E + D)*M

        x = (-b - (b**2 - 4*a*c))/(2*a)
        y = (M - D*x)/E

    # Infinite curvature case
    else:
        x = xl/2
        y = yl/2

    pc = x, y

    return pc
        
## Uses sum from numpy
def OC_switch(window, threshold=0.1):
    """Decided which OC function to use"""
    from numpy import sum
    # Get the window size
    yl = len(window)
    xl = len(window)

    # Calculate vector weighting parameter
    r = [n*window[n,m,0] + m*window[n,m,1] for n in range(yl) \
         for m in range(xl)]

    # Calculate curvature paramaters
    A = sum(window[:,:,0]**2)
    B = sum(window[:,:,1]**2)
    C = sum(window[:,:,0]*window[:,:,1])
    D = sum(window[:,:,0]*r)
    E = sum(window[:,:,1]*r)
    M = sum(r**2)

    # Calculate decision eigenvalue
    l = (A + B - ((A - B)**2 + 4*C**2)**0.5)/(2*yl*xl)

    # Decide which OC to use
    if l > threshold:
        pc = OC1(A, B, C, D, E, xl, yl)
    else:
        pc = OC2(A, B, C, D, E, M, xl, yl)

    return pc
