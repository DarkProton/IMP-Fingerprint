#! /usr/bin/python3

## Basic fragments of code that will probably be useful

## Uses axes, imshow, colorbar, title, savefig, close, figure
## from matplotlib.pyplot
from matplotlib.pyplot import axes,imshow,colorbar,title,savefig,close,figure
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
            n1 = (-image[n,m] + image[n-1,m] + image[n-1,m-1] - image[n,m-1])/4
            n2 = (-image[n,m] - image[n-1,m] + image[n-1,m-1] + image[n,m-1])/4

            normals[n-1 ,  m-1] = array([n1, n2])

    return normals

## Uses zeros, sum, sqrt from numpy
def atd(image, window=9):
    from numpy import zeros, sum, sqrt
    """Returns the averaged tangent diraction of a normal array"""
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
    tangents = zeros((y_pix-1, x_pix-1))

    # Loop through and calculate the normal vector
    for n in range(0, y_pix):
        for m in range(0, x_pix):
            # Extract window from image
            window = image[n-k1:n+k2, m-k1:m+k2]

            # Calculate normalised ATD
            A = sum(window[:,:,0]**2)
            B = sum(window[:,:,1]**2)
            C = sum(window[:,:,0]*window[:,:,1])

            D = (B - A)/(2*C) - sqrt(((B - S)/(2*C))**2 + 1)
            u1 = sqrt(1/(1 + D**2))
            u2 = u1*D

            tangents[n,m] = [u1, u2]

    return tangents
