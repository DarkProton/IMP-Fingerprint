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
            n1 = (-image[n,m] + image[n-1,m] + image[n-1,m-1] \
                  - image[n,m-1])/4
            n2 = (-image[n,m] - image[n-1,m] + image[n-1,m-1] \
                  + image[n,m-1])/4

            normals[n-1,  m-1] = array([n1, n2])

    return normals

## Uses zeros, sum, sqrt from numpy
def atd(image, window=9):
    from numpy import zeros, sum, sqrt, array
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
    tangents = zeros((y_pix-1, x_pix-1,2))

    normals = pnd(image)

    for n in range(1,y_pix-(k1+1)):
            print(n)
            for m in range(1,x_pix-(k1+1)):
                    A = 0
                    B = 0
                    C = 0
                    for ky in range(-3,4):
                            for kx in range(-3,4):
                                    if (n + ky) < 0 or (n + ky) > y_pix-5 or (m + kx) < 0 or (m + kx) > x_pix-5:
                                            continue
                                            #This will probably need to be changed.
                                    v = normals[n + ky][m + kx]
                                    A += v[0] ** 2
                                    B += v[1] ** 2
                                    C += v[0] * v[1]
                    if C == 0:
                            if A < B:
                                    u1 = 1
                                    u2 = 0
                            else:
                                    u1 = 0
                                    u2 = 1
                    else:
                        D = (B - A)/(2*C) - sqrt( 1 + ( (B - A) / (2 * C) ) ** 2) 

                        u1 = (1 + D**2) ** (-0.5)
                        u2 = u1 * D
                    tangents[n,m,:] = array([u1,u2])

    return tangents
