#! /usr/bin/python3

## Basic fragments of code that will probably be useful

## Uses axes, imshow, colorbar, title, savefig, close, figure
## from matplotlib.pyplot
def show_pic(image, plot_name, colourmap=None):
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

## Uses zeros from numpy
def pnd(image):
    """Returns the point normal determination"""
    y_pix = len(image)
    x_pix = len(image[0])
    normals = zeros((y_pix-1, x_pix-1))

    for n in range(1, y_pix):
        for m in range(1, x_pix):
            n1 = (-image[n,m] + image[n-1,m] + image[n-1,m-1] - image[n,m-1])/4
            n2 = (-image[n,m] - image[n-1,m] + image[n-1,m-1] + image[n,m-1])/4

            normals[n-1, m-1] = [n1, n2]

    return normals
