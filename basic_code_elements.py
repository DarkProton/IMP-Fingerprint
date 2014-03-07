#! /usr/bin/python3

## Basic functions for the construction of a ridge
## following system

from numpy import sqrt, zeros, array, sum, sqrt, mean, \
     sign, shape, where, amax, angle, argmin
from math import cos, sin, pi

def normalise(array):
    """Normalises a 2D array of 2D vectors"""
    
    # Get array dimensions
    yl = len(array)
    xl = len(array[0])
    
    for n in range(yl):
        for m in range(xl):
            # Get vector components
            y = array[n,m,0]
            x = array[n,m,1]

            # check is the vector is zero
            if x == 0 and y == 0:
                array[n,m,:] = 0, 0
            else:
                array[n,m,:] = array[n,m,:]* \
                               (x**2 + y**2)**(-0.5)

    return array

def pnd(image):
    """Returns the point normal determination"""

    # Get image dimensions
    y_pix = len(image)
    x_pix = len(image[0])

    # Make array to hold normals
    normals = zeros((y_pix-1, x_pix-1,2))

    for n in range(1, y_pix):
        for m in range(1, x_pix):
            # Fit plane and calculate normal
            n1 = (-image[n-1, m] + image[n-1, m-1] + \
                  image[n, m-1] - image[n, m])/4
            n2 = (-image[n-1, m] - image[n-1, m-1] + \
                  image[n, m-1] + image[n, m])/4

            normals[n-1, m-1, :] = n1, n2

    # normalise the normals
    normals = normalise(normals)

    return normals

def atd(image, window=3):
    """Averages the tangent diraction of a PND array"""
    
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
    inImageTest = lambda ky,kx: not \
                  ((ky < 0) and (kx < 0) and \
                  (ky > y_pix - k1 - 1) and \
                  (kx > x_pix - k1 - 1))

    # Loop through and calculate ATD
    for n in range(1,y_pix-(k1+1)):
        for m in range(1,x_pix-(k1+1)):
            # Extract window from image
            window = [normals[ky,kx] for ky in \
                      range(n-k1,n+k2) for kx in \
                      range(m-k1,m+k2) \
                      if inImageTest(ky,kx)]

            # Resize for easier indexing
            window = array(window)

            # Calculate normalised ATD
            A = sum(window[:,0]**2)
            B = sum(window[:,1]**2)
            C = sum(window[:,0]*window[:,1])
            mn1 = mean(window[:,0])

            # Diagonal matrix cases
            if C == 0 and A < B:
                u1 = 1
                u2 = 0

            elif C == 0 and A > B:
                u1 = 0
                u2 = 1

            # Plateau vase
            elif ((A+B) + sqrt((A-B)**2 + 4*C**2)) == \
                  ((A+B) - sqrt((A-B)**2 + 4*C**2)):
                u1 = 0
                u2 = 0

            # Non-diagonal case
            else:
                D = (B - A)/(2*C) - \
                    sqrt(1 + ((B - A)/(2*C))**2)
                # Ensure that u1 orents properly
                if mn1 > 0:
                    u1 = 1/sqrt(1 + D**2)
                else:
                    u1 = -1/sqrt(1 + D**2)
                    
                u2 = u1*D

            tangents[n,m,:] = u1, u2

    return tangents

def followRidge(tangents, cX, cY, img, mu=5, rad=3):
    """
       Given an list of tangents and an input starting
       position, returns the list of points on a ridge
    """

    # Get size of the image.
    ySize, xSize, zSize  = shape(tangents)

    # Create array to hold traversed pixel
    visited = array([[False]*ySize]*xSize)
    double_visied = array([[False]*ySize]*xSize)

    # Calculate angles
    angels = angle(tangents[:,:,1] + 1j*tangents[:,:,0])

    # Calculate line length
    ll = rad*2+1

    # Array for holding line greylevels, angles and
    # coordinates
    ang_line = zeros((ll,2))
    grey_line = zeros(ll)
    coord_line = zeros((ll,2)).astype('int')

    # Set starting angle
    psi_s = pi/4

    # Record starting position
    uC = [[cX,cY]]

    # Only visit pixels that haven't been visited before
    while cX - rad - mu >= 0 and cY - rad - mu >= 0 \
          and cX + rad + mu < xSize \
          and cY + rad + mu < ySize \
          and not visited[cX,cY]:

            # Mark current pixel as visited
            visited[cX,cY] = True

            # Extract line for processing
            for sig in range(-rad, rad+1):
                # Calculate coords of pixels in line
                lx = cX + round(sig*cos(psi_s + pi/2))
                ly = cY + round(sig*sin(psi_s + pi/2))

                # Extract greylevel of line
                grey_line[sig+rad] = img[ly, lx]

                # Store coordinates of pixels in line
                coord_line[sig] = ly, lx

            # Move to ridge peak in line
            cY, cX = coord_line[argmin(grey_line)]

            # Reset line angle counters and angles
            ang_c1 = 0
            ang1 = 0
            ang_c2 = 0
            ang2 = 0

            # Extract line for processing at new coords
            for sig in range(-rad, rad+1):
                # Calculate coords of pixels in line
                lx = cX + round(sig*cos(psi_s + pi/2))
                ly = cY + round(sig*sin(psi_s + pi/2))

                # Extract angle
                if sig < 0 \
                   and tangents[ly,lx,0] != 0 \
                   and tangents[ly,lx,1] != 0:
                    ang1 += angels[ly, lx]
                    ang_c1 += 1

                elif sig > 0 \
                     and tangents[ly,lx,0] != 0 \
                     and tangents[ly,lx,1] != 0:
                    ang2 += angels[ly, lx]
                    ang_c2 += 1

                # Store coordinates of pixels in line
                coord_line[sig] = ly, lx
                
            # Calculate greyscale tangent
            if ang_c1:
                ang1 = ang1/ang_c1
            else:
                ang1 = psi_s

            if ang_c2:
                ang2 = ang2/ang_c2
            else:
                ang2 = psi_s

            # Add/substract 90deg to minimise distance
            # from phi_c
            if abs(psi_s - ang1 - pi/2) < \
               abs(psi_s - ang1 + pi/2):
                ang1 += pi/2
            else:
                ang1 -= pi/2

            if abs(psi_s - ang2 - pi/2) < \
               abs(psi_s - ang2 + pi/2):
                ang2 += pi/2
            else:
                ang2 -= pi/2

            # Set new value of psi
            psi_s = (ang1 + ang2)/2

            # Move forward in the phi_s diraction
            cX += round(cos(psi_s) * mu)
            cY += round(sin(psi_s) * mu)

            # Record new coordinates
            uC.append([cX,cY])

    return uC
