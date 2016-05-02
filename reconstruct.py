""" Weighted back projection algorithm. """
from __future__ import print_function

import sys
import struct
import timeit
import scipy as sp
from scipy.interpolate import griddata
import numpy as np
import numpy.fft as fft
import os.path
from collections import namedtuple
import math

############ THE ALGORITHM ###############
# Given N images (I_1 ... I_N), viewing directions (F_1 ... F_N), and distance D
# find an approximation to rho defined as follow:
# rho(x,y,z) = IFFT (B(x,y,z) / H(x,y,z)).
#
# B(x,y,z) = sum i from 1 to N of FFT(b_i * [(x,y,z) dot a_i + (x,y,z) dot b_i] + (x,y,z) dot c_i)
# (a_i b_i c_i are the columns of rotation matrix F_i)
# H(x,y,z) = sum i from 1 to N of D * sinc(D * pi * [(x,y,z) dot c_i])
# 
# b_i(x_i,y_i,z_i) = IFFT( FFT(I) * FFT(l_i) ) (this b_i is a separate function, not b_i as above)
# FFT(l_i) = FFT(delta(x,y) * rect(z)) = D * sinc(D * pi * z_i) by some clever math.
# sinc(x) = sin(x) / x
# rect(z) = 1 if -D/2 < z < D/2 and 0 otherwise
# x_i = (x,y,z) dot a_i, and likewise for y_i and z_i
##########################################

# NOTE: Whenever we use fft, we must perform a fourier shift fftshift. I.e. Shifting the 
# zero-frequency component to the center of the spectrum (matrix or vector in our case)

def fourier_l(z_i, D):
    return D * sinc(D * math.pi * z_i)

def fourier_image(image):
    my_image = fft.fftn(image)
    my_image = fft.fftshift(my_image)
    my_image = my_image[:None]
    return my_image

def b(x_i, y_i, z_i, image, D): #image = 1 image
    fft_image = fourier_image(x_i, y_i, image)
    fft_l = fourier_l(z_i, D)
    result = fft.ifftn(fft_image * fft_l)

    return fft.ifftshift(result)

def B(x, y, z, images, orientations, D):
    N = len(images)
    sum = 0.0
    for i in range(0, N):
        a_i = [row[0] for row in orientations[i]]
        b_i = [row[1] for row in orientations[i]]
        c_i = [row[2] for row in orientations[i]]

        v1 = [x, y, z]
        v2 = [a_i, b_i, c_i]
        dotted = dot(v1, v2)

        x_i = dot(v1, a_i) # these are scalars
        y_i = dot(v1, b_i)
        z_i = dot(v1, c_i)
        my_b = b(x_i, y_i, z_i, images[i], D)
        result = fft.fftn(my_b * dotted)
        sum += fft.fftshift(result)
    return sum

def sinc(x):
    return math.sin(x) / x

def H(x, y, z, orientations, D):
    N = len(orientations)
    sum = 0.0
    for j in range(0, N):
        c = [row[2] for row in orientations[j]]
        v1 = [x, y, z]
        dotted = dot(v1, c)
        sum += D * sinc(D * math.pi * dotted)
    return sum

def dot(v1, v2):
    assert len(v1) == 3, 'v1 is not of length 3, is that right?'
    assert len(v2) == 3, 'v2 is not of length 3, is that right?'
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

def reconstruct(images, orientations, D):
    assert len(images) == len(orientations), '# of images not the same as # of orientations'
    size = len(images)

    rho = np.zeros((size, size, size))
    for x in range(size):
        for y in range(size):
            for z in range(size):
                B = B(x, y, z, images, orientations, D)
                H = H(x, y, z, orientations, D)
                result = fft.ifftn(B / H)
                rho[x][y][z] = fft.ifftshift(result)
    return rho

def test():
    r = [[0.25581,  -0.77351,   0.57986], 
     [-0.85333,  -0.46255,  -0.24057], 
     [0.45429,  -0.43327,  -0.77839]]
    a_i = [row[0] for row in r]
    b_i = [row[1] for row in r]
    c_i = [row[2] for row in r]

    rho = np.zeros((3, 3))
    rho = np.expand_dims(rho, axis=2)
    print(rho.shape)
    rho = np.repeat(rho, 3, axis=2)
    print(rho.shape)


test()