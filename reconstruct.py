""" Weighted back projection algorithm. """
from __future__ import print_function
from MRCFile import MRCFile

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
import matplotlib.pyplot as plt

############ THE ALGORITHM ###############
# Given N images (I_1 ... I_N), viewing directions (F_1 ... F_N), and distance D = 153
# find an approximation to rho defined as follow:
# rho(x,y,z) = IFFT (B(x,y,z) / H(x,y,z)). 
#
# Here B is the back projection of the images, and H is some filtering.
#
# To compute B, for each image, compute its back projection as done in b(image) below, and sum them.
# To compute H, for each orientation, compute D * sinc(D pi (x,y,z) dot c_i) where c_i = last
# column of the orientation matrix. Then sum them all. 
##########################################

# CAVEATS 
# 1) Whenever we use FFTN (or IFFTN), we must perform a fourier shift fftshift. I.e. Shifting the 
# zero-frequency component to the center of the spectrum (2D or 3D matrix in our case)

# after reconstruct, write an MRC file, and load that in the kimera to look (from UCSF)?

def b(image): 
    image_hat = fft.fftn(image) / 153.0 / 153.0
    image_hat = fft.fftshift(image_hat) 
    image_hat = np.expand_dims(image_hat, axis=2)
    image_hat = np.repeat(image_hat, 153, axis=2) #turn 2D image to 153 copies of 2D image

    z = np.linspace(-1, 1, 153)
    s = 2 * np.sinc(2 * math.pi * z)
    s = np.expand_dims(s, axis=0)
    s = np.expand_dims(s, axis=0)
    s = np.repeat(s, 153, axis=0)
    s = np.repeat(s, 153, axis=1)

    result = fft.ifftshift(image_hat * s)
    return fft.ifftn(result)

def H(orientations): # This is filtering
    N = len(orientations)
    H_matrix = np.zeros((153, 153, 153), dtype=complex)

    x = np.linspace(-1, 1, 153)
    y = np.linspace(-1, 1, 153)
    z = np.linspace(-1, 1, 153)
    xx, yy, zz = np.meshgrid(x, y, z)

    for i in range(N):
        last_col = orientations[i][:, 2]
        H_matrix += 153.0 * np.sinc(153.0 * math.pi * (xx * last_col[0] + yy * last_col[1] + zz * last_col[2]) )
    return H_matrix

    # for i in range(153):
    #     print("Computing H: currently on the " + str(i) + "th loop")
    #     for j in range(153):
    #         for k in range(153):
    #             sum = 0.0

    #             for l in range(N):
    #                 c = orientations[l][:, 2]
    #                 v1 = [x[i], y[j], z[k]]
    #                 dotted = np.dot(v1, c)
    #                 sum += 153 * np.sinc(153 * math.pi * dotted)

    #             H_matrix[i][j][k] = 1.0 / sum

    # return H_matrix

def make_B_matrix(images):
    rho_matrix_sum = np.zeros((153, 153, 153), dtype=complex)
    for image in images:
        rho_matrix_sum += b(image) 
        print("Computed some b for some image! Continue running....")

    return rho_matrix_sum

def reconstruct(images, orientations): # input should be lists of matrices
    B_matrix = make_B_matrix(images)
    H_matrix = H(orientations)
    result = fft.ifftshift(B_matrix / H_matrix)

    return fft.ifft(result)

def main():
    f = MRCFile('zika_153.mrc')
    f.load_all_slices()

    r = np.array([[0.25581,  -0.77351,   0.57986], 
     [-0.85333,  -0.46255,  -0.24057], 
     [0.45429,  -0.43327,  -0.77839]])


    molecule = f.slices[:, :, :]
    middle_slice = fft.ifftn(fft.ifftshift(fft.fftshift(fft.fftn(molecule))[:, :, 76]))
    # middle_slice = middle_slice.real
    # print(middle_slice)
    # plt.imshow(middle_slice)
    # plt.show()


    middle_back_proj = b(middle_slice)
    middle_back_proj = middle_back_proj.real[:, :, 76]
    plt.imshow(middle_back_proj)
    plt.show()


    # display_image(f.slices[:,:,0])


    # s = r[:,2]
    # print(s)

    # matrix_dotted = np.full((3,3,3), 0.17)
    # print(matrix_dotted * r)

    # r = np.expand_dims(r, axis=2)
    # r = np.repeat(r, 153, axis=2) #turn 2D image to 153 copies of 2D image
    # print(r.shape)
    # print(r[:, :, 133] - r[:, :, 46])

    # print(f.slices[:,:,0].shape)
    # plt.imshow(f.slices[:, :, 68], cmap='hot')
    # plt.show()

    # fdsa = b(f.slices[:, :, 76])
    # fdsa = fdsa.real
    # plt.imshow(fdsa[:, :, 100], cmap='hot')
    # plt.show()

    # orientations = []
    # for i in range(153):
    #     orientations.append(r)

    # images = []
    # for i in range(153):
    #     images.append(f.slices[:, :, i])

    # fdsa = reconstruct(images, orientations)
    # fdsa = fdsa.real
    # plt.imshow(fdsa[:, :, 100], cmap='hot')
    # plt.show()


main()