import numpy as np

import matplotlib.pyplot as plt

from tqdm import tqdm

import os

import sys

def calc_mean(I):

    nz = np.where((I[:, :] > 0.0) & (I[:, :] < 1.0) & (~np.isnan(I[:, :])))

    Im = 0.0

    for i in range(len(nz[0])):

        Im += I[nz[0][i], nz[1][i]]

    Im /= len(nz[0])

    return Im

def read_slice(cube, snapshot, mu, slice_num, Ny, Nw):

    I_slice = np.zeros(Ny * Nw)

    f = './' + cube + '/' + snapshot + '/1357/spec/' + str(mu) + '.' + str(slice_num)

    if os.path.exists(f): I_slice = np.genfromtxt(f)

    return I_slice.reshape(Ny, Nw)

def read_spec(cube, snapshot, read_mu, Nx, Ny, Nw):

    I = np.zeros((Nx, Ny, Nw))

    I_averaged = np.zeros((9, Nw))

    for m, mu in enumerate(read_mu):

        for i in tqdm(range(Nx), desc = cube + ', ' + snapshot + ', mu = ' + str(mu)):
#        for i in range(Nx):

            I[i, :, :] = read_slice(cube, snapshot, mu, i + 1, Ny, Nw)

        for j in tqdm(range(Nw), desc = cube + ', ' + snapshot + ', mu = ' + str(mu)):
#        for j in range(Nw):

            I_averaged[mu - 1, j] = calc_mean(I[:, :, j])

    return I_averaged

#s0m0 = ['230777', '233100', '235427']
#s0m1 = ['004661', '007049', '009383']
#sSm0 = ['237495', '244258', '251568']
#sSm1 = ['536452', '542081', '547390']
#s3m0 = ['308000']
#s3m1 = ['182000']

Nx = 512
Ny = 512
Nw = 1221

read_mu = np.array([1, 3, 5, 7]).astype(int)

cube_snapshot = sys.argv[1]

if not os.path.exists(cube_snapshot):

    print('Directory ' + cube_snapshot + ' does not exist.')

    sys.exit(1)

cube = cube_snapshot.split('/')[0]
snapshot = cube_snapshot.split('/')[1]

I = read_spec(cube, snapshot, read_mu, Nx, Ny, Nw)

np.savez('./npz/' + cube + '_' + snapshot, I = I)