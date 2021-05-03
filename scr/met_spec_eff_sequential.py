import numpy as np

import matplotlib.pyplot as plt

from tqdm import tqdm

import os

import sys

def calc_mean(I, m, j):

    nz = np.where((I[m, :, :, j] > 0.0) & (I[m, :, :, j] < 1.0) & (~np.isnan(I[m, :, :, j])))

    Im = 0.0

    for i in range(len(nz[0])):

        Im += I[m, nz[0][i], nz[1][i], j]

    Im /= len(nz[0])

    return Im

def read_slice(cube, snapshot, mu, slice_num, Ny, Nw):

    I_slice = np.zeros(Ny * Nw)

    f = './' + cube + '/' + snapshot + '/1357/spec/' + str(mu) + '.' + str(slice_num)

    if os.path.exists(f): I_slice = np.genfromtxt(f)

    return I_slice.reshape(Ny, Nw)

def read_spec(cube, snapshots, read_mu, Nx, Ny, Nw):

    Nm = len(read_mu)

    I = np.zeros((Nm, Nx, Ny, Nw))

    I_averaged = np.zeros((9, Nw))

    Ns = len(snapshots)

    for k in range(Ns):

        for m, mu in enumerate(read_mu):

#            for i in tqdm(range(Nx), desc = cube + ', reading k = ' + str(k + 1) + ', mu = ' + str(mu)):
            for i in range(Nx):

                I[m, i, :, :] = read_slice(cube, snapshots[k], mu, i + 1, Ny, Nw)

#            for j in tqdm(range(Nw), desc = cube + ', averaging'):
            for j in range(Nw):

                I_averaged[mu - 1, j] += calc_mean(I, m, j)

    I_averaged /= Ns

    return I_averaged

s0m0 = ['230777', '233100', '235427']
s0m1 = ['004661', '007049', '009383']
sSm0 = ['237495', '244258', '251568']
sSm1 = ['536452', '542081', '547390']
s3m0 = ['308000']
s3m1 = ['182000']

Nx = 512
Ny = 512
Nw = 1221

read_mu = np.array([1, 3, 5, 7]).astype(int)
#read_mu = np.array([1]).astype(int)

print('I0m0')
I0m0 = read_spec('mh0-hyd', s0m0, read_mu, Nx, Ny, Nw)
np.savez('I0m0', I = I0m0)

print('I0m1')
I0m1 = read_spec('mh1-hyd', s0m1, read_mu, Nx, Ny, Nw)
np.savez('I0m1', I = I0m1)

print('ISm0')
ISm0 = read_spec('mh0-ssd', sSm0, read_mu, Nx, Ny, Nw)
np.savez('ISm0', I = ISm0)

print('ISm1')
ISm1 = read_spec('mh1-ssd', sSm1, read_mu, Nx, Ny, Nw)
np.savez('ISm1', I = ISm1)

print('I3m0')
I3m0 = read_spec('mh0-300', s3m0, read_mu, Nx, Ny, Nw)
np.savez('I3m0', I = I3m0)

print('I3m1')
I3m1 = read_spec('mh1-300', s3m1, read_mu, Nx, Ny, Nw)
np.savez('I3m1', I = I3m1)
