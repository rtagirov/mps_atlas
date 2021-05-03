import numpy as np

import matplotlib.pyplot as plt

from tqdm import tqdm

import os

import sys

def read_spec(read_mu, Nx, Ny, Nw):

#    I0m0 = np.zeros((3, 9, Nx, Ny, Nw))
#    I0m1 = np.zeros((3, 9, Nx, Ny, Nw))
#    ISm0 = np.zeros((3, 9, Nx, Ny, Nw))
#    ISm1 = np.zeros((3, 9, Nx, Ny, Nw))
#    I3m0 = np.zeros((3, 9, Nx, Ny, Nw))
#    I3m1 = np.zeros((3, 9, Nx, Ny, Nw))

    I0m0 = np.zeros((9, Nx, Ny, Nw))
    I0m1 = np.zeros((9, Nx, Ny, Nw))
    ISm0 = np.zeros((9, Nx, Ny, Nw))
    ISm1 = np.zeros((9, Nx, Ny, Nw))
    I3m0 = np.zeros((9, Nx, Ny, Nw))
    I3m1 = np.zeros((9, Nx, Ny, Nw))

    s0m0 = ['230777', '233100', '235427']
    sSm0 = ['237495', '244258', '251568']
    s3m0 = ['308000']

    s0m1 = ['004661', '007049', '009383']
    sSm1 = ['536452', '542081', '547390']
    s3m1 = ['182000']

    Ns = 3

    for m in read_mu:

        for i in tqdm(range(1, Nx + 1), desc = 'hyd, ssd: mu = ' + str(m)):

            I0m0_f = np.zeros(Ny * Nw)
            I0m1_f = np.zeros(Ny * Nw)
            ISm0_f = np.zeros(Ny * Nw)
            ISm1_f = np.zeros(Ny * Nw)

            for k in range(Ns):

                f = './mh0-hyd/' + s0m0[k] + '/1357/spec/' + str(m) + '.' + str(i)

                if os.path.exists(f): I0m0_f += np.genfromtxt(f) / Ns

            I0m0[m - 1, i - 1, :, :] = I0m0_f.reshape(Ny, Nw)

            for k in range(Ns):

                f = './mh1-hyd/' + s0m1[k] + '/1357/spec/' + str(m) + '.' + str(i)

                if os.path.exists(f): I0m1_f += np.genfromtxt(f) / Ns

            I0m1[m - 1, i - 1, :, :] = I0m1_f.reshape(Ny, Nw)

            for k in range(Ns):

                f = './mh0-ssd/' + sSm0[k] + '/1357/spec/' + str(m) + '.' + str(i)

                if os.path.exists(f): ISm0_f += np.genfromtxt(f) / Ns

            ISm0[m - 1, i - 1, :, :] = ISm0_f.reshape(Ny, Nw)

            for k in range(Ns):

                f = './mh1-ssd/' + sSm1[k] + '/1357/spec/' + str(m) + '.' + str(i)

                if os.path.exists(f): ISm1_f += np.genfromtxt(f) / Ns

            ISm1[m - 1, i - 1, :, :] = ISm1_f.reshape(Ny, Nw)

    Ns = 1

    for m in read_mu:

        for i in tqdm(range(1, Nx + 1), desc = '300: mu = ' + str(m)):

            I3m0_f = np.zeros(Ny * Nw)
            I3m1_f = np.zeros(Ny * Nw)

            for k in range(Ns):

                f = './mh0-300/' + s3m0[k] + '/1357/spec/' + str(m) + '.' + str(i)

                if os.path.exists(f): I3m0_f += np.genfromtxt(f) / Ns

            I3m0[m - 1, i - 1, :, :] = I3m0_f.reshape(Ny, Nw)

            for k in range(Ns):

                f = './mh1-300/' + s3m1[k] + '/1357/spec/' + str(m) + '.' + str(i)

                if os.path.exists(f): I3m1_f += np.genfromtxt(f) / Ns

            I3m1[m - 1, i - 1, :, :] = I3m1_f.reshape(Ny, Nw)

    return I0m0, I0m1, ISm0, ISm1, I3m0, I3m1

def calc_mean(I, m, j):

    nz = np.where((I[m, :, :, j] > 0.0) & (I[m, :, :, j] < 1.0) & (~np.isnan(I[m, :, :, j])))

    Im = 0.0

    for i in range(len(nz[0])):

        Im += I[m, nz[0][i], nz[1][i], j]

    Im /= len(nz[0])

    return Im

Nx = 512
Ny = 512
Nw = 1221

#read_mu = np.array([1, 3, 5, 7]).astype(int)
read_mu = np.array([1]).astype(int)

#I0m0, I0m1, ISm0, ISm1, I3m0, I3m1 = read_spec(read_mu, Nx, Ny, Nw)

#np.savez('met_spec', I0m0 = I0m0, I0m1 = I0m1, ISm0 = ISm0, ISm1 = ISm1, I3m0 = I3m0, I3m1 = I3m1)

print('I0m0')
#I0m0 = np.load('met_spec.npz')['I0m0']
I0m0 = np.load('I0m0.npz')['I']
#np.savez('I0m0', I = I0m0)

print('I0m1')
#I0m1 = np.load('met_spec.npz')['I0m1']
I0m1 = np.load('I0m1.npz')['I']
#np.savez('I0m1', I = I0m1)

print('ISm0')
#ISm0 = np.load('met_spec.npz')['ISm0']
ISm0 = np.load('ISm0.npz')['I']
#np.savez('ISm0', I = ISm0)

print('ISm1')
#ISm1 = np.load('met_spec.npz')['ISm1']
ISm1 = np.load('ISm1.npz')['I']
#np.savez('ISm1', I = ISm1)

print('I3m0')
#I3m0 = np.load('met_spec.npz')['I3m0']
I3m0 = np.load('I3m0.npz')['I']
#np.savez('I3m0', I = I3m0)

print('I3m1')
#I3m1 = np.load('met_spec.npz')['I3m1']
I3m1 = np.load('I3m1.npz')['I']
#np.savez('I3m1', I = I3m1)

#sys.exit()

I0m0_m = np.zeros((9, Nw))
I0m1_m = np.zeros((9, Nw))
ISm0_m = np.zeros((9, Nw))
ISm1_m = np.zeros((9, Nw))
I3m0_m = np.zeros((9, Nw))
I3m1_m = np.zeros((9, Nw))

for m in read_mu:

    for j in tqdm(range(Nw), desc = 'mu = ' + str(m)):

        I0m0_m[m - 1, j] = calc_mean(I0m0, m - 1, j)
        I0m1_m[m - 1, j] = calc_mean(I0m1, m - 1, j)

        ISm0_m[m - 1, j] = calc_mean(ISm0, m - 1, j)
        ISm1_m[m - 1, j] = calc_mean(ISm1, m - 1, j)

        I3m0_m[m - 1, j] = calc_mean(I3m0, m - 1, j)
        I3m1_m[m - 1, j] = calc_mean(I3m1, m - 1, j)

np.savez('met_spec_m', I0m0_m = I0m0_m, I0m1_m = I0m1_m, ISm0_m = ISm0_m, ISm1_m = ISm1_m, I3m0_m = I3m0_m, I3m1_m = I3m1_m)
