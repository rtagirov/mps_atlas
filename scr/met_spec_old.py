import numpy as np

import matplotlib.pyplot as plt

from tqdm import tqdm

import os

def read_spec(read_mu, Nx, Ny, Nw):

    I0m0 = np.zeros((3, 9, Nx, Ny, Nw))
    I0m1 = np.zeros((3, 9, Nx, Ny, Nw))
    ISm0 = np.zeros((3, 9, Nx, Ny, Nw))
    ISm1 = np.zeros((3, 9, Nx, Ny, Nw))
    I3m0 = np.zeros((3, 9, Nx, Ny, Nw))
    I3m1 = np.zeros((3, 9, Nx, Ny, Nw))

#    I0m0 = np.zeros((9, Nx, Ny, Nw))
#    I0m1 = np.zeros((9, Nx, Ny, Nw))
#    ISm0 = np.zeros((9, Nx, Ny, Nw))
#    ISm1 = np.zeros((9, Nx, Ny, Nw))
#    I3m0 = np.zeros((9, Nx, Ny, Nw))
#    I3m1 = np.zeros((9, Nx, Ny, Nw))

    s0m0 = ['230777', '233100', '235427']
    sSm0 = ['237495', '244258', '251568']
    s3m0 = ['308000']

    s0m1 = ['004661', '007049', '009383']
    sSm1 = ['536452', '542081', '547390']
    s3m1 = ['182000']

    Ns = 3

    for m in read_mu:

        for i in tqdm(range(1, Nx + 1), desc = 'hyd, ssd: mu = ' + str(m)):

            for k in range(Ns):

                f = './mh0-hyd/' + s0m0[k] + '/1357/spec/' + str(m) + '.' + str(i)

                if os.path.exists(f):

                    I0m0_f = np.genfromtxt(f)

                    I0m0[k, m - 1, i - 1, :, :] = I0m0_f.reshape(Ny, Nw)

                f = './mh1-hyd/' + s0m1[k] + '/1357/spec/' + str(m) + '.' + str(i)

                if os.path.exists(f):

                    I0m1_f = np.genfromtxt(f)

                    I0m1[k, m - 1, i - 1, :, :] = I0m1_f.reshape(Ny, Nw)

                f = './mh0-ssd/' + sSm0[k] + '/1357/spec/' + str(m) + '.' + str(i)

                if os.path.exists(f):

                    ISm0_f = np.genfromtxt(f)

                    ISm0[k, m - 1, i - 1, :, :] = ISm0_f.reshape(Ny, Nw)

                f = './mh1-ssd/' + sSm1[k] + '/1357/spec/' + str(m) + '.' + str(i)

                if os.path.exists(f):

                    ISm1_f = np.genfromtxt(f)

                    ISm1[k, m - 1, i - 1, :, :] = ISm1_f.reshape(Ny, Nw)

    Ns = 1

    for m in read_mu:

        for i in tqdm(range(1, Nx + 1), desc = '300: mu = ' + str(m)):

            for k in range(Ns):

                f = './mh0-300/' + s3m0[k] + '/1357/spec/' + str(m) + '.' + str(i)

                if os.path.exists(f):

                    I3m0_f = np.genfromtxt(f)

                    I3m0[k, m - 1, i - 1, :, :] = I3m0_f.reshape(Ny, Nw)

                f = './mh1-300/' + s3m1[k] + '/1357/spec/' + str(m) + '.' + str(i)

                if os.path.exists(f):

                    I3m1_f = np.genfromtxt(f)

                    I3m1[k, m - 1, i - 1, :, :] = I3m1_f.reshape(Ny, Nw)

    return I0m0, I0m1, ISm0, ISm1, I3m0, I3m1

Nx = 512
Ny = 512
Nw = 1221

#read_mu = np.array([1, 3, 5, 7]).astype(int)
read_mu = np.array([1]).astype(int)

I0m0, I0m1, ISm0, ISm1, I3m0, I3m1 = read_spec(read_mu, Nx, Ny, Nw)

np.savez('met_spec', I0m0 = I0m0, I0m1 = I0m1, ISm0 = ISm0, ISm1 = ISm1, I3m0 = I3m0, I3m1 = I3m1)
