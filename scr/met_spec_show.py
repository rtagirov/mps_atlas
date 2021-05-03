import numpy as np

import matplotlib.pyplot as plt

from tqdm import tqdm

import os

print('1')
I0m0 = np.load('met_spec.npz')['I0m0']
print('2')
I0m1 = np.load('met_spec.npz')['I0m1']
print('3')
ISm0 = np.load('met_spec.npz')['ISm0']
print('4')
ISm1 = np.load('met_spec.npz')['ISm1']
print('5')
I3m0 = np.load('met_spec.npz')['I3m0']
print('6')
I3m1 = np.load('met_spec.npz')['I3m1']

#I3m0[1, :, :, :] = I3m0[0, :, :, :]
#I3m0[2, :, :, :] = I3m0[0, :, :, :]
#I3m1[1, :, :, :] = I3m0[0, :, :, :]
#I3m1[2, :, :, :] = I3m0[0, :, :, :]

I0m0_m = np.zeros((9, Nw_atl))
I0m1_m = np.zeros((9, Nw_atl))
ISm0_m = np.zeros((9, Nw_atl))
ISm1_m = np.zeros((9, Nw_atl))
I3m0_m = np.zeros((9, Nw_atl))
I3m1_m = np.zeros((9, Nw_atl))

#read_mu = np.array([1, 3, 5, 7]).astype(int)
read_mu = np.array([1]).astype(int)

Nw = 1221

for m in read_mu:

   for j in range(Nw):

       nz0a = np.where((I0a[m - 1, :, :, j] > 0.0) & (~np.isnan(I0a[m - 1, :, :, j])))[0]
       nzSa = np.where((ISa[m - 1, :, :, j] > 0.0) & (~np.isnan(ISa[m - 1, :, :, j])))[0]
       nz3a = np.where((I3a[m - 1, :, :, j] > 0.0) & (~np.isnan(I3a[m - 1, :, :, j])))[0]
