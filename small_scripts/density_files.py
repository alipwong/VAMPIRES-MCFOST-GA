# creates density arrays
# described by where you put the 1's

import numpy as np
from astropy.io import fits
import os

def power_grid(p = -0.5):
    print(os.getcwd())
    fn = "../make-measurements/grid.fits.gz"
    hdul = fits.open(fn)
    # r = hdul[0].data[0,0,:,:]
    # x = np.power(r, p)
    r = hdul[0].data[:,:,:,:]
    # x = np.expand_dims(x, axis=0)
    # return np.expand_dims(x, axis=0)
    return r

def r_power_grid():
    fn = "../tests/grid_files/2D_default/grid.fits.gz"

    # fn = "../tests/grid_files/ngrains2/grid.fits.gz"
    hdul = fits.open(fn)
    x = hdul[0].data[0,0,:,:]
    y = hdul[0].data[1,0,:,:]
    r = np.power(np.power(x, 2) + np.power(y, 2), 0.5)
    r = np.power(r, -2)
    r = np.expand_dims(r, axis=0)
    r = np.expand_dims(r, axis=0)
    return r

def grid_3d():
    fn = "../tests/grid_files/3D_default_naz_200_ngrains_5/grid.fits.gz"
    hdul = fits.open(fn)
    r = hdul[0].data[:,:,:,:]
    x = hdul[0].data[0,:,:,:]
    y = hdul[0].data[1,:,:,:]
    r = np.power(np.power(x, 2) + np.power(y, 2), 0.5)
    # r = hdul[0].data[0,0,:,:]
    # x = np.power(r, p)
    r = np.power(r, -0.5)
    r = np.expand_dims(r, axis=0)
    # r = np.array([r, r, r, r, r])
    return r

def dim_3d():
    ones = np.ones((1, 2, 140, 100))
    return ones

def full(nz, n_az, n_grains):
    ones = np.ones((1, nz, n_az, n_grains ))
    return ones

def power(n_az, n_grains, min = 1, max = 2, power = -0.5):

    r = np.linspace(min, max, n_grains)


    row = list(map(lambda x: pow((x + 1), power), r))
    x = np.array(row * n_az).reshape(n_az, n_grains)
    x = np.power(r, power)
    x = np.expand_dims(x, axis=0)
    return np.expand_dims(x, axis=0)

def power_shell(n_az, n_grains, min = 1, max = 2, power = -0.5):
    zeros = int(np.floor(min / max * n_grains))
    nonzeros = n_grains - zeros
    r = np.array(np.linspace(min, max, nonzeros)).reshape(1, nonzeros)
    row = np.concatenate((np.zeros((1, zeros)), r), axis = 1)
    x = np.repeat(row, n_az, axis=0)
    x = np.expand_dims(x, axis=0)
    return np.expand_dims(x, axis=0)

def bottom(n_rad, nz, n_az, n_grains):
    ones = np.ones((n_rad, nz, int(np.ceil(n_az/2)), n_grains))
    zeros = np.zeros((n_rad, nz, int(np.floor(n_az/2)), n_grains))
    return np.concatenate((zeros, ones), axis=2)

def top(n_rad, nz, n_az, n_grains):
    ones = np.ones((n_rad, nz, int(np.ceil(n_az/2)), n_grains))
    zeros = np.zeros((n_rad, nz, int(np.floor(n_az/2)), n_grains))
    return np.concatenate((ones, zeros), axis=2)

def right(n_az, n_grains):
    # this shouldn't affect anything because the 1's and 0's will just define which grains there are
    zeros = np.zeros((n_az, int(np.floor(n_grains/2))))
    ones = np.ones((n_az, int(np.ceil(n_grains/2))))
    x = np.hstack((zeros, ones))
    x = np.expand_dims(x, axis = 0)
    return np.expand_dims(x, axis = 0)

def left(n_az, n_grains):
    # this shouldn't affect anything because the 1's and 0's will just define which grains there are
    zeros = np.zeros((n_az, int(np.floor(n_grains/2))))
    ones = np.ones((n_az, int(np.ceil(n_grains/2))))
    x = np.hstack((ones, zeros))
    x = np.expand_dims(x, axis = 0)
    return np.expand_dims(x, axis = 0)

def layer_1(n_az, n_grains):
    zeros = np.zeros((n_az, n_grains))
    ones = np.ones((n_az, n_grains))
    nz = np.array([ones, zeros])
    return np.expand_dims(nz, axis=0)

def layer_2(n_az, n_grains):
    zeros = np.zeros((n_az, n_grains))
    ones = np.ones((n_az, n_grains))
    twos = np.ones((n_az, n_grains)) * 2
    nz = np.array([zeros, ones, twos])
    return np.expand_dims(nz, axis=0)

def cross_section_top(n_az, n_grains):
    ones = np.ones((int(np.ceil(n_az / 2)), n_grains))
    zeros = np.zeros((int(np.floor(n_az / 2)), n_grains))
    top_layer = np.concatenate((ones, zeros), axis=0)

    #
    zeros_layer = np.zeros((n_az, n_grains))
    nz = np.array([top_layer, zeros_layer])
    return np.expand_dims(nz, axis=0)

def cross_section_bottom(n_az, n_grains):
    ones = np.ones((int(np.ceil(n_az / 2)), n_grains))
    zeros = np.zeros((int(np.floor(n_az / 2)), n_grains))
    top_layer = np.concatenate((zeros, ones), axis=0)

    #
    zeros_layer = np.zetros((n_az, n_grains))
    nz = np.array([top_layer, zeros_layer])
    return np.expand_dims(nz, axis=0)

def grain_size(n_az, n_grains):
    zeros = np.zeros((n_az, n_grains))
    ones = np.ones((n_az, n_grains))
    grain_1 = np.array([ones, zeros])
    grain_2 = np.array([zeros, ones])

    return  np.array([grain_1, grain_2])

def ring(n_az, n_grains, thickness = 0.1, height = 0.1):
    zeros = np.zeros((n_az, n_grains))
    h = int(np.ceil(n_az * height))
    t = int(np.ceil(n_grains * thickness))
    for i in range(h):
        for j in range(t):
            zeros[i][-(j + 1)] = 1
    # nz = np.array([zeros, zeros])
    nz = np.array([zeros])
    return np.expand_dims(nz, axis=0)



