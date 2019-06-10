# creates density arrays
# described by where you put the 1's

import numpy as np

def full(nz, n_az, n_grains):
    ones = np.ones((1, nz, n_az, n_grains ))
    return ones

def power(n_az, n_grains, min = 1, max = 2, power = -0.5):
    r = np.linspace(min, max, n_grains)
    row = list(map(lambda x: pow((x + 1), power), r))
    x = np.array(row * n_az).reshape(n_az, n_grains)
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
    zeros_layer = np.zeros((n_az, n_grains))
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
    nz = np.array([zeros, zeros])
    return np.expand_dims(nz, axis=0)