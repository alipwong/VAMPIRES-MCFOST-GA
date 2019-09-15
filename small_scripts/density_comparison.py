# determining whether you specify a number density or mass density
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

folder_1 = "default"
folder_2 = "density_file"
#
# folder_1 = "3_10_-5e-1-8"
# folder_2 = "3_10_-5e-1-8-m1000"

# set image boundaries
x1 = 475
x2 = 525

img1 = fits.open("../tests/logic/{}/RT.fits".format(folder_1))[0].data
# default_dm = fits.open("../{}/dust_mass_density.fits".format(folder_1))[0].data

img2 = fits.open("../tests/logic/{}/RT.fits".format(folder_2))[0].data
# ones_dm = fits.open("../{}/dust_mass_density.fits".format(folder_2))[0].data
# for i in range(3):
#     img1[i][0][0] = img1[i][0][0]/np.max(img1[i][0][0])
#     img2[i][0][0] = img2[i][0][0]/np.max(img2[i][0][0])
# print(np.max(img1))

plt.figure(figsize=(12,9))
for i in range(3):
    plt.subplot(3, 4, i + 1)
    plt.imshow(img1[i][0][0][x1:x2, x1:x2])
    plt.colorbar()
for i in range(3):
    plt.subplot(3, 4, i + 5)
    plt.imshow(img2[i][0][0][x1:x2, x1:x2])
    plt.colorbar()
for i in range(3):
    plt.subplot(3, 4, i + 9)
    plt.imshow((img1[i][0][0] - img2[i][0][0])[x1:x2, x1:x2])
    plt.colorbar()
# plt.subplot(3, 4, 4)
# plt.imshow(default_dm)
# plt.colorbar()
# plt.subplot(3, 4, 8)
# plt.imshow(ones_dm)
# plt.colorbar()
# plt.subplot(3, 4, 12)
# plt.imshow(abs(default_dm-ones_dm))
# plt.colorbar()
plt.savefig("../tests/logic/images/{}-{}.png".format(folder_1, folder_2))
