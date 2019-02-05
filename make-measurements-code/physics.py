import numpy as np
from scipy.fftpack import fft2

#--------------
#	Image Manipulation
#--------------

def power_spectrum(image):
	''' Returns the power spectrum function of an image'''
	image = image.astype(float)
	ps_img = abs(pow(fft2(image), 2))

	return ps_img

def normalise(image):
	return image/np.amax(image)

def center(image):
	''' Re-centers an image from the origin to the middle of the display image'''
	size = image.shape
	half = int(np.floor(size[0]/2))
	image = np.roll(np.roll(image, half, 0), half, 1)
	return image

#--------------
#	Unit conversions
#--------------

def pc_to_m(x):
	''' converts from parsecs to meters'''
	return x * 3.08567758149137e16

def AU_to_m(x):
	''' converts from AU to meters'''
	return x * 1.49597870700e11

def m_to_AU(x):
	''' converts from meters to AU'''
	return x / 1.49597870700e11

def pc_to_AU(x):
	''' converts from parsecs to AU'''
	return x * 206264.80749673

def rad_to_mas(x):
	''' converts from radians to milliarcseconds'''
	return x / np.pi * 180 * 60 * 60 * 1000

def mas_to_rad(x):
	''' converts from radians to milliarcseconds'''
	return x * np.pi / 180 / 60 / 60 / 1000

def as_to_rad(x):
	''' converts from milliarcseconds to radians'''
	return x * np.pi / 180 / 60 / 60

def mum_to_m(x):
	''' coverts from micrometers to meters'''
	return x * 1e-6

def rad_to_deg(x):
	''' converts from radians to degrees'''
	return x / np.pi * 180

def sr_to_AU(x):
	''' converts from solar radius to AU'''
	return x * 0.0046524726373787	# Asteriks on this value: using what 3 conversion sites use and not what wikipedia gives.

def AU_to_sr(x):
	''' convers from AU to solar radius'''
	return x / 0.0046524726373787

#--------------
#	Unit conversions
#--------------

def angular_size(size, distance):
	''' calculates the angular size (diameter) of an object given the size of the object and the distance to it. Note that these have to be in the same units. The returned value is in radians'''

	trig_ratio = 0.5 * float(size) / distance
	angular_size = 2 * np.arctan(trig_ratio)

	return angular_size





