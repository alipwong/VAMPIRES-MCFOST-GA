'''
Get a data object with
    diffData = vampDiffdata(dataFilename, cubeInfoFilename)
Which has the following data attributes:
    .vhhv       - polarised differential visibilities for instrument Stokes Q
    .vhhverr    - 1-sigma errors for polarised differential visibilities for instrument Stokes Q
    .vhhvu      - polarised differential visibilities for instrument Stokes U
    .vhhvuerr   - 1-sigma errors for polarised differential visibilities for instrument Stokes U
    .diffCP     - polarised differential closure-phases for instrument Stokes Q
    .diffCPerr  - 1-sigma errors for polarised differential closure-phases for instrument Stokes Q
    .diffCPu    - polarised differential closure-phases for instrument Stokes U
    .diffCPuerr - 1-sigma errors for polarised differential closure-phases for instrument Stokes U
    .blengths   - baseline lengths
    .bazims     - baseline azimuths (in instrument coords)
    .u_coords   - u coordinates (in instrument coords)
    .v_coords   - v coordinates (in instrument coords)

It also has attributes containing geometry data, as follows:
    .BL2H_IX     - The 'baseline-to-hole indices'. A (numBLs x 2) array which specifies the pair of hole
                        numbers corresponding to each baseline.
                        I.e., given a baseline, bl2h_ix gives the 2 holes that go to make it up
    .H2BL_IX     - The 'hole-to-baseline indices'. An (numHoles x numHoles) array which baselines correspond
                        to any pair of holes.
                        I.e., given a pair of holes i,j H2BL_IX[i,j] gives the number of the baseline
    .BS2BL_IX    - The 'bispectrum-to-baseline' indices. Relates bispectrum/closure phase triangle indices
                        to baselines. I.e., given a point in the bispectrum, BS2BL_IX gives the 3 baselines
                        which make the triangle.
    .BL2BS_IX    - The 'baseline-to-bispectrum indexes'. BL2BS_IX gives the index of all points in the
                        bispectrum containing a given baseline
    .FILTER      - 2-element array giving the centre wavelength and bandpass of the filter

and metadata:
    .UTCs        - UTC times for each observation
    .ras, .decs  - telescope RA and dec for each observation
    .mask        - name of aperture mask used
    .emgains     - EM gain of camera for each observation
    .mffile      - the name of the matched-filter file used in data reduction. Thie file contains
                   information on u,v sampling, correspondence of closure phase triangles to baselines, etc.
    .pkflux      - peak flux in frame for each observation
    .totflux     - total flux in frame for each observation

'''

from scipy import io
import numpy as np

class VAMPIRES_data:

	def __init__(self, filename, cube_info_filename):
		''' Reads in the polarised-differential calibrated VAMPIRES data produced by the IDL pipeline:param filename: The full filename of the data file. E.g. diffdata_vega_.......idlvar '''

		data_object = io.readsav(filename, python_dict=False, verbose=False)
		self.vhvv = np.copy(data_object.vhvv)
		self.vhvverr = data_object.vhvverr
		self.vhvvu = np.copy(data_object.vhvvu)
		self.vhvvuerr = data_object.vhvvuerr
		self.blengths = data_object.blengths
		self.bazims = data_object.bazims
		self.inFilename = filename
		try:
			self.diffCP = data_object.cp
			self.diffCPerr = data_object.cperr
			self.diffCPu = data_object.cpu
			self.diffCPuerr = data_object.cpuerr
			self.BL2H_IX = data_object.BL2H_IX
			self.H2BL_IX = data_object.H2BL_IX
			self.BL2BS_IX = data_object.BL2BS_IX
			self.BS2BL_IX = data_object.BS2BL_IX
		except:
			print("Couldn't find diff CP data for " + filename)
		try:
			self.u_coords = data_object.u_coords
			self.v_coords = data_object.v_coords
		except:
			# Some older files didn't have these saved
			self.u_coords = ()
			self.v_coords = ()
		del (data_object)

		# Get useful metadata from cubeinfo file
		cube_info_object = io.readsav(cube_info_filename, python_dict=False, verbose=False)
		self.UTCs = cube_info_object.olog.utc[0]
		self.filters = cube_info_object.olog.filter[0]
		self.ras = cube_info_object.olog.ra[0]
		self.decs = cube_info_object.olog.dec[0]
		self.mask = cube_info_object.olog.mask[0]
		self.adate = cube_info_object.olog.adate[0]
		self.emgains = cube_info_object.olog.emgain[0]
		self.mffile = cube_info_object.plog.mf_file[0]
		self.pkflux = cube_info_object.framestats.pkflx[0]
		self.totflux = cube_info_object.framestats.totflx[0]
		self.cubename = cube_info_object.olog.cube_fname[0][0]
		del (cube_info_object)










