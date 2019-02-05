from physics import *
from run_MCFOST import *
from scipy import interpolate
import matplotlib.pyplot as plot
import matplotlib.cm as cm


class Star:

    def __init__(self, free_parameters, default_params, WAVELENGTH, MCFOST_path, obs_data, load=False, verbose=False, scale_fact = False):

        self.free_parameters = free_parameters
        self.n_free_parameters = len(free_parameters)
        self.default_parameters = default_params
        self.WAVELENGTH = WAVELENGTH
        self.MCFOST_path = MCFOST_path
        self.obs_data = obs_data
        self.verbose = verbose
        self.scale_fact = scale_fact

        if load:
            self.data = fits.getdata(self.MCFOST_path + 'RT.fits')
        else:
            self.data, self.time = build_model(self.free_parameters, self.default_parameters, WAVELENGTH, MCFOST_path,
                                    verbose=self.verbose, scale_fact = self.scale_fact)

        if self.data is None:
            print("...\nStar not made\n...")
            self.made = False
        else:
            self.made = True

        if self.made:
            self.parameters = pickle.load(open(self.MCFOST_path + "star_parameters", "rb"))

            self.observe()

            # calculate some other useful constants
            self.distance_m = pc_to_m(float(self.parameters.distance))
            self.radius_m = AU_to_m(sr_to_AU(float(self.parameters.radius)))
            self.rin_m = AU_to_m(float(self.parameters.Rin))
            self.rout_m = AU_to_m(float(self.parameters.Rout))

            self.diameter_mas = rad_to_mas(angular_size(2 * self.radius_m, self.distance_m))
            self.din_mas = rad_to_mas(angular_size(2 * self.rin_m, self.distance_m))
            self.dout_mas = rad_to_mas(angular_size(2 * self.rout_m, self.distance_m))

            self.mas_per_px = rad_to_mas(self.rad_per_px)

            self.resolution_mask_soft = rad_to_mas(mum_to_m(self.WAVELENGTH) / max(self.obs_data.blengths))
            self.resolution_mask_hard = rad_to_mas(1.22 * mum_to_m(self.WAVELENGTH) / max(self.obs_data.blengths))

    #
    # Functions
    #

    def observe(self, load_IQUV = True):
        ''' When we "look" at the star, we produce the interferometric measurements.'''
        if load_IQUV:
            self.__load_IQUV()
        self.__obtain_polarised_images()
        self.__obtain_power_spectrums()
        self.__obtain_PVRs()

        self.__calculate_rad_per_px()
        self.__calculate_max_baseline()
        self.__load_uv_coords()
        self.__sample_PVRs()

    def display_PVRs(self, overplot=True, scale="auto", interpolation=None, color_scale="auto"):
        ''' Displays the polarised visibility ratios '''

        # scale the colormap
        if color_scale == "auto":
            max_val = np.amax(np.concatenate((self.data_Q, self.data_U))) - 1
        else:
            max_val = color_scale

        # scales the image to zoom in on the baseline samples
        if scale == "off":
            ext = [-self.max_baseline, self.max_baseline, -self.max_baseline, self.max_baseline]
            PVR_Q_img = self.PVR_Q
            PVR_U_img = self.PVR_U
        else:
            # calculates how far out the maximum baseline of the mask will sit from the center of the image in pixels
            max_baseline_mask_px = int(
                np.ceil(max(self.obs_data.blengths) / self.max_baseline * float(self.parameters.grid_nx) * 0.5))
            margin_px = int(np.ceil(0.1 * max_baseline_mask_px))
            half_size_px = max_baseline_mask_px + margin_px
            center_px = np.ceil(float(self.parameters.grid_nx) / 2)
            c1 = int(center_px - half_size_px)  # first cut off (left and top cut off pixel)
            c2 = int(center_px + half_size_px)  # second cut off (right and bottom cut off pixel)
            baseline_per_px = 2 * self.max_baseline / float(
                self.parameters.grid_nx)  # meters per px that correspond to the baseline

            bl_extent = baseline_per_px * half_size_px  # with of the image in baseline units (meters)

            ext = [-bl_extent, bl_extent, -bl_extent, bl_extent]
            PVR_Q_img = self.PVR_Q[c1:c2, c1:c2]
            PVR_U_img = self.PVR_U[c1:c2, c1:c2]

        plot.figure(figsize=(12, 5))

        plot.subplot(1, 2, 1)
        plot.imshow(PVR_Q_img, cmap="magma", interpolation=interpolation, extent=ext)
        plot.title("Polarised visibility ratio (Q)")
        plot.xlabel("u (m)")
        plot.ylabel("v (m)")
        plot.colorbar()
        plot.clim(1 - max_val, 1 + max_val)

        if overplot == True:
            plot.plot(self.u_meters, self.v_meters, 'o', alpha=0.7, mew=0.0)
            plot.plot(self.u_neg_meters, self.v_neg_meters, 'o', alpha=0.7, mew=0.0)

        plot.subplot(1, 2, 2)
        plot.imshow(PVR_U_img, cmap="magma", interpolation=interpolation, extent=ext)
        plot.title("Polarised visibility ratio (U)")
        plot.xlabel("u (m)")
        plot.ylabel("v (m)")
        plot.colorbar()
        plot.clim(1 - max_val, 1 + max_val)

        if overplot == True:
            plot.plot(self.u_meters, self.v_meters, 'o', alpha=0.7, mew=0.0)
            plot.plot(self.u_neg_meters, self.v_neg_meters, 'o', alpha=0.7, mew=0.0)

    def display_data(self, option=6, ylims="auto", epsilon=0.001):
        ''' Displays the sampled values of the PVR in a 2D plot.
		Option 1: Simulated data
		Option 2: Observed data
		Option 3: Simulated and observed data in Stokes Q (separate plots)
		Option 4: Simulated and observed data in Stokes U (separate plots)
		Option 5: Simulated and observed data in Stokes Q and U (separate plots)
		Option 6: Simulated and observed data in Stokes Q and U (Q's and U's on the same plot)'''

        if ylims != False and ylims != "auto":
            ylims = ylims

        sim_vals = np.append(self.data_Q, self.data_U)
        if option != 1:
            obs_vals = np.append(self.obs_data.vhvv, self.obs_data.vhvvu)
            Q_vals = np.append(self.obs_data.vhvv, self.data_Q)
            U_vals = np.append(self.obs_data.vhvvu, self.data_U)
            all_vals = np.append(Q_vals, U_vals)

        if option == 1:

            if ylims == "auto":
                ylims = [min(sim_vals) - epsilon, max(sim_vals) + epsilon]

            fig = plot.figure(figsize=(15, 6))
            ax = fig.add_subplot(121)
            self.__display_sim_data(ax, Stokes="Q", ylims=ylims)
            ax = fig.add_subplot(122)
            self.__display_sim_data(ax, Stokes="U", ylims=ylims)

        if option == 2:

            if ylims == "auto":
                ylims = [min(obs_vals) - epsilon, max(obs_vals) + epsilon]

            fig = plot.figure(figsize=(15, 7))
            ax = fig.add_subplot(121)
            self.__display_observed_data(ax, Stokes="Q", ylims=ylims)
            ax = fig.add_subplot(122)
            self.__display_observed_data(ax, Stokes="U", ylims=ylims)

        if option == 3:

            if ylims == "auto":
                ylims = [min(Q_vals) - epsilon, max(Q_vals) + epsilon]

            fig = plot.figure(figsize=(15, 7))
            ax = fig.add_subplot(121)
            self.__display_observed_data(ax, Stokes="Q", ylims=ylims)
            ax = fig.add_subplot(122)
            self.__display_sim_data(ax, Stokes="Q", ylims=ylims)

        if option == 4:

            if ylims == "auto":
                ylims = [min(U_vals) - epsilon, max(U_vals) + epsilon]

            fig = plot.figure(figsize=(15, 7))
            ax = fig.add_subplot(121)
            self.__display_observed_data(ax, Stokes="U", ylims=ylims)
            ax = fig.add_subplot(122)
            self.__display_sim_data(ax, Stokes="U", ylims=ylims)

        if option == 5:

            if ylims == "auto":
                ylims = [min(all_vals) - epsilon, max(all_vals) + epsilon]

            fig = plot.figure(figsize=(15, 10))
            ax = fig.add_subplot(221)
            self.__display_observed_data(ax, Stokes="Q", ylims=ylims)
            ax = fig.add_subplot(222)
            self.__display_observed_data(ax, Stokes="U", ylims=ylims)
            ax = fig.add_subplot(223)
            self.__display_sim_data(ax, Stokes="Q", ylims=ylims)
            ax = fig.add_subplot(224)
            self.__display_sim_data(ax, Stokes="U", ylims=ylims)

        if option == 6:

            if ylims == "auto":
                ylims = [min(U_vals) - epsilon, max(U_vals) + epsilon]
            fig = plot.figure(figsize=(15, 7))
            ax = fig.add_subplot(121)
            self.__display_sim_data(ax, Stokes="Q", ylims=ylims, colorbar=False, marker='o')
            self.__display_observed_data(ax, Stokes="Q", ylims=ylims)
            ax.set_title("Stokes Q")
            ax.legend()
            ax = fig.add_subplot(122)
            self.__display_sim_data(ax, Stokes="U", ylims=ylims, colorbar=False, marker='o')
            self.__display_observed_data(ax, Stokes="U", ylims=ylims)
            ax.set_title("Stokes U")
            ax.legend()

    def print_parameters(self):
        for parameter in sorted(self.parameters.__dict__.keys()):
            print(parameter, ":\t", self.parameters.__dict__[parameter])

    #
    #	Testing Functions
    #

    def display_IQUV(self, scale="auto", interpolation=None):
        ''' Displays the images corresponding to the polarised images. These are the polarised intensities, NOT frational polarisation [I, Q/I, U/I, V/I] '''
        plot.figure(figsize=(9, 6))
        titles = ['I', 'Q', 'U', 'V']
        if scale == "off":
            half_size_mas = 0.5 * float(self.parameters.grid_nx) * self.mas_per_px
            images_IQUV = [self.I, self.Q, self.U, self.V]
        else:

            print('Auto scaling IQUV. To turn this off, use the option scale = "off"')
            # Calculate the radius of the outer shell in pixels
            rout_px = np.ceil(0.5 * self.dout_mas / self.mas_per_px)  # round up to a whole number of pixels
            if rout_px > 0.5 * int(self.parameters.grid_nx):
                print("Envelope is larger than image. Do not autoscale.")
            margin_px = np.ceil(0.1 * rout_px)  # add a 10 percent margin

            half_size_px = rout_px + margin_px
            center_px = np.ceil(float(self.parameters.grid_nx) / 2)
            c1 = int(center_px - half_size_px)  # first cut off (left and top cut off pixel)
            c2 = int(center_px + half_size_px)  # second cut off (right and bottom cut off pixel)
            half_size_mas = half_size_px * self.mas_per_px

            images_IQUV = [self.I[c1:c2, c1:c2], self.Q[c1:c2, c1:c2], self.U[c1:c2, c1:c2], self.V[c1:c2, c1:c2]]

        for i in range(4):
            plot.subplot(2, 2, i + 1)
            plot.subplots_adjust(hspace=0.4)
            plot.imshow(images_IQUV[i], cmap="magma",
                        extent=[-half_size_mas, half_size_mas, -half_size_mas, half_size_mas],
                        interpolation=interpolation)
            plot.xlabel("($mas$)")
            plot.ylabel("($mas$)")
            plot.title(titles[i])
            plot.colorbar().set_label("$W.m^2.pixel^{-1}$")

    def display_IQU(self, scale="auto", interpolation=None):
        ''' Displays the images corresponding to the polarised images. These are the polarised intensities, NOT frational polarisation [I, Q/I, U/I, V/I] '''
        plot.figure(figsize=(9, 3))
        titles = ['I', 'Q', 'U']
        if scale == "off":
            half_size_mas = 0.5 * float(self.parameters.grid_nx) * self.mas_per_px
            images_IQUV = [self.I, self.Q, self.U]
        else:

            print('Auto scaling IQU. To turn this off, use the option scale = "off"')
            # Calculate the radius of the outer shell in pixels
            rout_px = np.ceil(0.5 * self.dout_mas / self.mas_per_px)  # round up to a whole number of pixels
            if rout_px > 0.5 * int(self.parameters.grid_nx):
                print("Envelope is larger than image. Do not autoscale.")
            margin_px = np.ceil(0.1 * rout_px)  # add a 10 percent margin

            half_size_px = rout_px + margin_px
            center_px = np.ceil(float(self.parameters.grid_nx) / 2)
            c1 = int(center_px - half_size_px)  # first cut off (left and top cut off pixel)
            c2 = int(center_px + half_size_px)  # second cut off (right and bottom cut off pixel)
            half_size_mas = half_size_px * self.mas_per_px

            images_IQU = [self.I[c1:c2, c1:c2], self.Q[c1:c2, c1:c2], self.U[c1:c2, c1:c2]]

        for i in range(3):
            plot.subplot(1, 3, i + 1)
            plot.subplots_adjust(wspace = 0.8, hspace=0.4)
            plot.imshow(images_IQU[i], cmap="magma",
                        extent=[-half_size_mas, half_size_mas, -half_size_mas, half_size_mas],
                        interpolation=interpolation)
            plot.xlabel("($mas$)")
            plot.ylabel("($mas$)")
            plot.title(titles[i])
            plot.colorbar().set_label("$W.m^2.pixel^{-1}$")

    def display_power_spectrums(self):

        plot.figure(figsize=(9, 7))
        titles = ['H', 'V', 'A', 'B']
        images_HVAB = [self.Hps, self.Vps, self.Aps, self.Bps]

        for i in range(4):
            plot.subplot(2, 2, i + 1)
            plot.imshow(np.log10(images_HVAB[i]), cmap="magma", interpolation="nearest",
                        extent=[-self.max_baseline, self.max_baseline, -self.max_baseline, self.max_baseline])
            plot.title(titles[i] + " (log 10 scale)")
            plot.colorbar()
            plot.xlabel("u (m)")
            plot.ylabel("v (m)")

    def checksum(self):
        ''' For monochromatic coherent radiation: I^2 = Q^2 + U^2 + V^2. If this is a very small number I think it's ok: machine precision?'''

        checksum = np.power(self.I, 2) - np.power(self.Q, 2) - np.power(self.U, 2) - np.power(self.V, 2)
        print(np.amax(checksum))

    def calculate_reduced_chi2error(self):
        ''' Calculates the reduced chi square error in Stokes Q and Stokes U '''

        self.n_data_points = len(self.obs_data.vhvv)

        self.chi2err_Q = np.sum(pow((self.obs_data.vhvv - self.data_Q) / self.obs_data.vhvverr, 2))
        self.chi2err_U = np.sum(pow((self.obs_data.vhvvu - self.data_U) / self.obs_data.vhvvuerr, 2))

        self.reduced_chi2err_Q = self.chi2err_Q / (self.n_data_points - self.n_free_parameters)
        self.reduced_chi2err_U = self.chi2err_U / (self.n_data_points - self.n_free_parameters)

        self.reduced_chi2err = 0.5 * (self.reduced_chi2err_Q + self.reduced_chi2err_U)

    def reduce(self):
        '''This is to reduce the object size so that it can be saved'''

        del self.data

        del self.Hp
        del self.Vp
        del self.Ap
        del self.Bp

        del self.Hps
        del self.Vps
        del self.Aps
        del self.Bps

        del self.PVR_Q
        del self.PVR_U


    def sanity_check_star(self, star_edge="hard", shell=False, prec=5):
        ''' Does a quick sanity check.
		INSTRUCTIONS:
			- Make the dust shell very faint by reducing the dust mass. E.g. 'dust_mass': 1e-20.
				(note: you cannot set the dust mass to 0)
			- Typically you want the star to have a hard edge. This means the formula for the resolution that we use is:
					resolution = 1.220 lambda / baseline
				and is the default setting for star_edge.
			- If the star has a soft edge, set star_edge = "soft", in which case the resolution is approximately:
					resolution = lambda / baseline
			- If you want to see calculations for the size of the shell, set shell = True.
			- Example parameter set:
				free_params = {"Rin": 5, "Rout": 30, "dust_mass": 1e-20, "radius": 1000, "size": 80}
			'''

        # Running message
        print("--- SANITY CHECK (star) ---")
        print()

        if star_edge.lower() == "hard":
            Rayleigh_criterion_const = 1.220
        elif star_edge.lower() == "soft":
            Rayleight_criterion_const = 1
        else:
            print("Specify star_edge as either 'hard' or 'soft'.")
            print("Assuming the star has a hard edge...")
            Rayleigh_criterion_const = 1.22

        print("Distance to the star (pc):\t\t{}".format(self.parameters.distance))
        print("Radius of the star (R_sun):\t\t{}".format(self.parameters.radius))
        if shell:
            print("Inner radius of dust shell (AU):\t{}".format(self.parameters.Rin))
            print("Outer radius of dust shell (AU):\t{}".format(self.parameters.Rout))
        print()

        print("Distance to the star (m): \t\t{:.{prec}}".format(self.distance_m, prec=prec))
        print("Radius of the star (m): \t\t{:.{prec}}".format(self.radius_m, prec=prec))
        if shell:
            print("Inner radius of dust shell (m):\t\t{:.{prec}}".format(self.rin_m, prec=prec))
            print("Outer radius of dust shell (m):\t\t{:.{prec}}".format(self.rout_m, prec=prec))
        print()

        print("Radius (diameter) of the star (mas):\t\t\t{:.{prec}}\t\t({:.{prec}})".format(0.5 * self.diameter_mas,
                                                                                            self.diameter_mas,
                                                                                            prec=prec))
        if shell:
            print("Inner radius (diameter) of the dust shell (mas):\t{:.{prec}}\t\t({:.{prec}})".format(
                0.5 * self.din_mas, self.din_mas, prec=prec))
            print("Outer radius (diamter) of the dust shell (mas):\t\t{:.{prec}}\t\t({:.{prec}})".format(
                0.5 * self.dout_mas, self.dout_mas, prec=prec))
        print()

        self.__sample_power_spectrum_I()
        self.__display_samples_power_spectrum_I()

        plot.ion()
        plot.show()

        print("Look for the baseline length where the visibility goes to 0.")
        min_zero_baseline = float(input("Give an UNDERESTIMATE for this baseline length (in meters): "))
        max_zero_baseline = float(input("Give an OVERESTIMATE for this baseline length (in meters): "))
        print()

        plot.close()
        plot.ioff()

        print("Wavelength of observation (microns): ", self.WAVELENGTH)
        max_resolution = Rayleigh_criterion_const * mum_to_m(self.WAVELENGTH) / min_zero_baseline
        min_resolution = Rayleigh_criterion_const * mum_to_m(self.WAVELENGTH) / max_zero_baseline
        print("The baseline at which the visibility is zero corresponds to a resolution range (mas):")
        print("{:.{prec}}  -  {:.{prec}}".format(rad_to_mas(min_resolution), rad_to_mas(max_resolution), prec=prec))
        print("The diameter of the star should lie in this resolution range. ")

    def sanity_check_shell(self, star_edge="hard", prec=5):
        ''' TIP: if the simulated data looks messy, try making the star and the shell smaller'''

        print("--- SANITY CHECK (shell) ---")
        print()

        if star_edge.lower() == "hard":
            Rayleigh_criterion_const = 1.220
        elif star_edge.lower() == "soft":
            Rayleigh_criterion_const = 1
        else:
            print("Specify star_edge as either 'hard' or 'soft'.")
            print("Assuming the star has a hard edge...")
            Rayleigh_criterion_const = 1.220

        print(
            "The maximum baseline for this mask is (meters)\t{:.{prec}}".format(max(self.obs_data.blengths), prec=prec))
        resolution = Rayleigh_criterion_const * mum_to_m(self.WAVELENGTH) / max(self.obs_data.blengths)
        print("The resolution of this mask is (mas):\t\t{:.{prec}}".format(rad_to_mas(resolution), prec=prec))
        print("Distance to the star (pc):\t\t\t{:.{prec}}".format(self.parameters.distance, prec=prec))
        print()

        # Calculate the diameter of a star just resolved with shell twice the size of the star
        radius_star_m = np.tan(0.5 * resolution) * self.distance_m
        radius_star_sr = AU_to_sr(m_to_AU(radius_star_m))
        print("A star that is just resolved has radius (solar radii):\t\t{:.{prec}}".format(radius_star_sr, prec=prec))
        print(
            "A shell that is twice the size of the star has radius (AU):\t{:.{prec}}".format(2 * m_to_AU(radius_star_m),
                                                                                             prec=prec))
        print()

        print("--- Make a star that is just resolved with a shell twice the size of the star ---")

        print("Radius of the star (R_sun):\t\t{}".format(self.parameters.radius))
        print("Inner radius of dust shell (AU):\t{}".format(self.parameters.Rin))
        print("Outer radius of dust shell (AU):\t{}".format(self.parameters.Rout))
        print()

        print("Distance to the star (m):\t{:.{prec}}".format(self.distance_m, prec=prec))
        print("Radius of the star (m):\t\t{:.{prec}}".format(self.radius_m, prec=prec))
        print("Inner radius of dust shell (m):\t{:.{prec}}".format(self.rin_m, prec=prec))
        print("Outer radius of dust shell (m):\t{:.{prec}}".format(self.rout_m, prec=prec))
        print()

        print("Radius (diameter) of the star (mas):\t\t\t{:.{prec}}\t({:.{prec}})".format(0.5 * self.diameter_mas,
                                                                                          self.diameter_mas, prec=prec))
        print("Inner radius (diameter) of the dust shell (mas):\t{:.{prec}}\t({:.{prec}})".format(0.5 * self.din_mas,
                                                                                                  self.din_mas,
                                                                                                  prec=prec))
        print("Outer radius (diamter) of the dust shell (mas):\t\t{:.{prec}}\t({:.{prec}})".format(0.5 * self.dout_mas,
                                                                                                   self.dout_mas,
                                                                                                   prec=prec))
        print()

        self.display_data(option=1)
        plot.ion()
        plot.show()

        print("You should see distinct sinusoidal shapes. If you do not, you may need to zoom in.")

        print("Look at the blue-green sinusoid where the amplitude is largest. The colour tells you baseline.")
        min_baseline = float(input("What is an UNDERESTIMATE for the baseline length of this sinusoid? "))
        max_baseline = float(input("What is an OVERESTIMATE for the baseline length of this sinusoid? "))
        print()

        plot.close()
        plot.ioff()

        print("Wavelength of observation (microns): {}".format(self.WAVELENGTH))
        max_size_shell = Rayleigh_criterion_const * mum_to_m(self.WAVELENGTH) / min_baseline
        min_size_shell = Rayleigh_criterion_const * mum_to_m(self.WAVELENGTH) / max_baseline

        print("The diameter of the SHELL should lie between (mas):")
        print("{:.{prec}}  -  {:.{prec}}".format(rad_to_mas(min_size_shell), rad_to_mas(max_size_shell), prec=prec))

    #
    #	Private functions
    #

    def __load_IQUV(self):
        ''' Unloads the fits file into the I, Q U and V images'''
        self.I = self.data[0][0][0]
        self.Q = self.data[1][0][0]
        self.U = self.data[2][0][0]
        self.V = self.data[3][0][0]

    def __obtain_polarised_images(self):
        ''' Obtain the polarised images'''
        self.Hp = self.I + self.Q  # horizontal
        self.Vp = self.I - self.Q  # vertical
        self.Ap = self.I + self.U
        self.Bp = self.I - self.U

    def __obtain_power_spectrums(self):
        ''' These power spectrums are the absoulate value of the square of the 2D fourier transformation of the original images. The have also been normalised and re-centered.'''
        self.Hps = normalise(center(power_spectrum(self.Hp)))
        self.Vps = normalise(center(power_spectrum(self.Vp)))
        self.Aps = normalise(center(power_spectrum(self.Ap)))
        self.Bps = normalise(center(power_spectrum(self.Bp)))

    def __obtain_PVRs(self):
        self.PVR_Q = self.Hps / self.Vps
        self.PVR_U = self.Aps / self.Bps

    def __calculate_rad_per_px(self):
        self.distance_AU = pc_to_AU(float(self.parameters.distance))
        self.angular_size_image_rad = angular_size(self.parameters.size, self.distance_AU)
        self.rad_per_px = self.angular_size_image_rad / float(self.parameters.grid_nx)

    def __calculate_max_baseline(self):
        ''' The resolution of a telescope is: resolution = wavelength / baseline length.
		Hence, the rad_per_px will correspond to a maximum baseline length the power spectrum can represent.
		This is what we will calculate. '''

        # Each pixel corresponds to a number of radians.
        # The smallest thing we can see, is 2 pixels across, so we need to double the radians per pixel

        self.resolution = 2 * self.rad_per_px

        # calculate the corresponding maximum baseline: baseline length = wavelength / resolution
        self.max_baseline = mum_to_m(float(self.WAVELENGTH)) / self.resolution

    def __load_uv_coords(self):

        ''' Returns the co-ordinates of where we should sample the polarised visibility ratios in units of meters.'''

        self.u_meters = [u * mum_to_m(self.WAVELENGTH) for u in self.obs_data.u_coords]
        self.v_meters = [v * mum_to_m(self.WAVELENGTH) for v in self.obs_data.v_coords]
        self.u_neg_meters = [-u for u in self.u_meters]
        self.v_neg_meters = [-v for v in self.v_meters]

        # Check that this is consistent with the baseline lengths
        inconsistent_baselines = []
        for i in range(len(self.obs_data.blengths)):

            uv_length = np.sqrt(pow(self.u_meters[i], 2) + pow(self.v_meters[i], 2))
            if round(uv_length, 4) != round(self.obs_data.blengths[i], 4):
                inconsistent_baselines.append((uv_length, self.obs_data.blengths[i]))

        # Print out a warning if the baseline lengths are inconsistent
        if len(inconsistent_baselines) != 0:
            print(
                "...\nWARNING: The baseline lengths calculated from the u,v co-ordinates do not match the stored baseline lengths.\nIt is likely there is a small rounding difference, or the wavelength is wrong.\nThe values that do not match are:")
            for i in inconsistent_baselines:
                print(i)
            print("...")

    def __sample_PVR(self, PVR):
        ''' To be used inside sample PVRs. A function that samples a given PVR, and returns the sampled values.'''

        f = interpolate.interp2d(self.U_vals, self.V_vals, PVR, kind="linear")
        sampled_PVR = []
        for i in range(len(self.u_meters)):
            sampled_PVR.append(f(self.u_meters[i], self.v_meters[i])[0])

        return np.array(sampled_PVR)

    def __sample_PVRs(self):
        ''' Samples the polarised visibility ratio of both Q and U'''

        # We assume the u and v are symmetric
        self.U_vals = np.linspace(-1, 1, self.parameters.grid_nx) * self.max_baseline
        self.V_vals = self.U_vals

        self.data_Q = self.__sample_PVR(self.PVR_Q)
        self.data_U = self.__sample_PVR(self.PVR_U)

    def __display_sim_data(self, ax, max_baseline=8, Stokes='Q', ylims=False, colorbar=True, marker='x'):
        ''' Displays the simulated data produced from the model. Specify whether you want the data from Stokes Q or Stokes U. The maximum baseline is the largest baseline that will be displayed.'''

        baselines = self.obs_data.blengths
        # only select the baselines that are smaller than the maximum baseline
        baseline_colours = baselines[baselines <= max_baseline]

        if Stokes == 'Q':
            data = self.data_Q
            title = "Simulated Data (Stokes Q)"
        else:
            Stokes = 'U'
            data = self.data_U
            title = "Simulated Data (Stokes U)"

        scatter_plot = ax.scatter(self.obs_data.bazims, data, s=25, marker=marker, label="simulated",
                                  c=baseline_colours)
        if colorbar:
            clb = plot.colorbar(scatter_plot, label="Baseline length (m)")

        if ylims:
            ax.set_ylim(ylims)

        ax.set_title(title)
        ax.set_xlabel("Baseline Azimuth Angle (rads)")
        ax.set_ylabel("Polarised Visibility Ratio")

    def __display_observed_data(self, ax, max_baseline=8, Stokes='Q', ylims=False, marker='x'):
        ''' Displays the data from the observation. Specify whether you want the data from Stokes Q or Stokes U. The maximum baseline is the largest baseline that will be displayed.'''

        baselines = self.obs_data.blengths
        # only select the baselines that are smaller than the maximum baseline
        baseline_colours = baselines[baselines <= max_baseline]

        if Stokes == 'Q':
            data = self.obs_data.vhvv
            errors = self.obs_data.vhvverr
            title = "Observed Data (Stokes Q)"
        else:
            Stokes = 'U'
            data = self.obs_data.vhvvu
            errors = self.obs_data.vhvvuerr
            title = "Observed Data (Stokes U)"

        # Plots the markers
        scatter_plot = ax.scatter(self.obs_data.bazims, data, s=25, marker=marker, label="observed", c=baseline_colours)
        clb = plot.colorbar(scatter_plot, label="Baseline length (m)")
        bar_colour = clb.to_rgba(baseline_colours)

        if ylims:
            ax.set_ylim(ylims)

        # Plots the errorbars
        ax.errorbar(self.obs_data.bazims, data, yerr=errors, marker='', linestyle='', alpha=0.8, capsize=0, zorder=0,
                    ecolor=bar_colour)

        ax.set_title(title)
        ax.set_xlabel("Baseline Azimuth Angle (rads)")
        ax.set_ylabel("Polarised Visibility Ratio")

    def __sample_power_spectrum_I(self):
        self.Ips = normalise(center(power_spectrum(self.I)))
        self.data_I = self.__sample_PVR(self.Ips)

    def __display_samples_power_spectrum_I(self):
        ''' The first plot is of the powerspectrum of I with the location of the samples over the image. The second is the samples sorted by baseline length. '''

        fig = plot.figure(figsize=(15, 7))

        ax = fig.add_subplot(121)
        plot.imshow(np.log10(self.Ips), cmap="magma", interpolation="nearest",
                    extent=[-self.max_baseline, self.max_baseline, -self.max_baseline, self.max_baseline],
                    shape=[10, 10])
        plot.title("Intentisty Power Spectrum (log 10 scale)")
        plot.xlabel("u (m)")
        plot.ylabel("v (m)")

        # over plotting
        plot.plot(self.u_meters, self.v_meters, 'o')
        plot.plot(self.u_neg_meters, self.v_neg_meters, 'o')

        ax = fig.add_subplot(122)
        plot.plot(self.obs_data.blengths, self.data_I, 'x')
        plot.title("Samples of the power spectrum of I")
        plot.xlabel("Baseline length (m)")
        plot.ylabel("Visibility I")
