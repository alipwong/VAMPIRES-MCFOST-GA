def write_MCFOST_parafile(parameters, MCFOST_path, file_name = "star.para"):

	f = open(MCFOST_path + file_name, 'w')

	f.write("3.0                      mcfost version\n\n")

	f.write("#Number of photon packages\n")
	f.write("  " + parameters.nbr_photons_eq_th + "                  nbr_photons_eq_th  : T computation\n")
	f.write("  " + parameters.nbr_photons_lambda + "	          nbr_photons_lambda : SED computation\n")
	f.write("  " + parameters.nbr_photons_image + "                  nbr_photons_image  : images computation\n\n")

	f.write("#Wavelength\n")
	f.write("  " + parameters.n_lambda + "  " + parameters.lambda_min + " " + parameters.lambda_max + "          n_lambda, lambda_min, lambda_max [mum] Do not change this line unless you know what you are doing\n")
	f.write("  " + parameters.compute_temp + " " + parameters.compute_sed + " " + parameters.default_wavelength + " 		  compute temperature?, compute sed?, use parameters wavelength grid for output ?\n")
	f.write("  " + parameters.wavelength_file + "		  wavelength file (if previous parameter is F)\n")
	f.write("  " + parameters.separation_of_different_contributions + " " + parameters.stokes_parameters + "			  separation of different contributions?, stokes parameters?\n\n")

	f.write("#Grid geometry and size\n")
	f.write("  " + parameters.geometry + "			  1 = cylindrical, 2 = spherical, 3 = Voronoi tessellation (this is in beta, please ask Christophe)\n")
	f.write("  " + parameters.n_rad + " " + parameters.nz + " " + parameters.n_az + " " + parameters.n_rad_in + "             n_rad (log distribution), nz (or n_theta), n_az, n_rad_in\n\n")

	f.write("#Maps\n")
	f.write("  " + parameters.grid_nx + " " + parameters.grid_ny + " " + parameters.size + "            grid (nx,ny), size [AU]\n")
	f.write("  " + parameters.imin + "  " + parameters.imax + "  " + parameters.n_incl + "  " + parameters.centered + "          RT: imin, imax, n_incl, centered ?\n")
	f.write("  " + parameters.az_min + "    " + parameters.az_max + "   " + parameters.n_az_angles + "             RT: az_min, az_max, n_az angles\n")
	f.write("  " + parameters.distance + "			  distance (pc)\n")
	f.write("  " + parameters.disk_PA + "			  disk PA\n\n")

	f.write("#Scattering method\n")
	f.write("  " + parameters.scattering_method + "	                  0=auto, 1=grain prop, 2=cell prop\n")
	f.write("  " + parameters.mie_hg + "	                  1=Mie, 2=hg (2 implies the loss of polarizarion)\n\n")

	f.write("#Symetries\n")
	f.write("  " + parameters.image_symmetry + "	                  image symmetry\n")
	f.write("  " + parameters.central_symmetry + "	                  central symmetry\n")
	f.write("  " + parameters.axial_symmetry + "	                  axial symmetry (important only if N_phi > 1)\n\n")

	f.write("#Disk physics\n")
	f.write("  " + parameters.dust_settling  + "     " + parameters.exp_strat + "  " + parameters.a_strat + "	  dust_settling (0=no settling, 1=parametric, 2=Dubrulle, 3=Fromang), exp_strat, a_strat (for parametric settling)\n")
	f.write("  " + parameters.dust_radial_migration + "                       dust radial migration\n")
	f.write("  " + parameters.sublimate_dust + "		  	  sublimate dust\n")
	f.write("  " + parameters.hydrostatic_equilibrium + "                       hydrostatic equilibrium\n")
	f.write("  " + parameters.viscous_heating + "  " + parameters.alpha_viscosity + "		  viscous heating, alpha_viscosity\n\n")

	f.write("#Number of zones : 1 zone = 1 density structure + corresponding grain properties\n")
	f.write("  " + parameters.number_of_zones + "\n\n")

	f.write("#Density structure\n")
	f.write("  " + parameters.zone_type + "                       zone type : 1 = disk, 2 = tappered-edge disk, 3 = envelope, 4 = debris disk, 5 = wall\n")
	f.write("  " + parameters.dust_mass + "    " + parameters.gas_to_dust_mass_ratio + "		  dust mass,  gas-to-dust mass ratio\n")
	f.write("  " + parameters.scale_height + "  " + parameters.reference_radius + "  " + parameters.vertical_profile_exponent + "           scale height, reference radius (AU), unused for envelope, vertical profile exponent (only for debris disk)\n")
	f.write("  " + parameters.Rin + "  " + parameters.edge + "    " + parameters.Rout + "  " + parameters.Rc + "  Rin, edge, Rout, Rc (AU) Rc is only used for tappered-edge & debris disks (Rout set to 8*Rc if Rout==0)\n")
	f.write("  " + parameters.flaring_exponent + "                   flaring exponent, unused for envelope\n")
	f.write("  " + parameters.surface_density_exponent + "  " + parameters.negative_gamma_exp + "    	          surface density exponent (or -gamma for tappered-edge disk or volume density for envelope), usually < 0, -gamma_exp (or alpha_in & alpha_out for debris disk)\n\n")

	f.write("#Grain properties\n")
	f.write("  " + parameters.number_of_species + "  Number of species\n")
	f.write("  " + parameters.grain_type + "  " + parameters.N_components + " " + parameters.mixing_rule + "  " + parameters.porosity + "  " + parameters.max_fraction + "  " + parameters.Vmax + " Grain type (Mie or DHS), N_components, mixing rule (1 = EMT or 2 = coating),  porosity, mass fraction, Vmax (for DHS)\n")
	f.write("  " + parameters.optical_indicies_file + "  " + parameters.volume_fraction + "  Optical indices file, volume fraction\n")
	f.write("  " + parameters.heating_method + "	                  Heating method : 1 = RE + LTE, 2 = RE + NLTE, 3 = NRE\n")
	f.write("  " + parameters.amin + "  " + parameters.amax + " " + parameters.aexp + " " + parameters.n_grains + " 	  amin, amax [mum], aexp, n_grains (log distribution)\n\n")

	f.write("#Molecular RT parameters\n")
	f.write("  " + parameters.lpop + " " + parameters.laccurate_pop + " " + parameters.LTE + " " + parameters.profile_width + "	          lpop, laccurate_pop, LTE, profile width (km.s^-1)\n")
	f.write("  " + parameters.v_turb + " 			  v_turb (delta)\n")
	f.write("  " + parameters.nmol + "			  nmol\n")
	f.write("  " + parameters.molecular_data_filename + " " + parameters.level_max + "           molecular data filename, level_max\n")
	f.write("  " + parameters.vmax + " " + parameters.n_speed + "     	  	  vmax (km.s^-1), n_speed\n")
	f.write("  " + parameters.cst_molecule_abundance + " " + parameters.abundance + " " + parameters.abundance_file + "   cst molecule abundance ?, abundance, abundance file\n")
	f.write("  " + parameters.ray_tracing + "  " + parameters.number_lines_in_RT + "                       ray tracing ?,  number of lines in ray-tracing\n")
	f.write("  " + parameters.transition_number_1 + " " + parameters.transition_number_2 + " " + parameters.transition_number_3 + "	 		  transition numbers\n\n")

	f.write("#Star properties\n")
	f.write("  " + parameters.number_of_stars + " Number of stars\n")
	f.write("  " + parameters.temp + "	" + parameters.radius + "	" + parameters.mass + "	" + parameters.x + "	" + parameters.y + "	" + parameters.z + "  " + parameters.is_blackbody + " Temp, radius (solar radius),M (solar mass),x,y,z (AU), is a blackbody?\n")
	f.write("  lte4000-3.5.NextGen.fits.gz\n")

	f.write("  " + parameters.fUV + "	" + parameters.slope_FUV + "  fUV, slope_fUV\n")
	f.close()