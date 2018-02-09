#Code for Astro301 Homework 1
#Astro301 Spring 2018
#Author: Chase Urasaki
#Created: 2.3.2018
#Updated: 2.7.2018

###################################################################

#Plot the four SED's onto one figure

#What we want to first load the necessary libraries

from matplotlib import pyplot
import numpy as np
from astropy.io import ascii
import re
from  mpl_toolkits.mplot3d import Axes3D
#Next, we need to read in the data 

elliptical = np.genfromtxt('CWW_E_ext.sed', dtype = float)
irregular = np.genfromtxt('CWW_Im_ext.sed', dtype = float)
Sbc = np.genfromtxt('CWW_Sbc_ext.sed', dtype = float)
Scd = np.genfromtxt('CWW_Scd_ext.sed', dtype = float)

#Pyplot is going to plot as x vs y. Thus, we need to differentiate the values from each of our information

elliptical_x = elliptical [:,0]
elliptical_y = elliptical [:,1]
irregular_x = irregular [:,0]
irregular_y = irregular [:,1]
Sbc_x = Sbc [:,0]
Sbc_y = Sbc [:,1]
Scd_x = Scd [:,0]
Scd_y = Scd [:,1]


#Now that the data is read in, time to start plotting
pyplot.plot(elliptical_x, elliptical_y, label = 'Elliptical Galaxy')
pyplot.plot(irregular_x, irregular_y, label = 'Irregular Galaxy')
pyplot.plot(Sbc_x, Sbc_y, label = 'Sbc-type Galaxy')
pyplot.plot(Scd_x, Scd_y, label = 'Scd-type Galaxy')

#Now, for the finishing touches for the plot
pyplot.xlim(0, 20000)
pyplot.xlabel('Wavelength (Angstroms)')
pyplot.ylabel('Relative Flux(unitless)')
pyplot.title('SED Spectra for Four Different Types of Galaxies')
pyplot.legend()
pyplot.show(); pyplot.close()

######################################################################

#Plot the filter curves

#Since all of the libraries are already imported, we don't have to worry about doing that. So we start off with reading the data in

Sloang = np.genfromtxt('Sloan_g.txt', dtype = float, skip_header = 1)
#g' = ECW - 475.4 nm; FWHM - 138.7 nm 
Sloani = np.genfromtxt('Sloan_i.txt', dtype = float, skip_header = 1)
#i' = ECW - 769.8 nm; FWHM - 130.0 nm
Sloanr = np.genfromtxt('Sloan_r.txt', dtype = float, skip_header = 1)
#r' = ECW - 620.4 nm; FWHM - 124.0 nm
Sloanu = np.genfromtxt('Sloan_u.txt', dtype = float, skip_header = 1)
#u' = ECW - 358.0 nm; FWHM - 33.9 nm
Sloanz = np.genfromtxt('Sloan_z.txt', dtype = float, skip_header = 1)
#z' = ECW - 966.5 nm; FWHM - 225.8 nm

#Now that we have read in all the data, we have to create the x and y values for each. The -1 says that we just revesed the order of the entries. 

Sloang_x = Sloang [:, 0]
Sloang_y = Sloang [:, 1]
Sloani_x = Sloani [:, 0]
Sloani_y = Sloani [:, 1]
Sloanr_x = Sloanr [:, 0]
Sloanr_y = Sloanr [:, 1]
Sloanu_x = Sloanu [:, 0]
Sloanu_y = Sloanu [:, 1]
Sloanz_x = Sloanz [:, 0]
Sloanz_y = Sloanz [:, 1]


#Reversing the order of the entries in the Sloan Stuff, 

Sloang_x = np.flip(Sloang_x, 0)
Sloang_y = np.flip(Sloang_y, 0)
Sloani_x = np.flip(Sloani_x, 0)
Sloani_y = np.flip(Sloani_y, 0)
Sloanr_x = np.flip(Sloanr_x, 0)
Sloanr_y = np.flip(Sloanr_y, 0)
Sloanu_x = np.flip(Sloanu_x, 0)
Sloanu_y = np.flip(Sloanu_y, 0)
Sloanz_x = np.flip(Sloanz_x, 0)
Sloanz_y = np.flip(Sloanz_y, 0)

#So that we can have an easier time later, we should normalize the curve now. Recall that we only need to normalize the "y - values"
#Note that this procedure normalizes for the respective max values. 

Sloang_y_norm = Sloang_y/max(Sloang_y)
Sloani_y_norm = Sloani_y/max(Sloani_y)
Sloanr_y_norm = Sloanr_y/max(Sloanr_y)
Sloanu_y_norm = Sloanu_y/max(Sloanu_y)
Sloanz_y_norm = Sloanz_y/max(Sloanz_y)

#We want to convert the x values in the Filters from nm to Angstroms
Sloang_x = np.array(Sloang_x, dtype = float)*10.
Sloani_x = np.array(Sloani_x, dtype = float)*10.
Sloanr_x = np.array(Sloanr_x, dtype = float)*10.
Sloanu_x = np.array(Sloanu_x, dtype = float)*10.
Sloanz_x = np.array(Sloanz_x, dtype = float)*10.

#Plot the above values 
pyplot.plot(Sloang_x, Sloang_y_norm, label = 'Sloan g')
pyplot.plot(Sloani_x, Sloani_y_norm, label = 'Sloan i')
pyplot.plot(Sloanr_x, Sloanr_y_norm, label = 'Sloan r')
pyplot.plot(Sloanu_x, Sloanu_y_norm, label = 'Sloan u')
pyplot.plot(Sloanz_x, Sloanz_y_norm, label = 'Sloan z')

#Add the finishing touches to the plots
pyplot.xlabel('Wavelength (Angstroms)')
pyplot.ylabel('Transmission  % (Normalized)')
pyplot.title('Filter Curves')
pyplot.legend()
pyplot.show();pyplot.close()

##########################################################################

#Now, we want to compute the galaxies color (assuming redshift 0)

#Now, we need to interpolate the filter data. Note that the interpolation will adjust the dimension of the array to suit everything :)

#This returns the y-values after you've interpolated the Sloan g transmittance with the wavelengths in the SED

#Compute the colors for the elliptical galaxy

Sloang_yinter_el = np.interp(elliptical_x, Sloang_x, Sloang_y_norm)
Sloani_yinter_el = np.interp(elliptical_x, Sloani_x, Sloani_y_norm)
Sloanr_yinter_el = np.interp(elliptical_x, Sloanr_x, Sloanr_y_norm)
Sloanu_yinter_el = np.interp(elliptical_x, Sloanu_x, Sloanu_y_norm)
Sloanz_yinter_el = np.interp(elliptical_x, Sloanz_x, Sloanz_y_norm)

#Multiply the values of the SED for the galaxy by the flux values

elliptical_flux_sloani = np.multiply(elliptical_y, Sloani_yinter_el)
elliptical_flux_sloang = np.multiply(elliptical_y, Sloang_yinter_el)
elliptical_flux_sloanr = np.multiply(elliptical_y, Sloanr_yinter_el)
elliptical_flux_sloanu = np.multiply(elliptical_y, Sloanu_yinter_el)
elliptical_flux_sloanz = np.multiply(elliptical_y, Sloanz_yinter_el)

#Now we are plotting the transmitted fluxes for the elliptical galaxy in each filter

pyplot.xlim(0, 20000)
pyplot.title('Elliptical Galaxy flux in Sloan Filters')
pyplot.xlabel('Wavelength (Angstroms)')
pyplot.ylabel('Relative Flux (unitless)')
pyplot.plot(elliptical_x, elliptical_flux_sloang, label = 'Sloan g')
pyplot.plot(elliptical_x, elliptical_flux_sloani, label = 'Sloan i')
pyplot.plot(elliptical_x, elliptical_flux_sloanr, label = 'Sloan r')
pyplot.plot(elliptical_x, elliptical_flux_sloanu, label = 'Sloan u')
pyplot.plot(elliptical_x, elliptical_flux_sloanz, label = 'Sloan z')
pyplot.legend()
pyplot.show();pyplot.close()

#Now that we have generated the plots, we can now find the flux in each filter, which will give us the magnitude in each filter, and from there, find the colors

elliptical_flux_g = np.trapz(elliptical_flux_sloang, elliptical_x)
elliptical_flux_i = np.trapz(elliptical_flux_sloani, elliptical_x)
elliptical_flux_r = np.trapz(elliptical_flux_sloanr, elliptical_x)
elliptical_flux_u = np.trapz(elliptical_flux_sloanu, elliptical_x)
elliptical_flux_z = np.trapz(elliptical_flux_sloanz, elliptical_x)

#Here is the magnitudes in each of the filters. 

El_g = -2.5*np.log10(elliptical_flux_g)
El_i = -2.5*np.log10(elliptical_flux_i)
El_r = -2.5*np.log10(elliptical_flux_r)
El_u = -2.5*np.log10(elliptical_flux_u)
El_z = -2.5*np.log10(elliptical_flux_z)

#Compute the colors

elliptical_u_g = (El_u) - (El_g)
elliptical_g_r = (El_g) - (El_r)
elliptical_r_i = (El_r) - (El_i)
elliptical_i_z = (El_i) - (El_z)

print('Colors for Elliptical Galaxy: u-g = {0:8.2f}, g-r = {1:8.2f}, r-i = {2:8.2f}, i-z = {3:8.2f}'.format(elliptical_u_g, elliptical_g_r, elliptical_r_i, elliptical_i_z))

#Calculate the colors of the irregular galaxy

Sloang_yinter_ir = np.interp(irregular_x, Sloang_x, Sloang_y_norm)
Sloani_yinter_ir = np.interp(irregular_x, Sloani_x, Sloani_y_norm)
Sloanr_yinter_ir = np.interp(irregular_x, Sloanr_x, Sloanr_y_norm)
Sloanu_yinter_ir = np.interp(irregular_x, Sloanu_x, Sloanu_y_norm)
Sloanz_yinter_ir = np.interp(irregular_x, Sloanz_x, Sloanz_y_norm)

#Multiply the values of the SED for the galaxy by the flux values

irregular_flux_sloani = np.multiply(irregular_y, Sloani_yinter_ir)
irregular_flux_sloang = np.multiply(irregular_y, Sloang_yinter_ir)
irregular_flux_sloanr = np.multiply(irregular_y, Sloanr_yinter_ir)
irregular_flux_sloanu = np.multiply(irregular_y, Sloanu_yinter_ir)
irregular_flux_sloanz = np.multiply(irregular_y, Sloanz_yinter_ir)

#Now we are plotting the transmitted fluxes for the elliptical galaxy in each filter

pyplot.xlim(0, 20000)
pyplot.title('Irregular Galaxy flux in Sloan Filters')
pyplot.xlabel('Wavelength (Angstroms)')
pyplot.ylabel('Relative Flux (Unitless)')
pyplot.plot(irregular_x, irregular_flux_sloang, label = 'Sloan g')
pyplot.plot(irregular_x, irregular_flux_sloani, label = 'Sloan i')
pyplot.plot(irregular_x, irregular_flux_sloanr, label = 'Sloan r')
pyplot.plot(irregular_x, irregular_flux_sloanu, label = 'Sloan u')
pyplot.plot(irregular_x, irregular_flux_sloanz, label = 'Sloan z')
pyplot.legend()
pyplot.show();pyplot.close()

#Now that we have generated the plots, we can now find the flux in each filter, which will give us the magnitude in each filter, and from there, find the colors

irregular_flux_g = np.trapz(irregular_flux_sloang, irregular_x)
irregular_flux_i = np.trapz(irregular_flux_sloani, irregular_x)
irregular_flux_r = np.trapz(irregular_flux_sloanr, irregular_x)
irregular_flux_u = np.trapz(irregular_flux_sloanu, irregular_x)
irregular_flux_z = np.trapz(irregular_flux_sloanz, irregular_x)

#Here is the magnitudes in each of the filters. 
Ir_g = -2.5*np.log10(irregular_flux_g)
Ir_i = -2.5*np.log10(irregular_flux_i)
Ir_r = -2.5*np.log10(irregular_flux_r)
Ir_u = -2.5*np.log10(irregular_flux_u)
Ir_z = -2.5*np.log10(irregular_flux_z)

#Compute the colors
irregular_u_g = (Ir_u) - (Ir_g)
irregular_g_r = (Ir_g) - (Ir_r)
irregular_r_i = (Ir_r) - (Ir_i)
irregular_i_z = (Ir_i) - (Ir_z)

print('Colors for Irregular Galaxy: u-g = {0:8.2f}, g-r = {1:8.2f}, r-i = {2:8.2f}, i-z = {3:8.2f}'.format(irregular_u_g, irregular_g_r, irregular_r_i, irregular_i_z))

#Compute the colors for the Sbc galaxy

Sloang_yinter_Sbc = np.interp(Sbc_x, Sloang_x, Sloang_y_norm)
Sloani_yinter_Sbc = np.interp(Sbc_x, Sloani_x, Sloani_y_norm)
Sloanr_yinter_Sbc = np.interp(Sbc_x, Sloanr_x, Sloanr_y_norm)
Sloanu_yinter_Sbc = np.interp(Sbc_x, Sloanu_x, Sloanu_y_norm)
Sloanz_yinter_Sbc = np.interp(Sbc_x, Sloanz_x, Sloanz_y_norm)

#Multiply the values of the SED for the galaxy by the flux values

Sbc_flux_sloani = np.multiply(Sbc_y, Sloani_yinter_Sbc)
Sbc_flux_sloang = np.multiply(Sbc_y, Sloang_yinter_Sbc)
Sbc_flux_sloanr = np.multiply(Sbc_y, Sloanr_yinter_Sbc)
Sbc_flux_sloanu = np.multiply(Sbc_y, Sloanu_yinter_Sbc)
Sbc_flux_sloanz = np.multiply(Sbc_y, Sloanz_yinter_Sbc)

#Now we are plotting the transmitted fluxes for the elliptical galaxy in each filter

pyplot.xlim(0, 20000)
pyplot.title('Sbc Galaxy flux in Sloan Filters')
pyplot.xlabel('Wavelength (Angstroms)')
pyplot.ylabel('Relative Flux (Unitless)')
pyplot.plot(Sbc_x, Sbc_flux_sloang, label = 'Sloan g')
pyplot.plot(Sbc_x, Sbc_flux_sloani, label = 'Sloan i')
pyplot.plot(Sbc_x, Sbc_flux_sloanr, label = 'Sloan r')
pyplot.plot(Sbc_x, Sbc_flux_sloanu, label = 'Sloan u')
pyplot.plot(Sbc_x, Sbc_flux_sloanz, label = 'Sloan z')
pyplot.legend()
pyplot.show();pyplot.close()

#Now that we have generated the plots, we can now find the flux in each filter, which will give us the magnitude in each filter, and from there, find the colors

Sbc_flux_g = np.trapz(Sbc_flux_sloang, Sbc_x)
Sbc_flux_i = np.trapz(Sbc_flux_sloani, Sbc_x)
Sbc_flux_r = np.trapz(Sbc_flux_sloanr, Sbc_x)
Sbc_flux_u = np.trapz(Sbc_flux_sloanu, Sbc_x)
Sbc_flux_z = np.trapz(Sbc_flux_sloanz, Sbc_x)

#Here is the magnitudes in each of the filters. 
Sbc_g = -2.5*np.log10(Sbc_flux_g)
Sbc_i = -2.5*np.log10(Sbc_flux_i)
Sbc_r = -2.5*np.log10(Sbc_flux_r)
Sbc_u = -2.5*np.log10(Sbc_flux_u)
Sbc_z = -2.5*np.log10(Sbc_flux_z)

#Compute the colors
Sbc_u_g = (Sbc_u) - (Sbc_g)
Sbc_g_r = (Sbc_g) - (Sbc_r)
Sbc_r_i = (Sbc_r) - (Sbc_i)
Sbc_i_z = (Sbc_i) - (Sbc_z)

print('Colors for Sbc Galaxy: u-g = {0:8.2f}, g-r = {1:8.2f}, r-i = {2:8.2f}, i-z = {3:8.2f}'.format(Sbc_u_g, Sbc_g_r, Sbc_r_i, Sbc_i_z))

#Finally we compute the colors for the Scd galaxy

Sloang_yinter_Scd = np.interp(Scd_x, Sloang_x, Sloang_y_norm)
Sloani_yinter_Scd = np.interp(Scd_x, Sloani_x, Sloani_y_norm)
Sloanr_yinter_Scd = np.interp(Scd_x, Sloanr_x, Sloanr_y_norm)
Sloanu_yinter_Scd = np.interp(Scd_x, Sloanu_x, Sloanu_y_norm)
Sloanz_yinter_Scd = np.interp(Scd_x, Sloanz_x, Sloanz_y_norm)

#Multiply the values of the SED for the galaxy by the flux values

Scd_flux_sloani = np.multiply(Scd_y, Sloani_yinter_Scd)
Scd_flux_sloang = np.multiply(Scd_y, Sloang_yinter_Scd)
Scd_flux_sloanr = np.multiply(Scd_y, Sloanr_yinter_Scd)
Scd_flux_sloanu = np.multiply(Scd_y, Sloanu_yinter_Scd)
Scd_flux_sloanz = np.multiply(Scd_y, Sloanz_yinter_Scd)

#Now we are plotting the transmitted fluxes for the elliptical galaxy in each filter

pyplot.xlim(0, 20000)
pyplot.title('Scd Galaxy flux in Sloan Filters')
pyplot.xlabel('Wavelength (Angstroms)')
pyplot.ylabel('Relative Flux (Unitless)')
pyplot.plot(Scd_x, Scd_flux_sloang, label = 'Sloan g')
pyplot.plot(Scd_x, Scd_flux_sloani, label = 'Sloan i')
pyplot.plot(Scd_x, Scd_flux_sloanr, label = 'Sloan r')
pyplot.plot(Scd_x, Scd_flux_sloanu, label = 'Sloan u')
pyplot.plot(Scd_x, Scd_flux_sloanz, label = 'Sloan z')
pyplot.legend()
pyplot.show();pyplot.close()

#Now that we have generated the plots, we can now find the flux in each filter, which will give us the magnitude in each filter, and from there, find the colors

Scd_flux_g = np.trapz(Scd_flux_sloang, Scd_x)
Scd_flux_i = np.trapz(Scd_flux_sloani, Scd_x)
Scd_flux_r = np.trapz(Scd_flux_sloanr, Scd_x)
Scd_flux_u = np.trapz(Scd_flux_sloanu, Scd_x)
Scd_flux_z = np.trapz(Scd_flux_sloanz, Scd_x)

#Here is the magnitudes in each of the filters. 
Scd_g = -2.5*np.log10(Scd_flux_g)
Scd_i = -2.5*np.log10(Scd_flux_i)
Scd_r = -2.5*np.log10(Scd_flux_r)
Scd_u = -2.5*np.log10(Scd_flux_u)
Scd_z = -2.5*np.log10(Scd_flux_z)

#Compute the colors
Scd_u_g = (Scd_u) - (Scd_g)
Scd_g_r = (Scd_g) - (Scd_r)
Scd_r_i = (Scd_r) - (Scd_i)
Scd_i_z = (Scd_i) - (Scd_z)

print('Colors for Scd Galaxy: u-g = {0:8.2f}, g-r = {1:8.2f}, r-i = {2:8.2f}, i-z = {3:8.2f}'.format(Scd_u_g, Scd_g_r, Scd_r_i, Scd_i_z))

#Create a table. You can print out each of them and create a table in Latex

#Which filter best distingshes the galaxy types at redshift zero?
#Look at the wavelengths of the spectral peaks?

#Color in terms of spectral populations. 
#Different stars give out different light colors, if, for instance red colors indicate that there are older stars and less star forming regions. Similarly,more blue colors are from younger more massive stars, but be weary, as these bluer, more massive stars are more luminous than their old, red counterparts, so you might not able to claim that younger stars hold more of the population. 

###############################################################################

#We want to find the colors at different redshifts and write it all to an ASCII file. Create the ASCII file first and then create the loops and at the end, instead of 'print', write it to the file. 

#variable = open ('filename', 'w+ = write with override"
outputascii = open('Homework1ascii' , 'w+')

#varialbe.write('writes string use \n to end line')
outputascii.write('Galaxy type | Redshift | u-g | g-r | r-i | i-z \n')

#Computing the colors of the Elliptical galaxy
#assign variable to redshift value

redshift = [ 0. , 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45, 1.50]

#for each redshift value in the array do ________

for x in redshift:

#Use the code you used in part II to find the colors, except this time, note the first line, where we change the values of the wavelengths for each iteration in the loop. 

	elliptical_x_red = np.array(elliptical_x, dtype = float)*(1 + x)
	Sloang_yinter_el_red = np.interp(elliptical_x_red, Sloang_x, Sloang_y_norm)
	Sloani_yinter_el_red = np.interp(elliptical_x_red, Sloani_x, Sloani_y_norm)
	Sloanr_yinter_el_red = np.interp(elliptical_x_red, Sloanr_x, Sloanr_y_norm)
	Sloanu_yinter_el_red = np.interp(elliptical_x_red, Sloanu_x, Sloanu_y_norm)
	Sloanz_yinter_el_red = np.interp(elliptical_x_red, Sloanz_x, Sloanz_y_norm)

	elliptical_flux_sloani_red = np.multiply(elliptical_y, Sloani_yinter_el_red)
	elliptical_flux_sloang_red = np.multiply(elliptical_y, Sloang_yinter_el_red)
	elliptical_flux_sloanr_red = np.multiply(elliptical_y, Sloanr_yinter_el_red)
	elliptical_flux_sloanu_red = np.multiply(elliptical_y, Sloanu_yinter_el_red)
	elliptical_flux_sloanz_red = np.multiply(elliptical_y, Sloanz_yinter_el_red)

	elliptical_flux_g_red = np.trapz(elliptical_flux_sloang_red, elliptical_x_red)
	elliptical_flux_i_red = np.trapz(elliptical_flux_sloani_red, elliptical_x_red)
	elliptical_flux_r_red = np.trapz(elliptical_flux_sloanr_red, elliptical_x_red)
	elliptical_flux_u_red = np.trapz(elliptical_flux_sloanu_red, elliptical_x_red)
	elliptical_flux_z_red = np.trapz(elliptical_flux_sloanz_red, elliptical_x_red)

	El_g_red = -2.5*np.log10(elliptical_flux_g_red)
	El_i_red = -2.5*np.log10(elliptical_flux_i_red)
	El_r_red = -2.5*np.log10(elliptical_flux_r_red)
	El_u_red = -2.5*np.log10(elliptical_flux_u_red)
	El_z_red = -2.5*np.log10(elliptical_flux_z_red)

	elliptical_u_g_red = (El_u_red) - (El_g_red)
	elliptical_g_r_red = (El_g_red) - (El_r_red)
	elliptical_r_i_red = (El_r_red) - (El_i_red)
	elliptical_i_z_red = (El_i_red) - (El_z_red)

	outputascii.write('Elliptical, {0:8.3}, {1:8.2f}, {2:8.2f}, {3:8.2f}, {4:8.2f} \n'.format(x, elliptical_u_g_red, elliptical_g_r_red, elliptical_r_i_red, elliptical_i_z_red)) 

#Now, find the colors for the Irregular galaxy

for x in redshift:
	irregular_x_red = np.array(irregular_x, dtype = float)*(1 + x)
	Sloang_yinter_ir_red = np.interp(irregular_x_red, Sloang_x, Sloang_y_norm)
	Sloani_yinter_ir_red = np.interp(irregular_x_red, Sloani_x, Sloani_y_norm)
	Sloanr_yinter_ir_red = np.interp(irregular_x_red, Sloanr_x, Sloanr_y_norm)
	Sloanu_yinter_ir_red = np.interp(irregular_x_red, Sloanu_x, Sloanu_y_norm)
	Sloanz_yinter_ir_red = np.interp(irregular_x_red, Sloanz_x, Sloanz_y_norm)

	irregular_flux_sloani_red = np.multiply(irregular_y, Sloani_yinter_ir_red)
	irregular_flux_sloang_red = np.multiply(irregular_y, Sloang_yinter_ir_red)
	irregular_flux_sloanr_red = np.multiply(irregular_y, Sloanr_yinter_ir_red)
	irregular_flux_sloanu_red = np.multiply(irregular_y, Sloanu_yinter_ir_red)
	irregular_flux_sloanz_red = np.multiply(irregular_y, Sloanz_yinter_ir_red)

	irregular_flux_g_red = np.trapz(irregular_flux_sloang_red, irregular_x_red)
	irregular_flux_i_red = np.trapz(irregular_flux_sloani_red, irregular_x_red)
	irregular_flux_r_red = np.trapz(irregular_flux_sloanr_red, irregular_x_red)
	irregular_flux_u_red = np.trapz(irregular_flux_sloanu_red, irregular_x_red)
	irregular_flux_z_red = np.trapz(irregular_flux_sloanz_red, irregular_x_red) 

	ir_g_red = -2.5*np.log10(irregular_flux_g_red)
	ir_i_red = -2.5*np.log10(irregular_flux_i_red)
	ir_r_red = -2.5*np.log10(irregular_flux_r_red)
	ir_u_red = -2.5*np.log10(irregular_flux_u_red)
	ir_z_red = -2.5*np.log10(irregular_flux_z_red)

	irregular_u_g_red = (ir_u_red) - (ir_g_red)
	irregular_g_r_red = (ir_g_red) - (ir_r_red)
	irregular_r_i_red = (ir_r_red) - (ir_i_red)
	irregular_i_z_red = (ir_i_red) - (ir_z_red)

	outputascii.write('Irregular, {0:8.3}, {1:8.2f}, {2:8.2f}, {3:8.2f}, {4:8.2f} \n'.format(x, irregular_u_g_red, irregular_g_r_red, irregular_r_i_red, irregular_i_z_red)) 

#Now, find the colors for the Sbc galaxy

for x in redshift:
	Sbc_x_red = np.array(Sbc_x, dtype = float)*(1 + x)
	Sloang_yinter_sbc_red = np.interp(Sbc_x_red, Sloang_x, Sloang_y_norm)
	Sloani_yinter_sbc_red = np.interp(Sbc_x_red, Sloani_x, Sloani_y_norm)
	Sloanr_yinter_sbc_red = np.interp(Sbc_x_red, Sloanr_x, Sloanr_y_norm)
	Sloanu_yinter_sbc_red = np.interp(Sbc_x_red, Sloanu_x, Sloanu_y_norm)
	Sloanz_yinter_sbc_red = np.interp(Sbc_x_red, Sloanz_x, Sloanz_y_norm)

	sbc_flux_sloani_red = np.multiply(Sbc_y, Sloani_yinter_sbc_red)
	sbc_flux_sloang_red = np.multiply(Sbc_y, Sloang_yinter_sbc_red)
	sbc_flux_sloanr_red = np.multiply(Sbc_y, Sloanr_yinter_sbc_red)
	sbc_flux_sloanu_red = np.multiply(Sbc_y, Sloanu_yinter_sbc_red)
	sbc_flux_sloanz_red = np.multiply(Sbc_y, Sloanz_yinter_sbc_red)

	sbc_flux_g_red = np.trapz(sbc_flux_sloang_red, Sbc_x_red)
	sbc_flux_i_red = np.trapz(sbc_flux_sloani_red, Sbc_x_red)
	sbc_flux_r_red = np.trapz(sbc_flux_sloanr_red, Sbc_x_red)
	sbc_flux_u_red = np.trapz(sbc_flux_sloanu_red, Sbc_x_red)
	sbc_flux_z_red = np.trapz(sbc_flux_sloanz_red, Sbc_x_red) 

	sbc_g_red = -2.5*np.log10(sbc_flux_g_red)
	sbc_i_red = -2.5*np.log10(sbc_flux_i_red)
	sbc_r_red = -2.5*np.log10(sbc_flux_r_red)
	sbc_u_red = -2.5*np.log10(sbc_flux_u_red)
	sbc_z_red = -2.5*np.log10(sbc_flux_z_red)

	sbc_u_g_red = (sbc_u_red) - (sbc_g_red)
	sbc_g_r_red = (sbc_g_red) - (sbc_r_red)
	sbc_r_i_red = (sbc_r_red) - (sbc_i_red)
	sbc_i_z_red = (sbc_i_red) - (sbc_z_red)

	outputascii.write('Sbc, {0:8.3}, {1:8.2f}, {2:8.2f}, {3:8.2f}, {4:8.2f} \n'.format(x, sbc_u_g_red, sbc_g_r_red, sbc_r_i_red, sbc_i_z_red)) 

#Now, find the colors for the Sbc galaxy

for x in redshift:
	Scd_x_red = np.array(Scd_x, dtype = float)*(1 + x)
	Sloang_yinter_scd_red = np.interp(Scd_x_red, Sloang_x, Sloang_y_norm)
	Sloani_yinter_scd_red = np.interp(Scd_x_red, Sloani_x, Sloani_y_norm)
	Sloanr_yinter_scd_red = np.interp(Scd_x_red, Sloanr_x, Sloanr_y_norm)
	Sloanu_yinter_scd_red = np.interp(Scd_x_red, Sloanu_x, Sloanu_y_norm)
	Sloanz_yinter_scd_red = np.interp(Scd_x_red, Sloanz_x, Sloanz_y_norm)

	scd_flux_sloani_red = np.multiply(Scd_y, Sloani_yinter_scd_red)
	scd_flux_sloang_red = np.multiply(Scd_y, Sloang_yinter_scd_red)
	scd_flux_sloanr_red = np.multiply(Scd_y, Sloanr_yinter_scd_red)
	scd_flux_sloanu_red = np.multiply(Scd_y, Sloanu_yinter_scd_red)
	scd_flux_sloanz_red = np.multiply(Scd_y, Sloanz_yinter_scd_red)

	scd_flux_g_red = np.trapz(scd_flux_sloang_red, Scd_x_red)
	scd_flux_i_red = np.trapz(scd_flux_sloani_red, Scd_x_red)
	scd_flux_r_red = np.trapz(scd_flux_sloanr_red, Scd_x_red)
	scd_flux_u_red = np.trapz(scd_flux_sloanu_red, Scd_x_red)
	scd_flux_z_red = np.trapz(scd_flux_sloanz_red, Scd_x_red) 

	scd_g_red = -2.5*np.log10(scd_flux_g_red)
	scd_i_red = -2.5*np.log10(scd_flux_i_red)
	scd_r_red = -2.5*np.log10(scd_flux_r_red)
	scd_u_red = -2.5*np.log10(scd_flux_u_red)
	scd_z_red = -2.5*np.log10(scd_flux_z_red)

	scd_u_g_red = (scd_u_red) - (scd_g_red)
	scd_g_r_red = (scd_g_red) - (scd_r_red)
	scd_r_i_red = (scd_r_red) - (scd_i_red)
	scd_i_z_red = (scd_i_red) - (scd_z_red)

	outputascii.write('Scd, {0:8.3}, {1:8.2f}, {2:8.2f}, {3:8.2f}, {4:8.2f} \n'.format(x, scd_u_g_red, scd_g_r_red, scd_r_i_red, scd_i_z_red))

outputascii.close() 

#Create a color - color diagram

#opened the ascii in a spreadsheet program and then changed it to a csv file 

#First read in the data from the csv file
data = np.genfromtxt('Photometric_redshifts.csv', dtype = float, delimiter = ',')

#Now, create arrays to use for the data
redshifts = data[1:31, 1]
elliptical_u_g_mag = data[1:31, 2]
elliptical_g_r_mag = data[1:31, 3]
elliptical_r_i_mag = data[1:31, 4]
elliptical_i_z_mag = data[1:31, 5]

irregular_u_g_mag = data[32:62, 2]
irregular_g_r_mag = data[32:62, 3]
irregular_r_i_mag = data[32:62, 4]
irregular_i_z_mag = data[32:62, 5]

Sbc_u_g_mag = data[63:93, 2]
Sbc_g_r_mag = data[63:93, 3]
Sbc_r_i_mag = data[63:93, 4]
Sbc_i_z_mag = data[63:93, 5]

Scd_u_g_mag = data[94:124, 2]
Scd_g_r_mag = data[94:124, 3]
Scd_r_i_mag = data[94:124, 4]
Scd_i_z_mag = data[94:124, 5]

#Start to plot

#Plot for u-g vs g-r

pyplot.title('Color vs Color Diagram for u-g vs g-r')
pyplot.scatter(elliptical_u_g_mag, elliptical_g_r_mag, label = 'Elliptical Galaxy')
pyplot.scatter(irregular_u_g_mag, irregular_g_r_mag, label = 'Irregular Galaxy')
pyplot.scatter(Sbc_u_g_mag, Sbc_g_r_mag, label = 'Sbc Galaxy')
pyplot.scatter(Scd_u_g_mag, Scd_g_r_mag, label = 'Scd Galaxy')
pyplot.xlabel('Magnitude (u-g)') 
pyplot.ylabel('Magnitude (g-r)')
pyplot.legend()
pyplot.grid()
pyplot.show();pyplot.close()

#Plot for g-r vs r-i

pyplot.title('Color vs Color Diagram for g-r vs r-i')
pyplot.scatter(elliptical_g_r_mag, elliptical_r_i_mag, label = 'Elliptical Galaxy')
pyplot.scatter(irregular_g_r_mag, irregular_r_i_mag, label = 'Irregular Galaxy')
pyplot.scatter(Sbc_g_r_mag, Sbc_r_i_mag, label = 'Sbc Galaxy')
pyplot.scatter(Scd_g_r_mag, Scd_r_i_mag, label = 'Scd Galaxy')
pyplot.xlabel('Magnitude (g-r)') 
pyplot.ylabel('Magnitude (r-i)')
pyplot.legend()
pyplot.grid()
pyplot.show();pyplot.close()

#Plot for r-i vs i-z
pyplot.scatter(elliptical_r_i_mag, elliptical_i_z_mag, label = 'Elliptical Galaxy')
pyplot.scatter(irregular_r_i_mag, irregular_i_z_mag, label = 'Irregular Galaxy')
pyplot.scatter(Sbc_r_i_mag, Sbc_i_z_mag, label = 'Sbc Galaxy')
pyplot.scatter(Scd_r_i_mag, Scd_i_z_mag, label = 'Scd Galaxy')
pyplot.xlabel('Magnitude (r-i)') 
pyplot.ylabel('Magnitude (i-z)')
pyplot.title('Colo r vs Color Diagram for r-i vs i-z')
pyplot.legend()
pyplot.grid()
pyplot.show();pyplot.close()

#Try creatiing a 3-D plot
fig = pyplot.figure()
ax = fig.add_subplot(111, projection = '3d')

#Create the u-g vs g-r vs redshift
ax.scatter(elliptical_u_g_mag, redshifts, elliptical_g_r_mag, label = 'Elliptical Galaxy')
ax.scatter(irregular_u_g_mag, redshifts, irregular_g_r_mag, label = 'Irregular Galaxy')
ax.scatter(Sbc_u_g_mag, redshifts, Sbc_g_r_mag, label = 'Sbc Galaxy')
ax.scatter(Scd_u_g_mag, redshifts, Scd_g_r_mag, label = 'Scd Galaxy')
#pyplot.legend()
#pyplot.xlabel('u-g (Magnitude)')
#pyplot.ylabel('g-r (Magnitude)')
pyplot.show();pyplot.close()


##############################################################################

###############################Scratch work###################################


#Starting off with the wavelengths for Sloan g'

#Sloan_g_cut = ([Sloan_x[1001:1502] [Sloan_y[1001:1502]])

#First scale the wavelength values of the plots by 0.1 (conversion from A to nm)

#pyplot.plot(elliptical_x_nm, elliptical_y)
#pyplot.show();pyplot.close()

#Since the we want to use the same wavelengths from the sloan filters, Since the units are different, we change the elliptical_x from Angstroms to n
#Starting with the elliptical galaxies, we want to interpolate the wavelengths so that they match with the Sloan ones. This creates an 1D array which we can multiply the transmission with the flux from the wavelengths. 

#interpolated_elliptical_fluxes= np.interp(Sloang_x, elliptical_x, elliptical_y)
#interpolated_irregular = np.interp(interpolatedwl, irregular_x, irregular_y)
#interpolated_Sbc = np.interp(interpolatedwl, Sbc_x, Sbc_y)
#interpolated_Scd = np.interp(interpolatedwl, Scd_x, Scd_y)

#pyplot.plot(Sloang_x, interpolated_elliptical_fluxes)
#pyplot.show();pyplot.close()

#Now that we have the interpolated values, we want to multiply the transmission functions from the Sloan filters by the new fluxes
#newellipticalflux_g = np.multiply(interpolated_elliptical[:,1], Sloang_y)
#pyplot.plot(interpolated_elliptical[:, 0], newellipticalflux_g[:, 1])
#pyplot.show();pyplot.close()

#For the interpolated flux in the galaxies

#flux_elliptical = np.trapz(elliptical_y, elliptical_x)
#flux_elliptical

#Export data into ascii file
#fp = open('Homework1_ascii', 'w')
#fp.write('Colors for Elliptical Galaxy for redshift {0:8.3}: u-g = {1:8.2f}, g-r = {2:8.2f}, r-i = {3:8.2f}, i-z = {4:8.2f}'.format(x, elliptical_u_g_red, elliptical_g_r_red, elliptical_r_i_red, elliptical_i_z_red)) \n
#fp.write('Colors for Irregular Galaxy for redshift {0:8.3}: u-g = {1:8.2f}, g-r = {2:8.2f}, r-i = {3:8.2f}, i-z = {4:8.2f}'.format(x, irregular_u_g_red, irregular_g_r_red, irregular_r_i_red, irregular_i_z_red)) \n
#fp.write('Colors for Sbc Galaxy for redshift {0:8.3}: u-g = {1:8.2f}, g-r = {2:8.2f}, r-i = {3:8.2f}, i-z = {4:8.2f}'.format(x, sbc_u_g_red, sbc_g_r_red, sbc_r_i_red, sbc_i_z_red)) \n
#fp.write('Colors for Scd Galaxy for redshift {0:8.3}: u-g = {1:8.2f}, g-r = {2:8.2f}, r-i = {3:8.2f}, i-z = {4:8.2f}'.format(x, scd_u_g_red, scd_g_r_red, scd_r_i_red, scd_i_z_red)) \n
#fp.close() 



