import numpy as np
from find_nearest import find_nearest
import os

#Default settings

# Where are the ABSCO tables kept?
ABSCO_TABLE_FOLDER = os.path.join(os.environ['HOME'], 'shared/ssim-ghg-data/retrieval_example')

#0.76 um O2 A-band, 1.60 um weak CO2 band, 1.65 um CH4 band
band_max_wn = np.array([13190.0,6285.0,6150.0]) #cm^-1
band_min_wn = np.array([12950.0,6166.0,5995.0]) #cm^-1

#Set the spectral resolution (FWHM) to be OCO-2-like
band_spectral_resolutions = np.array([0.042, 0.076, 0.076]) #nm

#What molecules do we care about for each band?
band_molecules = []
band_molecules.append(["o2","h2o"])
band_molecules.append(["co2","h2o"])
band_molecules.append(["ch4","h2o"])

#Create a "true" scene
sza_0 = 30.0 #solar zenith angle, degrees
sza   = 0.0 #sensor zenith angle, degrees
albedo_true = np.array([0.20, 0.25, 0.26]) #per-band albedo

#Signal-to-noise ratio
SNR = 400.0

#Aerosol properties
height_aerosol = 80000 #Assume the aerosol layer goes from the surface to this pressure in Pa
d_aerosol = 3.0e-6 #Particle diameter [m]
n_aerosol = 1.4 + 0.0003j #Refractive index
tau_aerosol_true = 0.0 #AOD in the first band (O2 A-band at ~0.76 um)


############################################
#True settings

p_true = np.array([ 0.01, 10000.,  20000.,  30000.,  40000.,  50000.,  60000.,  70000., 80000.,  90000., 100000.])
co2_true = np.array([0.0004032 , 0.00041564, 0.00041696, 0.00041686, 0.00041675, 0.00041665, 0.00041653, 0.00041641, 0.00041627, 0.00041612, 0.00041595]) #mol/mol
ch4_true = np.array([500e-9,1500e-9,1900e-9,1900e-9,1900e-9,1900e-9,1900e-9,1900e-9,1900e-9,1900e-9,1900e-9]) #mol/mol
T_true = np.array([205.59104919, 203.52311895, 220.08876161, 242.37260452, 255.89316285, 266.14411523, 271.82797304, 277.33967816, 279.68627102, 283.79536257, 291.02379557]) #K
q_true = np.array([3.17349759e-06, 3.60197244e-06, 5.41395304e-05, 6.18924547e-04, 1.51411285e-03, 8.10448167e-04, 4.89113692e-04, 1.82153489e-03, 7.49874173e-03, 8.96673750e-03, 1.03688360e-02]) #kg/kg


############################################
#A priori settings

#Let's say we have an a priori guess that's somewhat close to the truth, and that we can get back to the truth by adding/subtracting or multiplying/dividing the state vector
co2_prior = co2_true * 0.99 #Decrease the prior CO2 profile by 1%
ch4_prior = ch4_true * 1.05 #Increase the prior CH4 profile by 5%
T_prior = T_true + 2.0 #Add 2 K to each level of the true temperature profile
p_prior = p_true * 1.02 #Increase the true pressure profile by 2%
q_prior = q_true * 1.05 #Increase the true specific humidity profile by 5%
albedo_prior = [albedo_true[0] + 0.02,albedo_true[1] + 0.02,albedo_true[2] + 0.02] #Add 0.02 to the true albedo
tau_aerosol_prior = 0.05 #Start with a near-zero AOD. 


##############################################
#Retrieval settings
#Need to assign some 1-sigma a priori uncertainty to our state vector
co2_prior_uncert = 0.1 #10%
ch4_prior_uncert = 0.2 #20%
T_prior_uncert = 5 #K
p_prior_uncert = 0.1 #10%
q_prior_uncert = 0.1 #10%
albedo_uncert = [0.2,0.2,0.2]
tau_aerosol_prior_uncert = 0.5

#Create a prior covariance matrix (no aerosols by default)
S_prior = np.zeros((8,8))
np.fill_diagonal(S_prior,[co2_prior_uncert**2,ch4_prior_uncert**2,T_prior_uncert**2,p_prior_uncert**2,q_prior_uncert**2,albedo_uncert[0]**2,albedo_uncert[1]**2,albedo_uncert[2]**2])

#How good should our modeled fit of the radiances be before we stop the retrieval?
chisq_threshold = 1.01

#Used in the finite differencing calculation for T and p
perturbation = 0.001

##############################################
#Other spectral stuff that we only want to do once to save time

#Band limits but in micrometers. Min/max flipped because we're going from wavenumber to wavelength
band_min_um = 1e4/band_max_wn
band_max_um = 1e4/band_min_wn

#Using the spectral resolution in nm, convert to um and assume 3 channels per FWHM to get the number of spectral points in each band
band_spectral_points = np.empty((len(band_max_wn)),dtype=int)
for i in range(len(band_max_wn)):
    band_spectral_points[i] = int(1000.*3./band_spectral_resolutions[i] * (band_max_um[i] - band_min_um[i]))

#Create evenly spaced channels for each band for use in the spectral response function
band_wn = []
for i in range(len(band_max_wn)):
    band_wn.append(np.linspace(band_min_wn[i],band_max_wn[i],band_spectral_points[i]))

#Convert bands from wavenumber to wavelength in um
band_wl = []
for i in range(len(band_max_wn)):
    band_wl.append((1.e4/band_wn[i])[::-1])

#Channels in wavenumber but at ABSCO resolution (0.01 cm^-1)
band_absco_res_wn = []
for i in range(len(band_min_wn)):
    band_absco_res_wn.append(np.linspace(band_min_wn[i],band_max_wn[i],int((band_max_wn[i]-band_min_wn[i])/0.01) + 1))

#Calculate the resolving power of each band
resolving_power_band = []
for i in range(len(band_min_wn)):
    resolving_power_band.append(np.mean(band_wl[i])/band_spectral_resolutions[i]*1e3)

#Calculate standard deviations via FWHM of each band
sigma_band = []
for i in range(len(band_min_wn)):
    sigma_band.append(np.median(band_absco_res_wn[i])/resolving_power_band[i]/(2.0*((2.0*np.log(2.0))**(0.5))))

#Find the indicies of each instrument-res wavenumber in the high-res wavenumber array
band_wn_index = []
for i in range(len(band_min_wn)):
    band_wn_index.append(find_nearest(band_absco_res_wn[i],band_wn[i]))

#Calculate a big term used in the ILS function to save time
ILS_Gaussian_term = []
ILS_Gaussian_term_sum = []
for i in range(len(band_min_wn)):
  ILS_Gaussian_term.append(1.0/sigma_band[i]/((2.0*np.pi)**0.5) * np.exp(-(band_absco_res_wn[i][None,:]-band_absco_res_wn[i][band_wn_index[i][:,None]])**2.0 / 2.0 / sigma_band[i]**2.0))
  ILS_Gaussian_term_sum.append(np.sum(ILS_Gaussian_term[i],axis=1))

#Some default aerosol stuff
qext_aerosol = []
ssa_aerosol = []
P_aerosol = []
for i in range(len(band_min_wn)):
  qext_aerosol.append(np.zeros((len(band_absco_res_wn[i]))))
  ssa_aerosol.append(np.zeros((len(band_absco_res_wn[i]))))
  P_aerosol.append(np.zeros((len(band_absco_res_wn[i]))))

