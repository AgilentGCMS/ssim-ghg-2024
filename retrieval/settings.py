import numpy as np
import os

# Default settings

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
d_aerosol = 3.0e-6 #Dust-like particle diameter [m]
n_aerosol = 1.4 + 0.0003j #Dust-like refractive index
tau_aerosol_true = 0.1 #AOD in the first band (O2 A-band at ~0.76 um)


############################################
#True settings

co2_true = np.array([0.00040774, 0.00041397, 0.00041681, 0.00041669, 0.00041656, 0.00041641, 0.00041623, 0.00041603, 0.00041578, 0.00041557]) #mol/mol
ch4_true = np.full((len(co2_true)),1900e-9) #mol/mol
T_true = np.array([231.00832168, 221.11536418, 222.6548729 , 232.87167357, 243.33909691, 250.17586371, 256.20882439, 260.21844782, 265.51743697, 270.87057495]) #K
p_true = np.array([ 10000.,  20000.,  30000.,  40000.,  50000.,  60000.,  70000., 80000.,  90000., 100000.])
q_true = np.array([2.97593695e-06, 9.51235205e-06, 3.97622333e-05, 1.46971198e-04, 3.55363471e-04, 5.34367416e-04, 8.39863516e-04, 1.50594464e-03, 2.29157672e-03, 2.66968086e-03])


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

#Create a prior covariance matrix
S_prior = np.zeros((9,9))
np.fill_diagonal(S_prior,[co2_prior_uncert**2,ch4_prior_uncert**2,T_prior_uncert**2,p_prior_uncert**2,q_prior_uncert**2,albedo_uncert[0]**2,albedo_uncert[1]**2,albedo_uncert[2]**2,tau_aerosol_prior_uncert**2])

#How good should our modeled fit of the radiances be before we stop the retrieval?
chisq_threshold = 1.01

#Used in the finite differencing calculation for T and p
perturbation = 0.001


