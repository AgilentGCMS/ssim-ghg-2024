import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
import sys
import os
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
import settings as s
from absco_lookup import sigma_lookup
from copy import deepcopy
import time

#Constants
g=9.81 #m/s^2
M=.0289644 #kg/mol
Na=6.022e23 #molecule/mol^-1
sigma_sun = 6.794e-5 #sr
h = 6.62607015e-34 #m^2*kg/s
c = 2.99792458e8 #m/s
k = 1.381e-23 #JK^-1molecule^-1

#Constants for Rayleigh scattering calculations
a_rayleigh = 2.871e-4
b_rayleigh = 5.67e-3
rho_rayleigh = 0.0279 #depolarization factor
epsilon = 0.622

#Other variables that we're hard-coding.
ILS_width = 5.0


class ForwardFunction:
    def __init__(self,SNR=s.SNR,sza_0=s.sza_0,sza=s.sza,co2=s.co2_true,ch4=s.ch4_true,T=s.T_true,p=s.p_true,q=s.q_true,albedo=s.albedo_true,band_min_wn=s.band_min_wn,band_max_wn=s.band_max_wn,band_spectral_resolutions=s.band_spectral_resolutions,band_min_um=s.band_min_um,band_max_um=s.band_max_um,band_spectral_points=s.band_spectral_points,band_wn=s.band_wn,band_wl=s.band_wl,band_absco_res_wn=s.band_absco_res_wn,resolving_power_band=s.resolving_power_band,sigma_band=s.sigma_band,band_wn_index=s.band_wn_index,ILS_Gaussian_term=s.ILS_Gaussian_term,ILS_Gaussian_term_sum=s.ILS_Gaussian_term_sum,absco_data=None,band_molecules=s.band_molecules,P_aerosol=s.P_aerosol,ssa_aerosol=s.ssa_aerosol,qext_aerosol=s.qext_aerosol,height_aerosol=s.height_aerosol,tau_aerosol=None,measurement_error=False,jacobians=False):


        self.SNR = SNR
        self.sza_0 = sza_0
        self.sza = sza
        self.co2 = co2
        self.ch4 = ch4
        self.T = T
        self.p = p
        self.q = q
        self.albedo = albedo
        self.band_min_wn = band_min_wn
        self.band_max_wn = band_max_wn
        self.band_spectral_resolutions = band_spectral_resolutions
        self.band_molecules = band_molecules
        self.P_aerosol = P_aerosol
        self.ssa_aerosol = ssa_aerosol
        self.qext_aerosol = qext_aerosol
        self.height_aerosol = height_aerosol
        self.tau_aerosol = tau_aerosol
        self.measurement_error = measurement_error
        self.jacobians = jacobians
        self.band_min_um = band_min_um
        self.band_max_um = band_max_um
        self.band_spectral_points = band_spectral_points
        self.band_wn = band_wn
        self.band_wl = band_wl
        self.band_absco_res_wn = band_absco_res_wn
        self.resolving_power_band = resolving_power_band
        self.sigma_band = sigma_band
        self.band_wn_index = band_wn_index
        self.ILS_Gaussian_term = ILS_Gaussian_term
        self.ILS_Gaussian_term_sum = ILS_Gaussian_term_sum

        #Approximately calculate the solar irradiance at our bands using Planck's law (assuming the Sun is 5800 K), then account for the solid angle of the Sun and convert into per um instead of per m.
        self.band_solar_irradiances = np.empty((len(self.band_min_wn)))
        for i in range(len(self.band_max_wn)):
            self.band_solar_irradiances[i] = planck(5800.,np.mean(self.band_wl[i])*1e-6)*sigma_sun/1.e6 #W/m^2/um

        #Geometry
        self.mu_0=np.cos(np.deg2rad(self.sza_0)) #cosine of the solar zenith angle [deg]
        self.mu=np.cos(np.deg2rad(self.sza)) #cosine of the sensor zenith angle [deg]
        self.m = 1/self.mu_0+1/self.mu #airmass

        #Calculate the d(pressure) of a given layer
        self.p_diff=np.empty((len(self.p)-1))
        for i in range(len(self.p)-1):
            self.p_diff[i] = self.p[i+1]-self.p[i]

        #Layer calculations
        self.co2_layer = layer_average(self.co2)
        self.ch4_layer = layer_average(self.ch4)
        self.T_layer = layer_average(self.T)
        self.p_layer = layer_average(self.p)
        self.q_layer = layer_average(self.q)

        #Calculate the true Xgases from the gas profiles
        self.xco2, self.h = calculate_Xgas(self.co2, self.p, self.q)
        self.xch4, self.h = calculate_Xgas(self.ch4, self.p, self.q)

        self.xco2 *= 1e6 #Convert from mol/mol to ppm
        self.xch4 *= 1e9 #Convert from mol/mol to ppb

        ##########################################################
        #Calculate absorption cross sections for H2O, O2, and CO2 in m^2/molecule, as needed (i.e., there's no CO2 absorption in the OCO-2 O2 A-band)
        self.tau_star_band = [] #band_n total optical depths
        self.tau_above_aerosol_star_band = [] #band_n total optical depths above the aerosol layer

        #For analytical Jacobians. Calculate on each layer.
        self.tau_star_band_q = []
        self.tau_above_aerosol_star_band_q = []
        self.tau_star_band_co2 = []
        self.tau_above_aerosol_star_band_co2 = []
        self.tau_star_band_ch4 = []
        self.tau_above_aerosol_star_band_ch4 = []

        #Loop through the bands
        for i in range(len(self.band_min_um)):

            #Loop through the desired molecules for this band
            tau_star_temp = np.zeros((len(self.band_absco_res_wn[i])))
            tau_above_aerosol_star_temp = np.zeros((len(self.band_absco_res_wn[i])))
     
            #For analytical Jacobians (need q, co2, ch4)
            tau_star_temp_q = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_layer)))
            tau_above_aerosol_star_temp_q = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_layer)))
            tau_star_temp_co2 = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_layer)))
            tau_above_aerosol_star_temp_co2 = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_layer)))
            tau_star_temp_ch4 = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_layer)))
            tau_above_aerosol_star_temp_ch4 = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_layer)))

            for j in range(len(self.band_molecules[i])):
                molecules_sigma_temp = sigma_lookup(self.band_molecules[i][j],self.band_absco_res_wn[i],self.p_layer,self.T_layer,absco_data)

                #Tau is the cross section times number density times dz. We can assume hydrostatic equilibrium and the ideal gas law to solve it using this equation instead:
                if self.band_molecules[i][j] == 'o2':
                    tau_temp = np.tile(self.p_diff,(len(self.band_absco_res_wn[i]),1)) * 0.20935 * molecules_sigma_temp / M / g
                    if self.jacobians:
                        tau_temp_q = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_diff))) 
                        tau_temp_co2 = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_diff)))
                        tau_temp_ch4 = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_diff))) 

                elif self.band_molecules[i][j] == 'h2o': 
                    tau_temp = np.tile(self.p_diff,(len(self.band_absco_res_wn[i]),1)) * np.tile(self.q_layer,(len(molecules_sigma_temp),1)) * molecules_sigma_temp / M / g
                    if self.jacobians:
                        tau_temp_q = np.tile(self.p_diff,(len(self.band_absco_res_wn[i]),1)) * molecules_sigma_temp / M / g 
                        tau_temp_co2 = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_diff)))
                        tau_temp_ch4 = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_diff)))

                elif self.band_molecules[i][j] == 'co2': 
                    tau_temp = np.tile(self.p_diff,(len(self.band_absco_res_wn[i]),1)) * np.tile(self.co2_layer,(len(molecules_sigma_temp),1)) * molecules_sigma_temp / M / g
                    if self.jacobians:
                        tau_temp_q = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_diff)))
                        tau_temp_co2 = np.tile(self.p_diff,(len(self.band_absco_res_wn[i]),1)) * molecules_sigma_temp / M / g 
                        tau_temp_ch4 = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_diff)))

                elif self.band_molecules[i][j] == 'ch4': 
                    tau_temp = np.tile(self.p_diff,(len(self.band_absco_res_wn[i]),1)) * np.tile(self.ch4_layer,(len(molecules_sigma_temp),1)) * molecules_sigma_temp / M / g
                    if self.jacobians:
                        tau_temp_q = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_diff)))
                        tau_temp_co2 = np.zeros((len(self.band_absco_res_wn[i]),len(self.p_diff)))
                        tau_temp_ch4 = np.tile(self.p_diff,(len(self.band_absco_res_wn[i]),1)) * molecules_sigma_temp / M / g 

                else:
                    print("Choose a valid molecule!")
                    return

                #Sum the vertical profile dimension
                tau_star_temp += np.sum(tau_temp,axis=1)

                #Sum the vertical profile dimension but only above the aerosol layer
                tau_above_aerosol_star_temp += np.sum(tau_temp[:,self.p_layer < self.height_aerosol],axis=1)

                if self.jacobians:
                    #For analytic Jacobians, save each layer
                    tau_star_temp_q += tau_temp_q
                    tau_star_temp_co2 += tau_temp_co2
                    tau_star_temp_ch4 += tau_temp_ch4

                    #For analytic Jacobians, only add the layer aod above height_aerosol
                    tau_above_aerosol_star_temp_q[:,self.p_layer < self.height_aerosol] += tau_temp_q[:,self.p_layer < self.height_aerosol]
                    tau_above_aerosol_star_temp_co2[:,self.p_layer < self.height_aerosol] += tau_temp_co2[:,self.p_layer < self.height_aerosol]
                    tau_above_aerosol_star_temp_ch4[:,self.p_layer < self.height_aerosol] += tau_temp_ch4[:,self.p_layer < self.height_aerosol]


            #Also add the Rayleigh scattering optical depth  
            tau_rayleigh_band = np.empty((len(self.band_absco_res_wn[i]),len(self.p_layer)))

            #See Section 3.2.1.5 in the OCO-2 L2 ATBD
            n_s = 1 + a_rayleigh*(1. + b_rayleigh*np.mean(self.band_wl[i])**2.) #Just use Rayleigh in the middle of the O2 A-band (in um here)

            #Typo in the ATBD. Ns should be e-24 instead of e-20.
            rayleigh_sigma_band = (1.031e-24 * (n_s**2. - 1.)**2.)/(np.mean(self.band_wl[i])**4. * (n_s**2. + 2.)**2)*(6.+3.*rho_rayleigh)/(6.-7.*rho_rayleigh)
            #Simplify
            tau_rayleigh_band = self.p_diff * Na * rayleigh_sigma_band / M / g

            tau_star_temp += np.sum(tau_rayleigh_band)
            tau_above_aerosol_star_temp += np.sum(tau_rayleigh_band[self.p_layer < self.height_aerosol])

            #Append for the band we're on
            self.tau_star_band.append(tau_star_temp)
            self.tau_above_aerosol_star_band.append(tau_above_aerosol_star_temp)

            #For analytic Jacobians
            self.tau_star_band_q.append(tau_star_temp_q)
            self.tau_above_aerosol_star_band_q.append(tau_above_aerosol_star_temp_q)
            self.tau_star_band_co2.append(tau_star_temp_co2)
            self.tau_above_aerosol_star_band_co2.append(tau_above_aerosol_star_temp_co2)
            self.tau_star_band_ch4.append(tau_star_temp_ch4)
            self.tau_above_aerosol_star_band_ch4.append(tau_above_aerosol_star_temp_ch4)


        #####
        #Calculate radiances for each band
        self.R_band = []
        self.R_band_albedo = []
        self.R_band_aerosol = []
        self.R_band_q = []
        self.R_band_co2 = []
        self.R_band_ch4 = []

        #print("Calculating radiances...")
        for i in range(len(self.band_min_um)):

            if self.tau_aerosol != None: tau_aerosol_temp = np.full(len(self.band_absco_res_wn[i]),self.tau_aerosol)
            else: 
              tau_aerosol_temp = np.zeros(len(self.band_absco_res_wn[i]))

            I, I_albedo, I_aerosol, I_q, I_co2, I_ch4 = self.intensity(self.band_absco_res_wn[i],self.tau_star_band[i],self.tau_above_aerosol_star_band[i],self.tau_star_band_q[i],self.tau_above_aerosol_star_band_q[i],self.tau_star_band_co2[i],self.tau_above_aerosol_star_band_co2[i],self.tau_star_band_ch4[i],self.tau_above_aerosol_star_band_ch4[i],tau_aerosol_temp,self.ssa_aerosol[i],self.P_aerosol[i],self.qext_aerosol[0],self.qext_aerosol[i],self.mu,self.mu_0,self.m,self.albedo[i],self.band_solar_irradiances[i],self.jacobians)
            
            #Calculate the spectral response function (with and without multiplying by intensity)
            Sc_I_band, Sc_I_band_albedo, Sc_I_band_aerosol, Sc_I_band_q, Sc_I_band_co2, Sc_I_band_ch4 = self.spectral_response_function(self.band_wn_index[i],self.band_absco_res_wn[i],self.sigma_band[i],self.ILS_Gaussian_term[i],I,I_albedo,I_aerosol,I_q,I_co2,I_ch4,self.jacobians)

            #Calculate radiance (Rc) by integrating intensity times ILS, and reverse to plot in micrometers
            Rc_band = (Sc_I_band/self.ILS_Gaussian_term_sum[i])[::-1]

            #For analytic Jacobians
            if jacobians:
              Rc_band_albedo = (Sc_I_band_albedo/self.ILS_Gaussian_term_sum[i])[::-1]
              Rc_band_aerosol = (Sc_I_band_aerosol/self.ILS_Gaussian_term_sum[i])[::-1]
              Rc_band_q = (Sc_I_band_q/self.ILS_Gaussian_term_sum[i][:,None])[::-1,:]
              Rc_band_co2 = (Sc_I_band_co2/self.ILS_Gaussian_term_sum[i][:,None])[::-1,:]
              Rc_band_ch4 = (Sc_I_band_ch4/self.ILS_Gaussian_term_sum[i][:,None])[::-1,:]

            #Append for the band we're on
            self.R_band.append(Rc_band)
            if jacobians:
              #For analytic Jacobians
              self.R_band_albedo.append(Rc_band_albedo)
              self.R_band_aerosol.append(Rc_band_aerosol)
              self.R_band_q.append(Rc_band_q)
              self.R_band_co2.append(Rc_band_co2)
              self.R_band_ch4.append(Rc_band_ch4)

        self.y = np.concatenate(self.R_band)
        if jacobians:
          self.y_albedo = np.concatenate(self.R_band_albedo)
          self.y_aerosol = np.concatenate(self.R_band_aerosol)
          self.y_q = np.concatenate(self.R_band_q)
          self.y_co2 = np.concatenate(self.R_band_co2)
          self.y_ch4 = np.concatenate(self.R_band_ch4)

        noise = []
        for i in range(len(self.band_max_wn)):
          signal = self.R_band[i].max()
          sigma = signal/self.SNR
          noise_temp = np.random.normal(0,sigma,self.band_spectral_points[i])
          noise.append(noise_temp)

          #If we're adding noise:
          if self.measurement_error: self.R_band[i] = self.R_band[i] + noise_temp

        #Now combine into one spectra
        self.y = np.concatenate(self.R_band)

        #Calculate Sy even if we didn't add noise because we need something for Sy
        noise_std = []
        for i in range(len(self.band_max_wn)):
          noise_std.append(np.full((len(noise[i])),np.std(noise[i])**2.0))
        Sy = np.diag(np.concatenate(noise_std))

        #Calculate this ahead of time so we don't have to calculate it every retrieval iteration
        Sy_inv = np.zeros(Sy.shape)
        np.fill_diagonal(Sy_inv,1./Sy.diagonal())

        self.Sy_inv = Sy_inv


    #Calculate intensities for a single band
    def intensity(self,band,tau_star_band,tau_above_aerosol_star_band,tau_star_band_q,tau_above_aerosol_star_band_q,tau_star_band_co2,tau_above_aerosol_star_band_co2,tau_star_band_ch4,tau_above_aerosol_star_band_ch4,tau_aerosol,ssa_aerosol,P_aerosol,qext_aerosol_band_0,qext_aerosol,mu,mu_0,m,albedo,band_solar_irradiances,jacobians):

      I = np.zeros((len(band))) #wn
      I_albedo = np.zeros((len(band))) #wn x layers
      I_aerosol = np.zeros((len(band))) #wn x layers
      I_q = np.zeros((len(band),tau_star_band_q.shape[1])) #wn x layers
      I_co2 = np.zeros((len(band),tau_star_band_q.shape[1])) #wn x layers
      I_ch4 = np.zeros((len(band),tau_star_band_q.shape[1])) #wn x layers

      #Dealing with divide by zero issues
      if qext_aerosol_band_0[0] == 0:
        qext_scaling = np.zeros((len(qext_aerosol)))
      else:
        qext_scaling = qext_aerosol/qext_aerosol_band_0[0]

      #Direct exponential term
      exp_term = np.exp(-m*(tau_star_band + tau_aerosol*qext_scaling))

      #Scattering exponential term
      exp_term_above_aerosol = np.exp(-m*tau_above_aerosol_star_band)

      for i in range(len(band)):
        #Add an aerosol layer. Assume it scatters once.
        #Full qext scaling
        I[i] = band_solar_irradiances/np.pi * (albedo*mu_0*exp_term[i] + ssa_aerosol[i]*P_aerosol[i]*tau_aerosol[i]*qext_scaling[i]*exp_term_above_aerosol[i]/4./mu)

        #Calculate analytical Jacobians
        if jacobians:
          I_albedo[i] = band_solar_irradiances/np.pi * mu_0 * exp_term[i]

          I_aerosol[i] = band_solar_irradiances/np.pi * (-m*qext_scaling[i]*albedo*mu_0*exp_term[i] + ssa_aerosol[i]*P_aerosol[i]*qext_scaling[i]*exp_term_above_aerosol[i]/4./mu)

          #Full qext scaling:
          I_q[i,:] = band_solar_irradiances/np.pi * (albedo*mu_0*exp_term[i] * (-m) * tau_star_band_q[i,:] + ssa_aerosol[i]*P_aerosol[i]*tau_aerosol[i]*qext_scaling[i]*exp_term_above_aerosol[i]/4./mu * (-m) * tau_above_aerosol_star_band_q[i,:])
          I_co2[i,:] = band_solar_irradiances/np.pi * (-m*(albedo*mu_0*exp_term[i] * tau_star_band_co2[i,:] + ssa_aerosol[i]*P_aerosol[i]*tau_aerosol[i]*qext_scaling[i]*exp_term_above_aerosol[i]/4./mu * tau_above_aerosol_star_band_co2[i,:]))
          I_ch4[i,:] = band_solar_irradiances/np.pi * (-m*(albedo*mu_0*exp_term[i] * tau_star_band_ch4[i,:] + ssa_aerosol[i]*P_aerosol[i]*tau_aerosol[i]*qext_scaling[i]*exp_term_above_aerosol[i]/4./mu * tau_above_aerosol_star_band_ch4[i,:]))

      return I, I_albedo, I_aerosol, I_q, I_co2, I_ch4


    #Assume a Gaussian ILS
    def spectral_response_function(self,band_wn_index,band,sigma_band,ILS_Gaussian_term,I_band,I_band_albedo,I_band_aerosol,I_band_q,I_band_co2,I_band_ch4,jacobians):

        Sc_I_band = np.zeros((len(band_wn_index))) #wn instrument

        Sc_I_band_albedo = np.zeros((len(band_wn_index))) #wn instrument
        Sc_I_band_aerosol = np.zeros((len(band_wn_index))) #wn instrument
        Sc_I_band_q = np.zeros((len(band_wn_index),I_band_q.shape[1])) #wn instrument x layers
        Sc_I_band_co2 = np.zeros((len(band_wn_index),I_band_co2.shape[1])) #wn instrument x layers
        Sc_I_band_ch4 = np.zeros((len(band_wn_index),I_band_ch4.shape[1])) #wn instrument x layers

        round_term = round(sigma_band*ILS_width*100.0)

        for i in range(len(band_wn_index)):

            #Dealing with the starting edge of the band
            if band[band_wn_index[i]] <= band[0]+sigma_band*ILS_width:

                j_index_temp_lower = 0
                j_index_temp_upper = int(band_wn_index[i]+round_term)

            #Dealing with the trailing edge of the band
            elif band[band_wn_index[i]] >= band[len(band)-1]-sigma_band*ILS_width:

                j_index_temp_lower = int(band_wn_index[i]-round_term)
                j_index_temp_upper = len(band)

            #Most of the band
            else:
                j_index_temp_lower = int(band_wn_index[i]-round_term)
                j_index_temp_upper = int(band_wn_index[i]+round_term)

            Sc_I_band[i] = np.sum(I_band[j_index_temp_lower:j_index_temp_upper] * ILS_Gaussian_term[i,j_index_temp_lower:j_index_temp_upper])

            if jacobians:
                Sc_I_band_albedo[i] = np.sum(I_band_albedo[j_index_temp_lower:j_index_temp_upper] * ILS_Gaussian_term[i,j_index_temp_lower:j_index_temp_upper])
                Sc_I_band_aerosol[i] = np.sum(I_band_aerosol[j_index_temp_lower:j_index_temp_upper] * ILS_Gaussian_term[i,j_index_temp_lower:j_index_temp_upper])
                Sc_I_band_q[i,:] = np.sum(I_band_q[j_index_temp_lower:j_index_temp_upper,:] * ILS_Gaussian_term[i,j_index_temp_lower:j_index_temp_upper,None],axis=0)
                Sc_I_band_co2[i,:] = np.sum(I_band_co2[j_index_temp_lower:j_index_temp_upper,:] * ILS_Gaussian_term[i,j_index_temp_lower:j_index_temp_upper,None],axis=0)
                Sc_I_band_ch4[i,:] = np.sum(I_band_ch4[j_index_temp_lower:j_index_temp_upper,:] * ILS_Gaussian_term[i,j_index_temp_lower:j_index_temp_upper,None],axis=0)

        return Sc_I_band, Sc_I_band_albedo, Sc_I_band_aerosol, Sc_I_band_q, Sc_I_band_co2, Sc_I_band_ch4


class Retrieval:
    '''Retrieval object. Parameter defaults are set in settings.py
    band_max_wn               [array] sfddsfds
    band_min_wn               [array] sdfdsf
    band_spectral_resolution  [array] fdsf

    '''

    def __init__(self):

        self.iterations = 0
        self.chisq_reduced_previous = 9999999.0

    def run(self, x, model_prior, model_true, absco_data, chisq_threshold=s.chisq_threshold):

        time_total=time.time()

        #Set the initial guess to the prior
        x["ret"] = deepcopy(x["prior"])

        self.chisq_threshold = chisq_threshold
        done=False
        while not done:

            #Time the iterations
            time1=time.time()

            if self.iterations == 0:
                print("-----------------------------------------------------------------------------------------------------------------------")
                print("State vector:")
                print("Name                                              Prior  True")
                for i in range(len(x["prior"])):
                    print(x["names"][i].ljust(49)+" "+str(x["prior"][i]).ljust(5)+"  "+str(x["true"][i]).ljust(20))
                print("-----------")

            #Calculate y using the state vector, x
            model_ret = self.forward_model(x, model_prior, absco_data, jacobians=True)

            #Calculate chisq
            chisq = (((model_true.y-model_ret.y).dot(model_true.Sy_inv)).dot((model_true.y-model_ret.y).T))+\
                (((x["ret"]-x["prior"]).dot(LA.inv(x["S_prior"]))).dot((x["ret"]-x["prior"]).T))
            self.chisq_reduced = chisq/len(model_true.y)

            #Check if our chisq_reduced is good enough to stop
            print("Reduced chisq = ",self.chisq_reduced)
            if self.chisq_reduced < self.chisq_threshold: 
                print("Reduced chisq is less than "+str(self.chisq_threshold)+", so we're done!")
                done=True
                continue

            #Also stop if chisq is barely changing
            elif (self.chisq_reduced_previous - self.chisq_reduced) < 0.0001:
                print("Reduced chisq has stopped changing or increased, so we're done!")
                done=True
                continue

            else: print("Reduced chisq is not small enough yet, so update the state vector")
            self.chisq_reduced_previous = self.chisq_reduced

            #If not, loop through the state vector to calculate the Jacobian (K)
            self.K = np.zeros((len(model_true.y),len(x["ret"])))
            for i in range(len(x["ret"])):

                #Need to use finite differencing for T and p Jacobians
                if ("Temperature" in x["names"][i]) or ("Pressure" in x["names"][i]):
                    x_perturbed = deepcopy(x)
                    x_perturbed["ret"][i] += s.perturbation
                    model_perturbed = self.forward_model(x_perturbed, model_prior, absco_data, jacobians=False)
                    self.K[:,i] = ((model_perturbed.y - model_ret.y)/s.perturbation)

                #If we have the analytical derivative for the current state, use it!
                else: self.K[:,i] = model_ret.y_k[:,i]

            #Calculate the current error covariance matrix (S)
            self.S = LA.inv((self.K.T).dot(model_true.Sy_inv).dot(self.K) + LA.inv(x["S_prior"]))

            #If we haven't converged, determine how we want to perturb the state vector to try again
            x["dx"] = self.S.dot((((self.K.T).dot(model_true.Sy_inv).dot(model_true.y - model_ret.y)) - (LA.inv(x["S_prior"]).dot(x["ret"]-x["prior"]))))

            print("Updated state vector (iteration "+str(self.iterations+1)+"):")
            print("Name                                             Prior  True                  Current               Dx from prev. iteration:")
            for i in range(len(x["ret"])):
                print(x["names"][i].ljust(49)+str(x["prior"][i]).ljust(5)+"  "+str(x["true"][i]).ljust(20)+"  "+str(x["ret"][i] + x["dx"][i]).ljust(20)+"  "+str(x["dx"][i]))

            #Modify the state vector
            x["ret"] += x["dx"]

            #Print things. Index of CO2/CH4 depend on the forward model used.
            if "CO2 Profile Scale Factor" in x["names"] and "CH4 Profile Scale Factor" in x["names"]:
                print("-----------")
                print("Prior XCO2 =".ljust(31),'{:.5f}'.format(model_prior.xco2).rjust(10),"ppm")
                print("True XCO2 =".ljust(31),'{:.5f}'.format(model_true.xco2).rjust(10),"ppm")
                print("Current retrieved XCO2 =".ljust(31),'{:.5f}'.format(calculate_Xgas(model_prior.co2*x["ret"][0], model_prior.p*x["ret"][3], model_prior.q*x["ret"][4])[0] * 1e6).rjust(10),"ppm")
                print("XCO2 error (retrieved - true) =".ljust(31),'{:.5f}'.format(calculate_Xgas(model_prior.co2*x["ret"][0], model_prior.p*x["ret"][3], model_prior.q*x["ret"][4])[0] * 1e6 - model_true.xco2).rjust(10),"ppm")

                print("-----------")
                print("Prior XCH4 =".ljust(31),'{:.5f}'.format(model_prior.xch4).rjust(10),"ppb")
                print("True XCH4 =".ljust(31),'{:.5f}'.format(model_true.xch4).rjust(10),"ppb")
                print("Current retrieved XCH4 =".ljust(31),'{:.5f}'.format(calculate_Xgas(model_prior.ch4*x["ret"][1], model_prior.p*x["ret"][3], model_prior.q*x["ret"][4])[0] * 1e9).rjust(10),"ppb")
                print("XCH4 error (retrieved - true) =".ljust(31),'{:.5f}'.format(calculate_Xgas(model_prior.ch4*x["ret"][1], model_prior.p*x["ret"][3], model_prior.q*x["ret"][4])[0] * 1e9 - model_true.xch4).rjust(10),"ppb")
                print("-----------")

            elif "CO2 Profile Scale Factor" in x["names"] and "CH4 Profile Scale Factor" not in x["names"]:
                print("-----------")
                print("Prior XCO2 =".ljust(31),'{:.5f}'.format(model_prior.xco2).rjust(10),"ppm")
                print("True XCO2 =".ljust(31),'{:.5f}'.format(model_true.xco2).rjust(10),"ppm")
                print("Current retrieved XCO2 =".ljust(31),'{:.5f}'.format(calculate_Xgas(model_prior.co2*x["ret"][0], model_prior.p*x["ret"][2], model_prior.q*x["ret"][3])[0] * 1e6).rjust(10),"ppm")
                print("XCO2 error (retrieved - true) =".ljust(31),'{:.5f}'.format(calculate_Xgas(model_prior.co2*x["ret"][0], model_prior.p*x["ret"][2], model_prior.q*x["ret"][3])[0] * 1e6 - model_true.xco2).rjust(10),"ppm")
                print("-----------")

            elif "CO2 Profile Scale Factor" not in x["names"] and "CH4 Profile Scale Factor" in x["names"]:
                print("-----------")
                print("Prior XCH4 =".ljust(31),'{:.5f}'.format(model_prior.xch4).rjust(10),"ppb")
                print("True XCH4 =".ljust(31),'{:.5f}'.format(model_true.xch4).rjust(10),"ppb")
                print("Current retrieved XCH4 =".ljust(31),'{:.5f}'.format(calculate_Xgas(model_prior.ch4*x["ret"][0], model_prior.p*x["ret"][2], model_prior.q*x["ret"][3])[0] * 1e9).rjust(10),"ppb")
                print("XCH4 error (retrieved - true) =".ljust(31),'{:.5f}'.format(calculate_Xgas(model_prior.ch4*x["ret"][0], model_prior.p*x["ret"][2], model_prior.q*x["ret"][3])[0] * 1e9 - model_true.xch4).rjust(10),"ppb")
                print("-----------")

            else:
                print("Unexpected state vector setup!")

            print("Time for iteration",self.iterations+1,"=",'{:.2f}'.format(time.time()-time1), "s")
            print("-----------------------------------------------------------------------------------------------------------------------")

            #Add iteration
            self.iterations += 1

            if self.iterations > 5:
                print("More than 5 iterations, so we're stopping...")
                done=True
                continue

        print("-----------------------------------------------------------------------------------------------------------------------")
        print("Final reduced chisq = ",self.chisq_reduced)
        if "CO2 Profile Scale Factor" in x["names"] and "CH4 Profile Scale Factor" in x["names"]:
            xco2_ret_temp = calculate_Xgas(model_prior.co2*x["ret"][0], model_prior.p*x["ret"][3], model_prior.q*x["ret"][4])[0] * 1e6
            xch4_ret_temp = calculate_Xgas(model_prior.ch4*x["ret"][1], model_prior.p*x["ret"][3], model_prior.q*x["ret"][4])[0] * 1e9
            print("Final retrieved XCO2 =".ljust(23),'{:.5f}'.format(xco2_ret_temp).rjust(10),"+/-",'{:.5f}'.format(xco2_ret_temp*np.diagonal(self.S)[0]**0.5),"ppm")
            print("Final XCO2 error (retrieved - true) =".ljust(37),'{:.5f}'.format(xco2_ret_temp - model_true.xco2).rjust(8),"ppm")
            print("Final retrieved XCH4 =".ljust(23),'{:.5f}'.format(xch4_ret_temp).rjust(10),"+/-",'{:.5f}'.format(xch4_ret_temp*np.diagonal(self.S)[1]**0.5),"ppb")
            print("Final XCH4 error (retrieved - true) =".ljust(37),'{:.5f}'.format(xch4_ret_temp - model_true.xch4).rjust(8),"ppb")

        elif "CO2 Profile Scale Factor" in x["names"] and "CH4 Profile Scale Factor" not in x["names"]:
            xco2_ret_temp = calculate_Xgas(model_prior.co2*x["ret"][0], model_prior.p*x["ret"][2], model_prior.q*x["ret"][3])[0] * 1e6
            print("Final retrieved XCO2 =".ljust(23),'{:.5f}'.format(xco2_ret_temp).rjust(10),"+/-",'{:.5f}'.format(xco2_ret_temp*np.diagonal(self.S)[0]**0.5),"ppm")
            print("Final XCO2 error (retrieved - true) =".ljust(37),'{:.5f}'.format(xco2_ret_temp - model_true.xco2).rjust(8),"ppm")

        elif "CO2 Profile Scale Factor" not in x["names"] and "CH4 Profile Scale Factor" in x["names"]:
            xch4_ret_temp = calculate_Xgas(model_prior.ch4*x["ret"][0], model_prior.p*x["ret"][2], model_prior.q*x["ret"][3])[0] * 1e9
            print("Final retrieved XCH4 =".ljust(23),'{:.5f}'.format(xch4_ret_temp).rjust(10),"+/-",'{:.5f}'.format(xch4_ret_temp*np.diagonal(self.S)[0]**0.5),"ppb")
            print("Final XCH4 error (retrieved - true) =".ljust(37),'{:.5f}'.format(xch4_ret_temp - model_true.xch4).rjust(8),"ppb")

        else:
            print("Unexpected state vector setup!")

        print("Total retrieval time =",'{:.2f}'.format(time.time()-time_total),"s")


    def forward_model(self, x, model_prior, absco_data, jacobians=False):

        #Modify the prior state vector appropriately
        #Full-physics setup:
        if "CO2 Profile Scale Factor" in x["names"] and "CH4 Profile Scale Factor" in x["names"]:
            co2 = model_prior.co2 * x["ret"][0]
            ch4 = model_prior.ch4 * x["ret"][1]
            T = model_prior.T + x["ret"][2]
            p = model_prior.p * x["ret"][3]
            q = model_prior.q * x["ret"][4]
            albedo_band_1 = x["ret"][5]
            albedo_band_2 = x["ret"][6]
            albedo_band_3 = x["ret"][7]
            albedo = [albedo_band_1, albedo_band_2, albedo_band_3]

            #Set tau_aerosol appropriately
            if "Aerosol Optical Depth" in x["names"]: tau_aerosol = x["ret"][8]
            else: tau_aerosol = None

            #Call the foward function with info from the prior and the updated state vector elements
            model = ForwardFunction(SNR=model_prior.SNR,sza_0=model_prior.sza_0,sza=model_prior.sza,co2=co2,ch4=ch4,T=T,p=p,q=q,albedo=albedo,band_min_wn=model_prior.band_min_wn,band_max_wn=model_prior.band_max_wn,band_spectral_resolutions=model_prior.band_spectral_resolutions,band_min_um=model_prior.band_min_um,band_max_um=model_prior.band_max_um,band_spectral_points=model_prior.band_spectral_points,band_wn=model_prior.band_wn,band_wl=model_prior.band_wl,band_absco_res_wn=model_prior.band_absco_res_wn,resolving_power_band=model_prior.resolving_power_band,sigma_band=model_prior.sigma_band,band_wn_index=model_prior.band_wn_index,ILS_Gaussian_term=model_prior.ILS_Gaussian_term,ILS_Gaussian_term_sum=model_prior.ILS_Gaussian_term_sum,absco_data=absco_data,band_molecules=model_prior.band_molecules,P_aerosol=model_prior.P_aerosol,ssa_aerosol=model_prior.ssa_aerosol,qext_aerosol=model_prior.qext_aerosol,height_aerosol=model_prior.height_aerosol,tau_aerosol=tau_aerosol,jacobians=jacobians)

            #Calculate analytical derivative
            model.y_k = np.zeros((len(model.y),len(x["ret"])))

            if jacobians:
                model.y_k[:,0] = np.sum(model_prior.co2_layer[None,:] * model.y_co2, axis=1) #co2
                model.y_k[:,1] = np.sum(model_prior.ch4_layer[None,:] * model.y_ch4, axis=1) #ch4
                #model.y_k[:,2] #Calculate using finite differencing #T
                #model.y_k[:,3] #Calculate using finite differencing #p
                model.y_k[:,4] = np.sum(model_prior.q_layer[None,:] * model.y_q, axis=1) #q
                model.y_k[:,5] = np.concatenate((model.R_band_albedo[0],np.zeros((len(model.R_band_albedo[1]))),np.zeros((len(model.R_band_albedo[2]))))) #Band 1 albedo
                model.y_k[:,6] = np.concatenate((np.zeros((len(model.R_band_albedo[0]))),model.R_band_albedo[1],np.zeros((len(model.R_band_albedo[2]))))) #Band 2 albedo
                model.y_k[:,7] = np.concatenate((np.zeros((len(model.R_band_albedo[0]))),np.zeros((len(model.R_band_albedo[1]))),model.R_band_albedo[2])) #Band 3 albedo
                if "Aerosol Optical Depth" in x["names"]: model.y_k[:,8] = model.y_aerosol #tau_aerosol

        #CO2-only run
        elif "CO2 Profile Scale Factor" in x["names"] and "CH4 Profile Scale Factor" not in x["names"]:
            co2 = model_prior.co2 * x["ret"][0]
            T = model_prior.T + x["ret"][1]
            p = model_prior.p * x["ret"][2]
            q = model_prior.q * x["ret"][3]
            albedo_band_2 = x["ret"][4]
            albedo = [albedo_band_2]
            tau_aerosol = None

            #Call the foward function with info from the prior and the updated state vector elements
            model = ForwardFunction(SNR=model_prior.SNR,sza_0=model_prior.sza_0,sza=model_prior.sza,co2=co2,ch4=np.zeros(len(model_prior.ch4)),T=T,p=p,q=q,albedo=albedo,band_min_wn=model_prior.band_min_wn,band_max_wn=model_prior.band_max_wn,band_spectral_resolutions=model_prior.band_spectral_resolutions,band_min_um=model_prior.band_min_um,band_max_um=model_prior.band_max_um,band_spectral_points=model_prior.band_spectral_points,band_wn=model_prior.band_wn,band_wl=model_prior.band_wl,band_absco_res_wn=model_prior.band_absco_res_wn,resolving_power_band=model_prior.resolving_power_band,sigma_band=model_prior.sigma_band,band_wn_index=model_prior.band_wn_index,ILS_Gaussian_term=model_prior.ILS_Gaussian_term,ILS_Gaussian_term_sum=model_prior.ILS_Gaussian_term_sum,absco_data=absco_data,band_molecules=model_prior.band_molecules,P_aerosol=model_prior.P_aerosol,ssa_aerosol=model_prior.ssa_aerosol,qext_aerosol=model_prior.qext_aerosol,height_aerosol=model_prior.height_aerosol,tau_aerosol=tau_aerosol,jacobians=jacobians)

            #Calculate analytical derivative
            model.y_k = np.zeros((len(model.y),len(x["ret"])))

            if jacobians:
                model.y_k[:,0] = np.sum(model_prior.co2_layer[None,:] * model.y_co2, axis=1) #co2
                #model.y_k[:,1] #Calculate using finite differencing #T
                #model.y_k[:,2] #Calculate using finite differencing #p
                model.y_k[:,3] = np.sum(model_prior.q_layer[None,:] * model.y_q, axis=1) #q
                model.y_k[:,4] = model.R_band_albedo[0] #Band 2 albedo. 0th index in this case.

        #CH4-only run
        elif "CO2 Profile Scale Factor" not in x["names"] and "CH4 Profile Scale Factor" in x["names"]:
            ch4 = model_prior.ch4 * x["ret"][0]
            T = model_prior.T + x["ret"][1]
            p = model_prior.p * x["ret"][2]
            q = model_prior.q * x["ret"][3]
            albedo_band_3 = x["ret"][4]
            albedo = [albedo_band_3]
            tau_aerosol = None

            #Call the foward function with info from the prior and the updated state vector elements
            model = ForwardFunction(SNR=model_prior.SNR,sza_0=model_prior.sza_0,sza=model_prior.sza,co2=np.zeros(len(model_prior.co2)),ch4=ch4,T=T,p=p,q=q,albedo=albedo,band_min_wn=model_prior.band_min_wn,band_max_wn=model_prior.band_max_wn,band_spectral_resolutions=model_prior.band_spectral_resolutions,band_min_um=model_prior.band_min_um,band_max_um=model_prior.band_max_um,band_spectral_points=model_prior.band_spectral_points,band_wn=model_prior.band_wn,band_wl=model_prior.band_wl,band_absco_res_wn=model_prior.band_absco_res_wn,resolving_power_band=model_prior.resolving_power_band,sigma_band=model_prior.sigma_band,band_wn_index=model_prior.band_wn_index,ILS_Gaussian_term=model_prior.ILS_Gaussian_term,ILS_Gaussian_term_sum=model_prior.ILS_Gaussian_term_sum,absco_data=absco_data,band_molecules=model_prior.band_molecules,P_aerosol=model_prior.P_aerosol,ssa_aerosol=model_prior.ssa_aerosol,qext_aerosol=model_prior.qext_aerosol,height_aerosol=model_prior.height_aerosol,tau_aerosol=tau_aerosol,jacobians=jacobians)

            #Calculate analytical derivative
            model.y_k = np.zeros((len(model.y),len(x["ret"])))

            if jacobians:
                model.y_k[:,0] = np.sum(model_prior.ch4_layer[None,:] * model.y_ch4, axis=1) #ch4
                #model.y_k[:,1] #Calculate using finite differencing #T
                #model.y_k[:,2] #Calculate using finite differencing #p
                model.y_k[:,3] = np.sum(model_prior.q_layer[None,:] * model.y_q, axis=1) #q
                model.y_k[:,4] = model.R_band_albedo[0] #Band 3 albedo. 0th index in this case.

        else: 
          print("Unexpected state vector setup!")

        return model


#Function to calculate column-mean dry-air mole-fraction of a gas. Assume gravity is constant with height.
def calculate_Xgas(gas_level,p_level,q_level):

    #Calculate the d(pressure) of a given layer
    p_diff=np.empty((len(p_level)-1))
    for i in range(len(p_level)-1):
        p_diff[i] = p_level[i+1]-p_level[i]

    gas_layer=np.empty((len(gas_level)-1)) #Layer gas concentrations
    q_layer=np.empty((len(q_level)-1)) #Layer specific humidities
    for i in range(len(q_level)-1):
        gas_layer[i] = (gas_level[i+1]+gas_level[i])/2.0
        q_layer[i] = (q_level[i+1]+q_level[i])/2.0

    #Calculate the column density of dry air per unit pressure
    c_layer = (1.-q_layer)/g/M

    #Calcualte the pressure weighting function (PWF)
    h = (c_layer * p_diff) / np.sum(c_layer * p_diff)

    #Calculate the column-average dry-air mole-fraction of the gas
    xgas = np.sum(h * gas_layer)

    return xgas, h


#Find the layer average for a given array of levels
def layer_average(levels):
    layers=np.empty((len(levels)-1)) #Layer pressures
    for i in range(len(levels)-1):
        layers[i] = (levels[i+1]+levels[i])/2.0
    return layers


#Calculate the solar flux for each band's average wavelength. Ignoring all other complexities with the solar spectrum for now.
def planck(T_temp,wl_temp):
    B_temp=2.0*h*c**2.0/wl_temp**5.0/(np.exp(h*c/k/wl_temp/T_temp)-1.0)
    return B_temp


#Calculate the number of dry air molecules
def calculate_n_dry_air(p_diff,q_layer):
  n_dry_air = p_diff * Na / M / g * (1.-q_layer)
  return n_dry_air


