import h5py
import numpy as np
import sys
import os
import glob

Na=6.022e23 #molecule/mol^-1

#Function to load the v5.2 ABSCO tables, which are used in the OCO-2/3 B11 L2
def load_absco(data_path):

  absco_h2o_file_name = glob.glob(data_path+"h2o*.hdf")
  if len(absco_h2o_file_name) == 1 and os.path.isfile(absco_h2o_file_name[0]):
    absco_h2o_file = h5py.File(absco_h2o_file_name[0],'r')
  else:
    print("Couldn't find the H2O ABSCO file!")
    return
  absco_h2o_dataset = np.array(absco_h2o_file["Gas_01_Absorption"][:,:,0,:])
  absco_h2o_wavenumber = np.array(absco_h2o_file["Wavenumber"])
  absco_h2o_p=np.array((absco_h2o_file['Pressure']))
  absco_h2o_T=np.array((absco_h2o_file['Temperature']))
  absco_h2o_file.close()

  absco_co2_file_name = glob.glob(data_path+"co2*.hdf")
  if len(absco_co2_file_name) == 1 and os.path.isfile(absco_co2_file_name[0]):
    absco_co2_file = h5py.File(absco_co2_file_name[0],'r')
  else:
    print("Couldn't find the CO2 ABSCO file!")
    return
  absco_co2_dataset = np.array(absco_co2_file["Gas_02_Absorption"][:,:,0,:])
  absco_co2_wavenumber = np.array(absco_co2_file["Wavenumber"])
  absco_co2_p=np.array((absco_co2_file['Pressure']))
  absco_co2_T=np.array((absco_co2_file['Temperature']))
  absco_co2_file.close()

  absco_o2_file_name = glob.glob(data_path+"o2*.hdf")
  if len(absco_o2_file_name) == 1 and os.path.isfile(absco_o2_file_name[0]):
    absco_o2_file = h5py.File(absco_o2_file_name[0],'r')
  else:
    print("Couldn't find the O2 ABSCO file!")
    return
  absco_o2_dataset = np.array(absco_o2_file["Gas_07_Absorption"][:,:,0,:])
  absco_o2_wavenumber = np.array(absco_o2_file["Wavenumber"])
  absco_o2_p=np.array((absco_o2_file['Pressure']))
  absco_o2_T=np.array((absco_o2_file['Temperature']))
  absco_o2_file.close()

  absco_ch4_file_name = glob.glob(data_path+"ch4*.hdf")
  if len(absco_ch4_file_name) == 1 and os.path.isfile(absco_ch4_file_name[0]):
    absco_ch4_file = h5py.File(absco_ch4_file_name[0],'r')
  else:
    print("Couldn't find the CH4 ABSCO file!")
    return
  absco_ch4_dataset = np.array(absco_ch4_file["Gas_06_Absorption"][:,:,0,:])
  absco_ch4_wavenumber = np.array(absco_ch4_file["Wavenumber"])
  absco_ch4_p=np.array((absco_ch4_file['Pressure']))
  absco_ch4_T=np.array((absco_ch4_file['Temperature']))
  absco_ch4_file.close()

  #Create a dictionary with all the ABSCO data we need
  absco_dist = {'absco_h2o_dataset':absco_h2o_dataset,
                'absco_h2o_wavenumber':absco_h2o_wavenumber,
                'absco_h2o_p':absco_h2o_p,
                'absco_h2o_T':absco_h2o_T,
                'absco_co2_dataset':absco_co2_dataset,
                'absco_co2_wavenumber':absco_co2_wavenumber,
                'absco_co2_p':absco_co2_p,
                'absco_co2_T':absco_co2_T,
                'absco_o2_dataset':absco_o2_dataset,
                'absco_o2_wavenumber':absco_o2_wavenumber,
                'absco_o2_p':absco_o2_p,
                'absco_o2_T':absco_o2_T,
                'absco_ch4_dataset':absco_ch4_dataset,
                'absco_ch4_wavenumber':absco_ch4_wavenumber,
                'absco_ch4_p':absco_ch4_p,
                'absco_ch4_T':absco_ch4_T}

  return absco_dist

def load_absco_file(absco_file):
          
  absco_data = h5py.File(absco_file,"r")

  absco_dist = {'absco_h2o_dataset':np.array(absco_data['absco_h2o_dataset']),
                'absco_h2o_wavenumber':np.array(absco_data['absco_h2o_wavenumber']),
                'absco_h2o_p':np.array(absco_data['absco_h2o_p']),
                'absco_h2o_T':np.array(absco_data['absco_h2o_T']),
                'absco_co2_dataset':np.array(absco_data['absco_co2_dataset']),
                'absco_co2_wavenumber':np.array(absco_data['absco_co2_wavenumber']),
                'absco_co2_p':np.array(absco_data['absco_co2_p']),
                'absco_co2_T':np.array(absco_data['absco_co2_T']),
                'absco_o2_dataset':np.array(absco_data['absco_o2_dataset']),
                'absco_o2_wavenumber':np.array(absco_data['absco_o2_wavenumber']),
                'absco_o2_p':np.array(absco_data['absco_o2_p']),
                'absco_o2_T':np.array(absco_data['absco_o2_T']),
                'absco_ch4_dataset':np.array(absco_data['absco_ch4_dataset']),
                'absco_ch4_wavenumber':np.array(absco_data['absco_ch4_wavenumber']),
                'absco_ch4_p':np.array(absco_data['absco_ch4_p']),
                'absco_ch4_T':np.array(absco_data['absco_ch4_T'])}

  return absco_dist

#def load_absco_file_npy(absco_file):
#  #Order matters here!
#  with open(absco_file, 'rb') as f:
#    absco_h2o_dataset = np.load(f)
#    absco_h2o_wavenumber = np.load(f)
#    absco_h2o_p = np.load(f)
#    absco_h2o_T = np.load(f)
#    absco_co2_dataset = np.load(f)
#    absco_co2_wavenumber = np.load(f)
#    absco_co2_p = np.load(f)
#    absco_co2_T = np.load(f)
#    absco_o2_dataset = np.load(f)
#    absco_o2_wavenumber = np.load(f)
#    absco_o2_p = np.load(f)
#    absco_o2_T = np.load(f)
#    absco_ch4_dataset = np.load(f)
#    absco_ch4_wavenumber = np.load(f)
#    absco_ch4_p = np.load(f)
#    absco_ch4_T = np.load(f)


def interpol1(x, A, bound):
   # find the position of a single x within the A array.
   # A must be sorted (ascending)
   # Returns a single real number that shows the position of x within A.
   #    if BOUND is set, itarg will be between 1 and n_elements(A)

   # this function is exactly like interpol but works on a single x instead of an array.

   p = 0
   nA = len(A)

   while A[p]<x and p<nA-1:
     p+=1
   p-=1

   if (p<=0): p=0

   if x==A[p]:
     itarg=p
   else: itarg = p + (x - A[p])/(A[p+1]-A[p])

   if bound:
     itarg = min(max(itarg,0.0),nA)

   return itarg


def interpol_lookup(p,p_axis,T,T_axis):
  fp = interpol1(p,p_axis,False)
  ip = max(min(int(fp),len(p_axis)-2),0)
  fp = fp - ip

  ft1 = interpol1(T,T_axis[ip,:],False)
  it1 = max(min(int(ft1),len(T_axis[0,:])-2),0)
  ft1 = ft1 - it1

  ft2 = interpol1(T,T_axis[ip+1],False)
  it2 = max(min(int(ft2),len(T_axis[0,:])-2),0)
  ft2 = ft2 - it2

  return ip,it1,it2,fp,ft1,ft2

#Function to interpolate the ABSCO table along the temperature and pressure axes
def sigma_lookup(molecule,band_absco_res,p_ave,T_ave,absco_data):
  ip=np.empty((len(p_ave)))
  it1=np.empty((len(p_ave)))
  it2=np.empty((len(p_ave)))
  fp=np.empty((len(p_ave)))
  ft1=np.empty((len(p_ave)))
  ft2=np.empty((len(p_ave)))

  try:
    absco_dataset_temp = absco_data["absco_"+molecule+"_dataset"]
    absco_p_temp = absco_data["absco_"+molecule+"_p"]
    absco_T_temp = absco_data["absco_"+molecule+"_T"]

    wls = np.linspace(np.where(absco_data["absco_"+molecule+"_wavenumber"] == band_absco_res.min())[0][0],np.where(absco_data["absco_"+molecule+"_wavenumber"] == band_absco_res.max())[0][0],len(band_absco_res))

  except:
    print(molecule," in band ",band_absco_res," not found!")
    return

  for i in range(len(p_ave)):
    ip[i],it1[i],it2[i],fp[i],ft1[i],ft2[i] = interpol_lookup(p_ave[i],absco_p_temp,T_ave[i],absco_T_temp)

  ip_full = np.tile(ip,(len(wls),1)).astype(int)
  it1_full = np.tile(it1,(len(wls),1)).astype(int)
  it2_full = np.tile(it2,(len(wls),1)).astype(int)
  fp_full = np.tile(fp,(len(wls),1))
  ft1_full = np.tile(ft1,(len(wls),1))
  ft2_full = np.tile(ft2,(len(wls),1))
  wls_full = np.tile(wls,(len(ip),1)).T.astype(int)

  #Grab the absorption cross section and convert it from cm^-1 to molecules m^-2
  sigma = (absco_dataset_temp[ip_full,it1_full,wls_full] * (1-fp_full)*(1-ft1_full)\
                + absco_dataset_temp[ip_full,it1_full+1,wls_full] * (1-fp_full)*ft1_full\
                + absco_dataset_temp[ip_full+1,it2_full,wls_full] * fp_full * (1-ft2_full)\
                + absco_dataset_temp[ip_full+1,it2_full+1,wls_full] * fp_full * ft2_full)\
                  * Na / (100.0**2.)

  return sigma




