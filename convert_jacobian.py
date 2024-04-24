import rdata
import numpy as np
from netCDF4 import Dataset
import os, pickle
from collections import defaultdict

class Paths(object):

    def __init__(self):
        super(Paths, self).__init__()
        # define the order in which paths will be searched for input data
        input_folder_order = [
            '/home/jovyan/shared/SSIM-GHG', # on GHG Center's JupyterHub
            os.path.join(os.environ['HOME'], 'Code/Teaching/DA Summer School/2024/SSIM-GHG'), # Sourish's laptop
            ]
        output_folder_order = [
            os.path.join(os.environ['HOME'], 'inversion_output'), # GHG Center's JupyterHub
            os.path.join(os.environ['HOME'], 'Code/Teaching/DA Summer School/2024/output'), # Sourish's laptop
            ]
        for folder_i, folder_o in zip(input_folder_order, output_folder_order):
            if os.path.isdir(folder_i):
                self.data_root = folder_i
                self.output_root = folder_o
                if not os.path.isdir(self.output_root):
                    os.makedirs(self.output_root)
                break
        else:
            raise RuntimeError('No suitable folder found with input data')

        self.jacobi_rda = os.path.join(self.data_root, 'data/trunc_full_jacob_032624_with_dimnames_unit_pulse_4x5_mask_hour_timestamping.rda')
        self.jacobi_nc = os.path.join(self.data_root, 'data/trunc_full_jacob_032624_with_dimnames_unit_pulse_4x5_mask.nc')
        self.obs_rda = os.path.join(self.data_root, 'obs/obs_catalog_041724_unit_pulse_hour_timestamp_witherrors.rda')
        self.obs_nc = os.path.join(self.data_root, 'obs/obs_catalog_041724_unit_pulse_hour_timestamp_witherrors.nc')
        self.obs_by_dataset = os.path.join(self.data_root, 'obs/obs_by_dataset.pickle')
        self.bg_rda = os.path.join(self.data_root, 'data/jacob_bgd_021624.rda')
        self.bg_nc = os.path.join(self.data_root, 'data/jacob_bgd_021624.nc')

class Convert_RDA(Paths):

    def __init__(self, *args, **kwargs):
        super(Convert_RDA, self).__init__()
        self.error_scaling = { # IS data errors are in mol/mol, need to multiply by 10^6 to get ppm
            'IS': 1.0E6,
            }

    def convert_bg(self):
        rda_data = rdata.read_rda(self.bg_rda)
        bg = rda_data['jacob_bgd']
        nobs, ncomp = bg.shape
        with Dataset(self.bg_nc, 'w') as fid:
            fid.createDimension('n_obs', nobs)
            fid.createDimension('n_comp', ncomp)
            v = fid.createVariable('BG', bg.dtype, ('n_obs','n_comp'), zlib=True, shuffle=True, complevel=5)
            v[:] = bg

    def convert_jacobian(self):
        rda_data = rdata.read_rda(self.jacobi_rda)
        J = rda_data['jacob'].values
        nobs, nstate = J.shape
        with Dataset(self.jacobi_nc, 'w') as fid:
            fid.createDimension('n_obs', nobs)
            fid.createDimension('n_state', nstate)
            v = fid.createVariable('Jacobian', J.dtype, ('n_obs','n_state'), zlib=True, shuffle=True, complevel=5, chunksizes=(900,nstate))
            v[:] = J
            setattr(fid, 'original_file', self.jacobi_rda)

    def convert_obs(self):
        # first read the rda file
        rda_data = rdata.read_rda(self.obs_rda)
        obs_data = rda_data['obs_catalog']
        obs_by_dataset = defaultdict(list)
        comp_dict = {'zlib': True, 'complevel': 5, 'shuffle': True}
        with Dataset(self.obs_nc, 'w') as fid:
            fid.createDimension('nobs', len(obs_data))
            # write the coordinates first
            v = fid.createVariable('latitude', np.float32, ('nobs',), **comp_dict)
            v[:] = obs_data['LAT']
            v.units = 'degrees north'
            v = fid.createVariable('longitude', np.float32, ('nobs',), **comp_dict)
            v[:] = obs_data['LON']
            v.units = 'degrees east'
            v = fid.createVariable('time', np.int32, ('nobs',), **comp_dict)
            v[:] = obs_data['TIME']
            v.units = 'Minutes since 2000-01-01 00:00z'

            # convert the dataframe into array
            obs_data = obs_data.iloc[:].values

            # scale the errors to be in same units (hack to fix Andrew's file)
            data_errors = obs_data[:,5]
            for data_type, err_scale in self.error_scaling.items():
                data_idx = obs_data[:,0] == data_type
                data_errors[data_idx] = err_scale * data_errors[data_idx]
            obs_data[:,5] = data_errors

            # store the errors
            v = fid.createVariable('mip_mdm', np.float32, ('nobs',), **comp_dict)
            v[:] = obs_data[:,5]
            v.units = 'micromol/mol'
            v.description = 'Error on obs used in the OCO2 MIP'

            # now store the obs type
            possible_data_types, type_counts = np.unique(obs_data[:,0], return_counts=True)
            for data_type, count in zip(possible_data_types, type_counts):
                dim_name = 'n_%s'%data_type
                fid.createDimension(dim_name, count)

                data_idx = np.where(obs_data[:,0] == data_type)[0]
                var_name = '%s_idx'%data_type
                v = fid.createVariable(var_name, np.int32, (dim_name,), **comp_dict)
                v[:] = data_idx
                v.comment = 'Indices of the original array for the data type'

                var_name = '%s_id'%(data_type)
                if data_type == 'IS':
                    data_dtype = str
                    write_arr = obs_data[:,1][data_idx]
                    v = fid.createVariable(var_name, data_dtype, (dim_name,))
                    # a good place to segregate obs by dataset
                    for obspack_id, i in zip(write_arr, data_idx):
                        dataset = obspack_id.split('~')[1]
                        obs_by_dataset[dataset].append(i)
                else:
                    data_dtype = np.int64
                    write_arr = np.array([int(s) for s in obs_data[:,1][data_idx]], np.int64)
                    v = fid.createVariable(var_name, data_dtype, (dim_name,), **comp_dict)
                v[:] = write_arr

        with open(self.obs_by_dataset, 'wb') as fid:
            pickle.dump(obs_by_dataset, fid, pickle.HIGHEST_PROTOCOL)
