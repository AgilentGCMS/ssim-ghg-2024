import numpy as np
from netCDF4 import Dataset
import os, pickle, time, numbers, calendar
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from convert_jacobian import Paths
from scipy import optimize
from tqdm.auto import tqdm

class Timer(object):
    indent_levels = [] # data to be shared across all instances to keep track of indentation level
    def __init__(self, msg, **kwargs):
        self.msg = msg
        self.print_msg = kwargs['print'] if 'print' in kwargs else True
        self.addendum = {}

    def __enter__(self):
        self.t1 = time.time()
        self.indent_levels.append(1)
        return self.addendum

    def __exit__(self, ex_type, ex_value, ex_traceback):
        t2 = time.time()
        dt = t2-self.t1
        num_indents = len(self.indent_levels)
        if self.print_msg:
            if 'prefix' in self.addendum:
                print_line = "  "*num_indents + "%s %s %s"%(self.addendum['prefix'], self.msg, self.format_dt(dt))
            else:
                print_line = "  "*num_indents + "%s %s"%(self.msg, self.format_dt(dt))
            if 'postfix' in self.addendum:
                print_line = '%s %s'%(print_line, self.addendum['postfix'])
            print(print_line)
        _ = self.indent_levels.pop(0)

    def format_dt(self, dt):
        if dt > 60.0:
            minutes = dt // 60
            seconds = dt % 60.0
            return "%im %.2fs"%(minutes, seconds)
        else:
            return "%.2fs"%dt

class RunSpecs(Paths):

    def __init__(self, *args, **kwargs):
        super(RunSpecs, self).__init__()
        self.start_date = datetime(2014,9,1)
        self.end_date = datetime(2016,9,1) # exclusive
        self.n_region = 22 # there is an additional non-optimized region
        # set the number of months
        n_month = 0
        cur_date = self.start_date
        while cur_date < self.end_date:
            n_month += 1
            cur_date += relativedelta(months=1)
        self.n_month = n_month
        self.land_indices = np.s_[:11*self.n_month]
        self.ocean_indices = np.s_[11*self.n_month:22*self.n_month]

        self.transcom_regions = [ # because someone was not diligent about putting this in the data files
            'North American Boreal', 'North American Temperate', 'South American Tropical', 'South American Temperate',
            'Northern Africa', 'Southern Africa', 'Eurasian Boreal', 'Eurasian Temperate', 'Tropical Asia', 'Australia',
            'Europe', 'North Pacific Temperate', 'West Pacific Tropical', 'East Pacific Tropical', 'South Pacific Temperate',
            'Northern Ocean', 'North Atlantic Temperate', 'Atlantic Tropical', 'South Atlantic Temperate', 'Southern Ocean',
            'Indian Tropical', 'South Indian Temperate',
            ]

class Fluxes(RunSpecs):

    def __init__(self, *args, **kwargs):
        super(Fluxes, self).__init__()
        self.verbose = kwargs['verbose'] if 'verbose' in kwargs else True
        self.region_mask = os.path.join(self.data_root, 'transcom/TRANSCOM_mask_GEOS_Chem_4x5.nc')
        self.original_mask = os.path.join(self.data_root, 'transcom/TRANSCOM_mask_original_1x1.nc')
        self.Re = 6.371E+6
        self.dS = {'regular': {}, 'geos': {}}
        self.calculate_region_areas()

    def calculate_region_areas(self):
        with Dataset(self.original_mask, 'r') as fid:
            mask64 = fid.variables['mask64'][:]
        nlat, nlon = mask64.shape
        dS = self.surface_area(nlat, nlon, 'regular')
        self.transcom_region_areas = np.zeros(len(self.transcom_regions), np.float64)
        for i in range(1,len(self.transcom_regions)+1):
            region_mask = np.logical_and(mask64 > i-0.1, mask64 < i+0.1) # avoid floating point equality test
            self.transcom_region_areas[i-1] = dS[region_mask].sum() # m^2

    def surface_area(self, nlat, nlon, variant):
        if (nlat,nlon) not in self.dS[variant]:
            geom = self.geometry(nlat, nlon, variant)
            self.dS[variant][(nlat,nlon)] = geom['surface_area']
        return self.dS[variant][(nlat,nlon)]

    def geometry(self, nlat, nlon, variant='geos'):
        if variant == 'geos':
            dlat = 180.0/(nlat-1) # bulk of the sphere
            dlon = 360.0/nlon # bulk of the sphere
            lat_edges = np.zeros(nlat+1, dtype=np.float64)
            lat_edges[0] = -90.0
            lat_edges[-1] = 90.0
            lat_edges[1:nlat] = np.linspace(-90.0+dlat/2., 90.0-dlat/2., nlat-1)
            lon_edges = np.linspace(-180.0-dlon/2., 180.0-dlon/2., nlon+1)
        elif variant == 'regular':
            dlat = 180.0/nlat
            dlon = 360.0/nlon
            lat_edges = np.linspace(-90., 90., nlat+1)
            lon_edges = np.linspace(-180.,180.,nlon+1)
        else:
            raise RuntimeError("Variant is %s, not allowed"%variant)

        # now the surface area
        dlon = (np.pi/180.) * dlon
        dS = np.zeros((nlat+1, nlon), np.float64)
        lat_edges_rad = (np.pi/180.0) * lat_edges
        for i, lat in enumerate(lat_edges_rad):
            dS[i] = self.Re * self.Re * dlon * np.sin(lat)
        dS = np.diff(dS, axis=0)

        return {
            'lat_edges': lat_edges,
            'lon_edges': lon_edges,
            'surface_area': dS,
            }

    def read_ct2022_flux(self, year, month, var_name):
        file_name = os.path.join(self.data_root, 'fluxes/CT2022/CT2022.flux1x1.%04i%02i.nc'%(year, month))
        with Dataset(file_name, 'r') as fid:
            flux_1x1 = fid.variables[var_name][0] # mol/m^2/s
        return 1.0E6 * flux_1x1 # micromol/m^2/s

    def read_sib4_flux(self, year, month, var_name):
        file_name = os.path.join(self.data_root, 'fluxes/SiB4/isotm5_%04i-%02i.nc'%(year, month))
        with Dataset(file_name, 'r') as fid:
            flux_1x1 = fid.variables[var_name][0] # mol/m^2/s
        return 1.0E6 * flux_1x1 # micromol/m^2/s

    def convert_1x1_flux_to_4x5(self, flux_1x1):
        # The 4x5 GEOS Chem grid has
        #   lat edges -90, -88, -84, -80, ... 80, 84, 88, 90
        #   lon_edges -182.5, -177.5, -172.5, ... 172.5, 177.5
        # So I have to first put the input 1x1 flux onto a 1x0.5 grid
        flux_1x05 = np.repeat(flux_1x1, 2, axis=1)
        # flux_1x05 has lon edges at -180, -179.5, -179, etc., so we need to rotate it by 5 cells in the longitude direction
        flux_1x05 = np.roll(flux_1x05, 5, axis=1) # now left edge is at -182.5
        # convert to per cell
        flux_1x05 = flux_1x05 * self.surface_area(180, 720, 'regular')
        # now fill a 4x5 array using reduceat
        flux_4x5 = np.zeros((46,72), dtype=np.float64)
        flux_4x5[0,:] = np.add.reduceat(flux_1x05[0:2], np.arange(0, flux_1x05.shape[1], 10), axis=1).sum(axis=0)
        flux_4x5[45,:] = np.add.reduceat(flux_1x05[178:180], np.arange(0, flux_1x05.shape[1], 10), axis=1).sum(axis=0)
        flux_4x5[1:45,:] = np.add.reduceat(
            np.add.reduceat(flux_1x05[2:178,:], np.arange(0,flux_1x05.shape[0]-4,4), axis=0),
            np.arange(0, flux_1x05.shape[1], 10), axis=1)
        # convert back to per area
        if (46,72) not in self.dS['geos']:
            geom = self.geometry(46,72,'geos')
            self.dS['geos'][(46,72)] = geom['surface_area']
        flux_4x5 = flux_4x5 / self.surface_area(46, 72, 'geos')

        return flux_4x5

    def convert_4x5_flux_to_statevector(self, flux_4x5, year, month, **kwargs):
        state_vector = np.zeros(self.n_region*self.n_month, np.float64)
        # if year/month is outside our Jacobian simulation window, return zeros
        jd = datetime(year, month, 1)
        if jd < self.start_date or jd >= self.end_date:
            return state_vector

        dS = self.surface_area(46, 72, 'geos')
        region_masks = {}
        with Dataset(self.region_mask, 'r') as fid:
            for ireg in range(self.n_region):
                region_masks[ireg] = fid.variables['Region%02i'%(1+ireg)][:] # Region00 is unoptimized and not part of the Jacobian

        # get time index within one region
        for i_month in range(self.n_month):
            if self.start_date + relativedelta(months=i_month) == jd:
                break
        else:
            raise RuntimeError('Cannot find time index for %04i-%02i'%(year, month))

        if 'regions' in kwargs:
            region_range = kwargs['regions']
        else:
            region_range = range(self.n_region)
        for ireg in region_range:
            mask = region_masks[ireg]
            mean_flux = (flux_4x5 * dS * mask).sum() / (dS * mask).sum()
            state_vector[ireg*self.n_month+i_month] = mean_flux # in units of micromoles/m^2/s

        # state_vector is now in micromoles/m^2/s, need to convert to Kg CO2/m^2/s
        state_vector = 0.001 * 44.01 * 1.0E-6 * state_vector

        return state_vector

    def construct_state_vector_from_sib4(self, smush_regions=True, land_var='I2b', ocean_var='Fnetoce'):
        ym_tuples = []
        cur_month = self.start_date
        while cur_month < self.end_date:
            ym_tuples.append((cur_month.year, cur_month.month))
            cur_month += relativedelta(months=1)

        state_vec = 0.0
        for year, month in tqdm(ym_tuples, desc='Converting SiB4 to state vector'):
            flux_nee = self.read_sib4_flux(year, month, land_var)
            flux_oce = self.read_sib4_flux(year, month, ocean_var)
            if smush_regions:
                flux_4x5 = self.convert_1x1_flux_to_4x5(flux_nee+flux_oce)
                state_vec = state_vec + self.convert_4x5_flux_to_statevector(flux_4x5, year, month)
            else:
                flux_4x5 = self.convert_1x1_flux_to_4x5(flux_nee)
                state_vec = state_vec + self.convert_4x5_flux_to_statevector(flux_4x5, year, month, regions=range(11))
                flux_4x5 = self.convert_1x1_flux_to_4x5(flux_oce)
                state_vec = state_vec + self.convert_4x5_flux_to_statevector(flux_4x5, year, month, regions=range(11,22))

        return state_vec

    def construct_state_vector_from_ct2022(self, smush_regions=True):
        ym_tuples = []
        cur_month = self.start_date
        while cur_month < self.end_date:
            ym_tuples.append((cur_month.year, cur_month.month))
            cur_month += relativedelta(months=1)

        state_vec = 0.0
        for year, month in tqdm(ym_tuples, desc='Converting CT2022 to state vector'):
            flux_nee = self.read_ct2022_flux(year, month, 'bio_flux_opt')
            flux_oce = self.read_ct2022_flux(year, month, 'ocn_flux_opt')
            if smush_regions:
                flux_4x5 = self.convert_1x1_flux_to_4x5(flux_nee+flux_oce)
                state_vec = state_vec + self.convert_4x5_flux_to_statevector(flux_4x5, year, month)
            else:
                flux_4x5 = self.convert_1x1_flux_to_4x5(flux_nee)
                state_vec = state_vec + self.convert_4x5_flux_to_statevector(flux_4x5, year, month, regions=range(11))
                flux_4x5 = self.convert_1x1_flux_to_4x5(flux_oce)
                state_vec = state_vec + self.convert_4x5_flux_to_statevector(flux_4x5, year, month, regions=range(11,22))

        return state_vec

class Transport(RunSpecs):

    def __init__(self, *args, **kwargs):
        super(Transport, self).__init__()
        self.verbose = kwargs['verbose'] if 'verbose' in kwargs else True
        self.jacobian_file = {
            'GC': os.path.join(self.data_root, 'jacobians/trunc_full_jacob_032624_with_dimnames_unit_pulse_4x5_mask.nc'),
            'TM5': None,
            }
        self.background = {
            # 'GC': os.path.join(self.data_root, 'jacobians/jacob_bgd_021624.nc'), # probably not correct
            'GC': self.bg_nc,
            'TM5': None,
            }
        self.Jacobian = {}
        self.background_ppm = 400.0 # arbitrary choice

    def transport(self, state_vector, model='GC', add_bg=True):
        if model not in self.Jacobian:
            with Timer("Read Jacobian in ", print=self.verbose):
                with Dataset(self.jacobian_file[model], 'r') as fid:
                    self.Jacobian[model] = fid.variables['Jacobian'][:]

        if add_bg:
            with Timer("Added background in ", print=self.verbose):
                with Dataset(self.background[model], 'r') as fid:
                    bg = fid.variables['BG'][:]
                bg = bg[:,1:3].sum(axis=1) + self.background_ppm
        else:
            bg = 0.0

        with Timer("Simulated transport in ", print=self.verbose):
            obs = np.matmul(self.Jacobian[model], state_vector) + bg # ppm

        return obs

class Observations(RunSpecs):

    def __init__(self, *args, **kwargs):
        super(Observations, self).__init__()
        self.verbose = kwargs['verbose'] if 'verbose' in kwargs else True
        self.noaa_observatories = ['mlo', 'spo', 'smo', 'brw']
        self.mbl_sites = [
            'abp', 'alt', 'ams', 'asc', 'avi', 'azr', 'bme', 'bmw', 'brw', 'cba', 'cgo', 'chr', 'crz', 'gmi', 'hba',
            'ice', 'key', 'kum', 'mbc', 'mhd', 'mid', 'poc', 'psa', 'rpb', 'shm', 'smo', 'spo', 'stm', 'syo', 'zep',
            ]
        self.noaa_tall_tower_datasets = [
            'co2_amt_tower-insitu_1_allvalid-107magl', 'co2_amt_surface-pfp_1_allvalid-107magl',
            'co2_bao_tower-insitu_1_allvalid-300magl', 'co2_bao_surface-pfp_1_allvalid-300magl',
            'co2_lef_tower-insitu_1_allvalid-396magl',
            'co2_sct_tower-insitu_1_allvalid-305magl', 'co2_sct_surface-pfp_1_allvalid-305magl',
            'co2_wbi_tower-insitu_1_allvalid-379magl', 'co2_wbi_surface-pfp_1_allvalid-379magl',
            'co2_wgc_tower-insitu_1_allvalid-483magl',
            'co2_wkt_tower-insitu_1_allvalid-457magl']
        self.site_code_to_dataset = {
            'lef': ['co2_lef_tower-insitu_1_allvalid-396magl'],
            'amt': ['co2_amt_tower-insitu_1_allvalid-107magl'],
            'wkt': ['co2_wkt_tower-insitu_1_allvalid-457magl'],
            'wbi': ['co2_wbi_tower-insitu_1_allvalid-379magl', 'co2_wbi_surface-pfp_1_allvalid-379magl'],
            'sct': ['co2_sct_surface-pfp_1_allvalid-305magl', 'co2_sct_tower-insitu_1_allvalid-305magl'],
            'wgc': ['co2_wgc_tower-insitu_1_allvalid-483magl'],
            }
        self.unassim_mdm = 1.0E36 # ppm

    def convert_datetime_to_decimal_year(self, datetime_array):
        ret_arr = np.zeros(len(datetime_array), dtype=np.float64)
        for i,d in enumerate(datetime_array):
            seconds_elapsed = (d-datetime(d.year,1,1)).total_seconds()
            seconds_total = 86400 * (365 + calendar.isleap(d.year))
            ret_arr[i] = d.year + seconds_elapsed/seconds_total
        return ret_arr

    def get_datasets_from_site(self, site_code):
        with open(self.obs_by_dataset, 'rb') as fid:
            obs_dict = pickle.load(fid)
        ret_list = [k for k in obs_dict.keys() if k.startswith('co2_%s'%site_code.lower())]
        return ret_list

    def get_indices_from_datasets(self, datasets):
        with open(self.obs_by_dataset, 'rb') as fid:
            obs_dict = pickle.load(fid)
        ret_arr = []
        for k in datasets:
            ret_arr.extend(obs_dict[k])
        return np.array(ret_arr)

    def get_time_series(self, obs_vector, dataset_list):
        # first get the indices
        obs_indices = self.get_indices_from_datasets(dataset_list)
        # now read the obs times
        with Dataset(self.obs_nc, 'r') as fid:
            obs_times = fid.variables['time'][:]
        # subselect and sort by time
        obs = obs_vector[obs_indices]
        obs_times = obs_times[obs_indices]
        sort_order = np.argsort(obs_times)
        obs = obs[sort_order]
        obs_times = obs_times[sort_order]
        # now convert to datetime object
        obs_times = np.array([datetime(2000,1,1) + timedelta(minutes=int(m)) for m in obs_times])
        return obs_times, obs

    def save_timeseries_by_site(self, obs_vector, model_vector, mdm_vector, output_folder, site_codes, stage, append):
        # stage is 'apri' or 'apos'
        # append (to existing file) is True or False
        if not os.path.isdir(output_folder):
            os.makedirs(output_folder)

        # First get all times
        with Dataset(self.obs_nc, 'r') as fid:
            obs_times = fid.variables['time'][:]

        # Now loop over site
        for site in site_codes:
            if site in self.site_code_to_dataset:
                ds_ = self.site_code_to_dataset[site]
            else:
                ds_ = self.get_datasets_from_site(site)
                ds_ = [s for s in ds_ if 'flask' in s]
            obs_idx = self.get_indices_from_datasets(ds_)

            if len(obs_idx) == 0:
                continue

            site_times = obs_times[obs_idx]
            site_obs = obs_vector[obs_idx]
            site_model = model_vector[obs_idx]
            site_mdm = mdm_vector[obs_idx]
            # sort by time
            sort_order = np.argsort(site_times)
            site_times = [datetime(2000,1,1) + timedelta(minutes=int(m)) for m in site_times[sort_order]]
            site_obs = site_obs[sort_order]
            site_model = site_model[sort_order]
            site_mdm = site_mdm[sort_order]

            # now write to file
            comp_dict = {'zlib': True, 'complevel': 6, 'shuffle': True}
            fname = os.path.join(output_folder, 'timeseries_%s.nc'%site)
            file_mode = 'a' if append else 'w'

            with Dataset(fname, file_mode) as fid:
                if not append:
                    fid.createDimension('obs', len(obs_idx))
                    fid.createDimension('idate', 6)

                if append:
                    v = fid.variables['integer_times']
                else:
                    v = fid.createVariable('integer_times', np.int16, ('obs','idate'), **comp_dict)
                v[:] = np.array([d.timetuple()[:6] for d in site_times], dtype=np.int16)

                if append:
                    v = fid.variables['decimal_times']
                else:
                    v = fid.createVariable('decimal_times', np.float64, ('obs',), **comp_dict)
                v[:] = self.convert_datetime_to_decimal_year(site_times)

                if append:
                    v = fid.variables['observed']
                else:
                    v = fid.createVariable('observed', np.float64, ('obs',), **comp_dict)
                    v.units = 'ppm'
                v[:] = site_obs

                if append:
                    v = fid.variables['obs_error']
                else:
                    v = fid.createVariable('obs_error', np.float64, ('obs',), **comp_dict)
                    v.units = 'ppm'
                v[:] = site_mdm

                model_var_name = 'modeled_%s'%stage
                if append and model_var_name in fid.variables:
                    v = fid.variables[model_var_name]
                else:
                    v = fid.createVariable(model_var_name, np.float64, ('obs',), **comp_dict)
                    v.units = 'ppm'
                v[:] = site_model

                if append:
                    setattr(fid, 'modification_date', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
                else:
                    setattr(fid, 'creation_date', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

class Construct_Synthetic_Obs(RunSpecs):

    def __init__(self, true_flux='CT2022', transport='GC', *args, **kwargs):
        # Let's take CT2022 to be "truth", and SiB4 to be prior by default, unless otherwise specified
        super(Construct_Synthetic_Obs, self).__init__()
        self.verbose = kwargs['verbose'] if 'verbose' in kwargs else True
        self.true_flux = true_flux
        self.transport = transport
        self.synth_obs_file = os.path.join(self.data_root, 'obs/synthetic_obs_%s_%s.nc'%(true_flux, transport))

    def make_obs(self):
        fl = Fluxes()
        if self.true_flux.lower() == 'ct2022':
            state_vec = fl.construct_state_vector_from_ct2022(smush_regions=False)
        elif self.true_flux.lower() == 'sib4':
            state_vec = fl.construct_state_vector_from_sib4(smush_regions=False)
        else:
            raise RuntimeError('Unknown flux type %s'%self.true_flux)

        t = Transport()
        obs = t.transport(state_vec, model=self.transport, add_bg=True)

        with Dataset(self.synth_obs_file, 'w') as fid:
            fid.createDimension('observations', len(obs))
            v = fid.createVariable('co2', np.float64, ('observations',), zlib=True, complevel=6, shuffle=True)
            v[:] = obs
            setattr(v, 'units', 'micromol/mol')

class Var4D_Components(RunSpecs):

    def __init__(self, project, *args, **kwargs):
        super(Var4D_Components, self).__init__()
        self.project = project
        self.verbose = kwargs['verbose'] if 'verbose' in kwargs else True
        self.store_intermediate_states = kwargs['store_intermediate'] if 'store_intermediate' in kwargs else False
        self.output_dir = os.path.join(self.output_root, project)
        # create a bunch of instances of the classes required
        self.trans_op = Transport(verbose=self.verbose)
        self.flux_cons = Fluxes(verbose=self.verbose)
        self.obs_cons = Observations(verbose=self.verbose)

    def setup_obs_errors(self, **kwargs):
        # kwargs is a dictionary such as {'oco2': True, 'tccon': True, 'is': True, 'mbl': True, 'sites': ['mlo','spo','brw']}
        # by default nothing is assimilated
        # first get the total number of obs
        with Dataset(self.obs_nc, 'r') as fid:
            nobs = len(fid.dimensions['nobs'])
            oco2_idx = fid.variables['OCO2_idx'][:]
            is_idx = fid.variables['IS_idx'][:]
            tccon_idx = fid.variables['TCCON_idx'][:]
            mip_mdm = fid.variables['mip_mdm'][:]

        # by default assimilate nothing
        obs_err = self.obs_cons.unassim_mdm * np.ones(nobs, dtype=np.float64)

        if 'oco2' in kwargs and kwargs['oco2']:
            obs_err[oco2_idx] = mip_mdm[oco2_idx]

        if 'tccon' in kwargs and kwargs['tccon']:
            obs_err[tccon_idx] = mip_mdm[tccon_idx]

        if 'is' in kwargs and kwargs['is']:
            obs_err[is_idx] = mip_mdm[is_idx]

        if 'mbl' in kwargs and kwargs['mbl']:
            # make list of MBL datasets
            mbl_datasets = []
            for site in self.obs_cons.mbl_sites:
                ds_ = self.obs_cons.get_datasets_from_site(site)
                ds_ = [s for s in ds_ if 'flask' in s]
                mbl_datasets.extend(ds_)
            mbl_idx = self.obs_cons.get_indices_from_datasets(mbl_datasets)
            obs_err[mbl_idx] = mip_mdm[mbl_idx]

        if 'noaa_towers' in kwargs and kwargs['noaa_towers']:
            tower_idx = self.obs_cons.get_indices_from_datasets(self.obs_cons.noaa_tall_tower_datasets)
            obs_err[tower_idx] = mip_mdm[tower_idx]

        if 'noaa_observatories' in kwargs and kwargs['noaa_observatories']:
            obs_datasets = []
            for site in self.obs_cons.noaa_observatories:
                ds_ = self.obs_cons.get_datasets_from_site(site)
                ds_ = [s for s in ds_ if 'flask' in s]
                obs_datasets.extend(ds_)
            obs_idx = self.obs_cons.get_indices_from_datasets(obs_datasets)
            obs_err[obs_idx] = mip_mdm[obs_idx]

        # can specify individual site codes (flask only)
        if 'sites' in kwargs:
            assim_datasets = []
            for site in kwargs['sites']:
                ds_ = self.obs_cons.get_datasets_from_site(site)
                ds_ = [s for s in ds_ if 'flask' in s]
                assim_datasets.extend(ds_)
            assim_idx = self.obs_cons.get_indices_from_datasets(assim_datasets)
            obs_err[assim_idx] = mip_mdm[assim_idx]

        # finally, can specify individual datasets to assimilate
        if 'datasets' in kwargs:
            assim_idx = self.obs_cons.get_indices_from_datasets(kwargs['datasets'])
            obs_err[assim_idx] = mip_mdm[assim_idx]

        # we have the option of specifying a single MDM for all assimilated obs
        if 'uniform_mdm' in kwargs:
            assim_idx = obs_err < 0.5*self.obs_cons.unassim_mdm
            obs_err[assim_idx] = kwargs['uniform_mdm']

        # print the number of obs assimilated
        print('%i of %i obs will be assimilated'%( (obs_err < 0.5*self.obs_cons.unassim_mdm).sum(), len(obs_err) ))

        return obs_err

    def add_obs_bias(self, bias_ppm, **kwargs):
        with Dataset(self.obs_nc, 'r') as fid:
            nobs = len(fid.dimensions['nobs'])
            # change oco2_idx and is_idx to be vectors of boolean (full length)
            oco2_idx = np.zeros(nobs, bool)
            oco2_idx[fid.variables['OCO2_idx'][:]] = True
            is_idx = np.zeros(nobs, bool)
            is_idx[fid.variables['IS_idx'][:]] = True
            obs_lat = fid.variables['latitude'][:]
            obs_lon = fid.variables['longitude'][:]
            time_origin = datetime.strptime(fid.variables['time'].units, 'Minutes since %Y-%m-%d %H:%Mz')
            obs_times = [time_origin + timedelta(minutes=int(m)) for m in fid.variables['time'][:]]

        # We define a few types of biases, students are encouraged to try to add additional bias types
        obs_bias = np.zeros_like(self.obs_vec)
        bias_idx = np.ones(nobs, bool) # all obs are biased by default

        if 'lat_min' in kwargs:
            bias_idx = np.logical_and(bias_idx, obs_lat >= kwargs['lat_min'])
        if 'lat_max' in kwargs:
            bias_idx = np.logical_and(bias_idx, obs_lat <= kwargs['lat_max'])
        if 'lon_min' in kwargs:
            bias_idx = np.logical_and(bias_idx, obs_lon >= kwargs['lon_min'])
        if 'lon_max' in kwargs:
            bias_idx = np.logical_and(bias_idx, obs_lon <= kwargs['lon_max'])
        if 'platform' in kwargs:
            if kwargs['platform'].lower() == 'is':
                bias_idx = np.logical_and(bias_idx, is_idx)
            if kwargs['platform'].lower() == 'oco2':
                bias_idx = np.logical_and(bias_idx, oco2_idx)
        if 'months' in kwargs:
            month_idx = np.array([d.month in kwargs['months'] for d in obs_times])
            bias_idx = np.logical_and(bias_idx, month_idx)

        obs_bias[bias_idx] = obs_bias[bias_idx] + bias_ppm
        print('Added %.1f ppm bias to %i out of %i obs'%(bias_ppm, bias_idx.sum(), len(bias_idx)))

        self.obs_vec = self.obs_vec + obs_bias

    def setup_obs(self, true_flux='CT2022', trans_model='GC', obs_to_assim={}):
        if true_flux.lower() == 'ct2022':
            state_vec = self.flux_cons.construct_state_vector_from_ct2022(smush_regions=False)
        elif true_flux.lower() == 'sib4':
            state_vec = self.flux_cons.construct_state_vector_from_sib4(smush_regions=False)
        else:
            raise RuntimeError('Unknown true flux %s specified'%true_flux)

        self.obs_vec = self.trans_op.transport(state_vec, model=trans_model, add_bg=True)
        self.obs_err = self.setup_obs_errors(**obs_to_assim)
        self.true_flux = state_vec

    def setup_corr(self, **kwargs):
        # the state vector is length 528, arranged as 'region 1 month 1', 'region 1 month 2', ... 'region 22 month 23', 'region 22 month 24'
        temp_corr_mat = np.zeros((self.n_month, self.n_month), np.float64)
        # temp_corr being a number N means exponentially decaying with an e-folding length of N months
        try:
            temp_corr = kwargs['temp_corr']
        except KeyError:
            raise RuntimeError('temp_corr must be specified')

        if isinstance(temp_corr, numbers.Number):
            if temp_corr == 0.0:
                temp_corr_mat = np.eye(self.n_month)
            else:
                for i in range(self.n_month):
                    for j in range(self.n_month):
                        temp_corr_mat[i,j] = np.exp(-np.abs(i-j)/temp_corr)
        elif temp_corr.lower() == 'clim':
            # each month is highly correlated with the same month in other years
            temp_corr_mat = np.eye(self.n_month) # the diagonal
            for i in range(self.n_month):
                for j in range(i+1, self.n_month):
                    if (i-j) % 12 == 0:
                        temp_corr_mat[i,j] = 0.9
                        temp_corr_mat[j,i] = 0.9
        else:
            raise RuntimeError('Please implement how to handle temp_corr choice %s'%temp_corr)

        hor_corr_mat = np.eye(self.n_region)
        if 'region_pairs' in kwargs:
            # syntax is {(region 1, region 2): corr, (region 3, region 4): corr}, etc.
            # regions can either be region names or numbers
            all_regions = [s.lower() for s in self.transcom_regions]
            for (reg1,reg2), r in kwargs['region_pairs'].items():
                if isinstance(reg1, int):
                    ireg1 = reg1
                else:
                    ireg1 = all_regions.index(reg1.lower())
                if isinstance(reg2, int):
                    ireg2 = reg2
                else:
                    ireg2 = all_regions.index(reg2.lower())

                hor_corr_mat[ireg1,ireg2] = r
                hor_corr_mat[ireg2,ireg1] = r

        # now multiply the two to get the full correlation matrix, remembering that in the state vector,
        # the lateral dimension (regions) varies slower than the time dimension (months)
        corr_matrix = np.kron(hor_corr_mat, temp_corr_mat)

        return corr_matrix

    def setup_cov(self, corr_matrix, unc_vector):
        # corr_matrix is n_state x n_state, unc_vector is n_state
        n_state = unc_vector.shape[0]
        assert (corr_matrix.shape[0] == n_state and corr_matrix.shape[1] == n_state), 'Shapes of correlation matrix and uncertainty vector are non-compliant'
        cov_matrix = np.zeros((n_state, n_state), np.float64)
        for i in range(n_state):
            for j in range(n_state):
                cov_matrix[i,j] = corr_matrix[i,j] * unc_vector[i] * unc_vector[j]

        return cov_matrix

    def var4d_setup(self, **kwargs):
        obs_to_assim = kwargs['obs_to_assim']  if 'obs_to_assim' in kwargs else {'sites': ['mlo','spo']}
        prior_unc_scale = kwargs['prior_unc_scale'] if 'prior_unc_scale' in kwargs else {'land': 0.25, 'ocean': 0.25}
        prior_unc_source = kwargs['prior_unc_source'] if 'prior_unc_source' in kwargs else 'nee'
        corr_structure = kwargs['corr_structure'] if 'corr_structure' in kwargs else {'temp_corr': 2.0}
        perturb_obs = kwargs['perturb_obs'] if 'perturb_obs' in kwargs else False
        perturb_flux = kwargs['perturb_flux'] if 'perturb_flux' in kwargs else False

        with Timer("Created true obs in ", print=self.verbose):
            # set up the obs with CT2022 as truth
            self.setup_obs(true_flux='CT2022', obs_to_assim=obs_to_assim)
            if perturb_obs:
                assim_idx = self.obs_err < 0.5*self.obs_cons.unassim_mdm
                self.obs_vec[assim_idx] = self.obs_vec[assim_idx] + np.random.standard_normal(assim_idx.sum()) * self.obs_err[assim_idx]
                print('Perturbed %i of %i obs'%(assim_idx.sum(), len(assim_idx)))

        with Timer("Prior fluxes and covariance setup in ", print=self.verbose):
            # set up the prior, which is SiB4
            self.state_prior = self.flux_cons.construct_state_vector_from_sib4()
            # set up prior uncertainty
            self.unc_prior = np.zeros_like(self.state_prior)
            if prior_unc_source.lower() == 'nee':
                self.unc_prior[self.land_indices] = prior_unc_scale['land'] * np.abs(self.state_prior[self.land_indices])
            elif prior_unc_source.lower() == 'gpp':
                state_unc = self.flux_cons.construct_state_vector_from_sib4(smush_regions=False, land_var='Fgpp_SiB4')
                self.unc_prior[self.land_indices] = prior_unc_scale['land'] * np.abs(state_unc[self.land_indices])
            elif prior_unc_source.lower() == 'reco':
                state_unc = self.flux_cons.construct_state_vector_from_sib4(smush_regions=False, land_var='I4')
                self.unc_prior[self.land_indices] = prior_unc_scale['land'] * np.abs(state_unc[self.land_indices])
            # not many choices to set ocean errors
            self.unc_prior[self.ocean_indices] = prior_unc_scale['ocean'] * np.abs(self.state_prior[self.ocean_indices])
            # set up prior covariance matrix, assuming a 2 month decay and no cross-region correlation
            prior_corr = self.setup_corr(**corr_structure)
            # make correlation expicitly symmetric to avoid roundoff errors
            self.prior_corr = 0.5*(prior_corr + prior_corr.T)
            prior_cov = self.setup_cov(self.prior_corr, self.unc_prior)
            # make covariance expicitly symmetric to avoid roundoff errors
            self.prior_cov = 0.5*(prior_cov + prior_cov.T)
            # calculate L such that LL^T = self.prior_cov
            self.L = np.linalg.cholesky(self.prior_cov)
            # set the prior in preconditioned space to be a vector of zeros
            n_state = self.state_prior.shape[0]
            self.state_prior_preco = np.zeros(n_state, dtype=np.float64)
            if perturb_flux:
                random_vector = np.random.standard_normal(n_state)
                print('Perturbed preconditioned state with mean %.4f and standard deviation %.4f'%(random_vector.mean(), random_vector.std()))
                self.state_prior = self.state_to_flux(random_vector)

        self.progress_bars = {}
        if self.store_intermediate_states:
            self.interim_states = {'x': [], 'fun': []}

    def var4d_done(self):
        output_file = os.path.join(self.output_dir, 'optim_summary.nc')
        comp_dict = {'zlib': True, 'complevel': 6, 'shuffle': True}
        res = self.optim_result

        if self.store_intermediate_states:
            for k,v in self.interim_states.items():
                self.interim_states[k] = np.array(v)

        with Dataset(output_file, 'w') as fid:
            fid.createDimension('cost_eval', len(self.optim_diags['cost_function']))
            fid.createDimension('grad_eval', len(self.optim_diags['gradient_norm']))
            fid.createDimension('n_state', res.x.shape[0])
            fid.createDimension('n_region', self.n_region)
            fid.createDimension('n_month', self.n_month)
            fid.createDimension('n_ymd', 3)

            v = fid.createVariable('cost_function_evals', np.float64, ('cost_eval',), **comp_dict)
            v[:] = np.array(self.optim_diags['cost_function'])
            setattr(v, 'comment', 'All cost function evaluations during optimization')

            v = fid.createVariable('gradient_norm_evals', np.float64, ('grad_eval',), **comp_dict)
            v[:] = np.array(self.optim_diags['gradient_norm'])
            setattr(v, 'comment', 'L2 norm of all gradient evaluations during optimization')

            setattr(fid, 'iterations', np.int32(res.nit))
            setattr(fid, 'converged', np.int32(res.success))
            setattr(fid, 'final_cost_function', res.fun)
            setattr(fid, 'nfev', np.int32(res.nfev))
            setattr(fid, 'njev', np.int32(res.njev))
            setattr(fid, 'nit', np.int32(res.nit))
            if 'nhev' in dir(res):
                setattr(fid, 'nhev', np.int32(res.nhev))

            v = fid.createVariable('x_final', np.float64, ('n_state',), **comp_dict)
            v[:] = res.x
            setattr(v, 'comment', 'Final state vector in preconditioned space')

            v = fid.createVariable('j_final', np.float64, ('n_state',), **comp_dict)
            v[:] = res.jac
            setattr(v, 'comment', 'Final Jacobian in preconditioned space')

            v = fid.createVariable('prior_cov', np.float64, ('n_state','n_state'), **comp_dict)
            v[:] = self.prior_cov
            setattr(v, 'comment', 'Prior covariance')
            setattr(v, 'units', '(Kg CO2/m^2/s)^2')

            if 'hess_inv' in dir(res):
                v = fid.createVariable('hess_inv', np.float64, ('n_state','n_state'), **comp_dict)
                if isinstance(res.hess_inv, optimize.LbfgsInvHessProduct):
                    v[:] = res.hess_inv.todense()
                else:
                    v[:] = res.hess_inv
                setattr(v, 'comment', 'Approximate inverse Hessian in preconditioned space')

                v = fid.createVariable('poste_cov', np.float64, ('n_state','n_state'), **comp_dict)
                inv_hess_preco = fid.variables['hess_inv'][:]
                v[:] = np.matmul(self.L, np.matmul(inv_hess_preco, self.L.T))
                setattr(v, 'comment', 'Approximate posterior covariance calculated from inverse hessian')
                setattr(v, 'units', '(Kg CO2/m^2/s)^2')

            if self.store_intermediate_states:
                fid.createDimension('iterations', res.nit)
                v = fid.createVariable('cost_function', np.float64, ('iterations',), **comp_dict)
                v[:] = self.interim_states['fun']
                setattr(v, 'comment', 'Cost function for each iteration of optimizer')

                v = fid.createVariable('states', np.float64, ('iterations', 'n_state'), **comp_dict)
                v[:] = self.interim_states['x']
                setattr(v, 'comment', 'Value of state vector in preconditioned space for each iteration of optimizer')

            v = fid.createVariable('prior_flux', np.float64, ('n_region','n_month'), **comp_dict)
            for i in range(self.n_region):
                v[i] = self.state_prior[i*self.n_month:(i+1)*self.n_month]
            setattr(v, 'comment', 'Prior flux')
            setattr(v, 'units', 'Kg CO2/m^2/s')

            v = fid.createVariable('poste_flux', np.float64, ('n_region','n_month'), **comp_dict)
            for i in range(self.n_region):
                v[i] = self.state_poste[i*self.n_month:(i+1)*self.n_month]
            setattr(v, 'comment', 'Posterior flux')
            setattr(v, 'units', 'Kg CO2/m^2/s')

            v = fid.createVariable('true_flux', np.float64, ('n_region','n_month'), **comp_dict)
            for i in range(self.n_region):
                v[i] = self.true_flux[i*self.n_month:(i+1)*self.n_month]
            setattr(v, 'comment', 'The flux that was used to create the obs')
            setattr(v, 'units', 'Kg CO2/m^2/s')

            assert self.n_region == len(self.flux_cons.transcom_regions), "Length of n_region does not match number of TRANSCOM regions"

            v = fid.createVariable('region_areas', np.float64, ('n_region',), **comp_dict)
            v[:] = self.flux_cons.transcom_region_areas
            setattr(v, 'units', 'm^2')
            setattr(v, 'region_names', ','.join(self.flux_cons.transcom_regions))

            ymd_tuples = []
            cur_date = self.start_date
            while cur_date < self.end_date:
                ymd_tuples.append(cur_date.timetuple()[:3])
                cur_date += relativedelta(months=1)
            v = fid.createVariable('time_coordinate', np.int16, ('n_month','n_ymd'), **comp_dict)
            v[:] = np.array(ymd_tuples, dtype=np.int16)
            setattr(v, 'description', 'year/month/day denoting the beginning of each month corresponding to monthly fluxes')

            setattr(fid, 'creation_date', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    def var4d_chain(self, **kwargs):
        max_iter = kwargs['max_iter'] if 'max_iter' in kwargs else 500
        max_gradnorm = kwargs['gradnorm'] if 'gradnorm' in kwargs else 1.0E-4
        optim_method = kwargs['optim_method'] if 'optim_method' in kwargs else 'BFGS'
        use_hessian = kwargs['hessian'] if 'hessian' in kwargs else False
        sites_to_output = set(self.obs_cons.site_code_to_dataset.keys())
        sites_to_output = sites_to_output.union(set(self.obs_cons.mbl_sites))
        sites_to_output = sites_to_output.union(set(self.obs_cons.noaa_observatories))

        # self.var4d_setup(obs_to_assim=obs_to_assim)

        # Print number of function and adjoint evaluations
        if not self.verbose:
            common_dict = dict(total=float('inf'), unit='it', bar_format='{l_bar}{bar}| {n_fmt} [{elapsed}, ' '{rate_fmt}{postfix}]')
            self.progress_bars['cost']  = tqdm(desc='Cost function evaluation', **common_dict)
            self.progress_bars['grad']  = tqdm(desc='Gradient evaluation', **common_dict)
            self.progress_bars['fwd']   = tqdm(desc='Forward transport', **common_dict)
            self.progress_bars['adj']   = tqdm(desc='Adjoint transport', **common_dict)
            self.progress_bars['hessp'] = tqdm(desc='Hessian product evaluation', **common_dict)

        self.optim_diags = { # after optimization, show the progress of the cost function and gradient norm
            'iter_grad': 0,
            'iter_cost': 0,
            'iter_hess': 0,
            'iter': 0,
            'cost_function': [],
            'gradient_norm': [],
            }

        # first, store the observations simulated with prior fluxes as model diagnostic
        obs_apri = self.forward_transport(self.state_prior)
        self.obs_cons.save_timeseries_by_site(self.obs_vec, obs_apri, self.obs_err, self.output_dir, sites_to_output, 'apri', False)

        # now optimize
        minimize_args = dict(
            method=optim_method, jac=self.calculate_gradient,
            options={'maxiter': max_iter, 'gtol': max_gradnorm, 'disp': self.verbose},
            )
        if use_hessian and optim_method in ['Newton-CG', 'trust-krylov', 'trust-ncg', 'trust-constr']:
            minimize_args['hessp'] = self.hessian_product
        if self.store_intermediate_states and optim_method not in ['TNC', 'SLSQP', 'COBYLA']:
            minimize_args['callback'] = self.loop_callback
        with Timer("Optimization finished in ", print=self.verbose):
            res = optimize.minimize(self.calculate_cost, self.state_prior_preco, **minimize_args)

        self.optim_result = res
        self.state_poste_preco = res.x
        self.state_poste = self.state_to_flux(res.x)
        obs_apos = self.forward_transport(self.state_poste)
        self.obs_cons.save_timeseries_by_site(self.obs_vec, obs_apos, self.obs_err, self.output_dir, sites_to_output, 'apos', True)

        self.var4d_done()

        if not self.verbose:
            for step, pbar in self.progress_bars.items():
                pbar.close()
            self.progress_bars.pop(step)

        print('End of 4DVAR loop')

    def gradient_test(self, init_cond, starting_alpha=1.0, alpha_step=0.1, max_iter=10, tolerance=1.0E-7):
        assert init_cond.lower() in ['ones', 'random'], 'Invalid init_cond specified'

        self.var4d_setup() # setup the preconditioner and prior covariance matrix
        self.optim_diags = { # after optimization, show the progress of the cost function and gradient norm
            'iter': 0,
            'cost_function': [],
            'gradient_norm': [],
            }

        if init_cond.lower() == 'ones':
            x0 = np.ones_like(self.state_prior_preco)
        elif init_cond.lower() == 'random':
            x0 = np.random.standard_normal(self.state_prior_preco.shape)

        # calculate cost at x0
        J0 = self.calculate_cost(x0)
        # calculate gradient at x0
        g0 = self.calculate_gradient(x0)
        # calculate DJ1
        DJ1 = np.matmul(g0.T, g0)

        # now start the gradient test loop
        converged = False
        alpha = starting_alpha
        current_iter = 1
        # write the header line
        print(
            'iter'.rjust(4) + 'alpha'.rjust(8) + 'J_bg'.rjust(15) + 'J_obs'.rjust(15) + 'J_tot'.rjust(15) +
            'DJ1'.rjust(15) + 'DJ2'.rjust(15) + 'DJ2/DJ1'.rjust(15) + '1-DJ2/DJ1'.rjust(15))

        while (current_iter <= max_iter) and (not converged):
            x = x0 - alpha*g0
            # calculate cost at x0-alpha*g0
            J = self.calculate_cost(x)
            DJ2 = (J0-J)/alpha
            print(
                "%4i"%current_iter + "%8.3g"%alpha + "%15.7e"%(0.5*sum(x0**2)) + "%15.7e"%(J-0.5*sum(x0**2)) + "%15.7e"%J +
                "%15.7e"%DJ1 + "%15.7e"%DJ2 + "%15.7g"%(DJ2/DJ1) + "%15.7e"%(1-DJ2/DJ1))

            if np.abs(1-DJ2/DJ1) < tolerance:
                converged = True
            else:
                alpha = alpha*alpha_step
                current_iter += 1

    def loop_callback(self, intermediate_result):
        self.interim_states['x'].append(intermediate_result.x)
        self.interim_states['fun'].append(intermediate_result.fun)

    def state_to_flux(self, state_vector):
        return self.state_prior + np.matmul(self.L, state_vector)

    def forward_transport(self, flux_vector):
        if 'fwd' in self.progress_bars:
            self.progress_bars['fwd'].update()
        return self.trans_op.transport(flux_vector, model='GC', add_bg=True)

    def calculate_mismatch(self, model_obs):
        return model_obs - self.obs_vec

    def adjoint_transport(self, obs_vector):
        if 'adj' in self.progress_bars:
            self.progress_bars['adj'].update()
        return np.matmul(self.trans_op.Jacobian['GC'].T, obs_vector)

    def hessian_product(self, p, x):
        # returns the product of the Hessian and vector x, without explicitly evaluating the Hessian
        # ignore p, the value of the state vector
        with Timer("Calculated product with Hessian in ", print=self.verbose) as tm:
            Lx = np.matmul(self.L, x)
            HLx = self.forward_transport(Lx)
            RinvHLx = HLx/(self.obs_err**2)
            HTRinvHLx = self.adjoint_transport(RinvHLx)
            LTHTRinvHLx = np.matmul(self.L.T, HTRinvHLx)

            self.optim_diags['iter_hess'] += 1
            tm['prefix'] = '[%i]'%self.optim_diags['iter_hess']

        if 'hessp' in self.progress_bars:
            self.progress_bars['hessp'].update()

        return LTHTRinvHLx

    def calculate_cost(self, state_vector):
        with Timer("Cost calculated in ", print=self.verbose) as tm:
            # convert to flux
            flux_vec = self.state_to_flux(state_vector)
            # simulate obs
            obs = self.forward_transport(flux_vec)
            # Hx-y
            self.mdm = self.calculate_mismatch(obs) # saving it in 'self' because we need this to calculate the gradient as well
            obs_cost = 0.5 * np.sum((self.mdm/self.obs_err)**2)
            bg_cost = 0.5 * np.sum(state_vector**2)
            total_cost = obs_cost + bg_cost
            self.optim_diags['cost_function'].append(total_cost)
            self.optim_diags['iter_cost'] += 1
            tm['prefix'] = '[%i]'%self.optim_diags['iter_cost']
            tm['postfix'] = '(J = %.4g)'%total_cost

        if 'cost' in self.progress_bars:
            self.progress_bars['cost'].update()

        return total_cost

    def calculate_gradient(self, state_vector):
        with Timer("Gradient calculated in ", print=self.verbose) as tm:
            # assume that self.mdm has already been calculated and is the most current
            if not 'mdm' in dir(self): # can happen e.g., for CG algorithms
                # convert to flux
                flux_vec = self.state_to_flux(state_vector)
                # simulate obs
                obs = self.forward_transport(flux_vec)
                self.mdm = self.calculate_mismatch(obs)

            adj_forcing = self.mdm/(self.obs_err**2)
            gradient = self.adjoint_transport(adj_forcing)
            gradient_preco = np.matmul(self.L.T, gradient)
            # delete self.mdm to prevent reusing old ones
            del self.mdm
            total_gradient = gradient_preco + state_vector
            self.optim_diags['iter'] += 1
            self.optim_diags['iter_grad'] += 1
            gradient_norm = np.sqrt(np.sum(total_gradient**2))
            self.optim_diags['gradient_norm'].append(gradient_norm)
            tm['prefix'] = '[%i]'%self.optim_diags['iter_grad']
            tm['postfix'] = u'(\u23b8\u2202J/\u2202\u03be\u23b9 = %.4g)'%gradient_norm

        if 'grad' in self.progress_bars:
            self.progress_bars['grad'].update()

        return total_gradient
