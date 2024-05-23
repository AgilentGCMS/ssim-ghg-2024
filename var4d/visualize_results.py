from netCDF4 import Dataset
from matplotlib import pyplot as plt
import os, copy, calendar, glob, tqdm
from datetime import datetime, timedelta
from convert_jacobian import Paths
import numpy as np
from matplotlib.ticker import MaxNLocator

class Visualize(Paths):

    def __init__(self, project):
        super(Visualize, self).__init__()
        self.output_dir = os.path.join(self.output_root, project)
        self.project = project
        self.plot_styles = {
            'obs': {'c': 'xkcd:puce', 'mec': 'xkcd:puce', 'mfc': 'xkcd:goldenrod', 'marker': 'o', 'ms': 5, 'mew': 0.5},
            'apri': {'c': 'xkcd:cerulean', 'mec': 'xkcd:cerulean', 'mfc': 'xkcd:powder blue', 'marker': 'd', 'ms': 5, 'mew': 0.5},
            'apos': {'c': 'xkcd:orange red', 'mec': 'xkcd:orange red', 'mfc': 'xkcd:baby pink', 'marker': 's', 'ms': 5, 'mew': 0.5},
            }
        self.legend_props = dict(fontsize=12, labelspacing=0.2, handlelength=1, handletextpad=0.25, borderpad=0.2)
        self.tick_font_property = dict(fontsize=12, family=self.figure_font)
        self.label_font_property = dict(fontsize=14, family=self.figure_font)
        self.error_dirs = {
            'only_noaa_observatories': 'only_noaa_observatories_mc/summary',
            }

class Monte_Carlo_avg(Paths):

    def __init__(self, project, glob_pattern='???'):
        super(Monte_Carlo_avg, self).__init__()
        self.output_dirs = glob.glob(os.path.join(self.output_root, project, glob_pattern))
        self.summary_dir = os.path.join(self.output_root, project, 'summary')

    def summarize_observations(self):
        if not os.path.isdir(self.summary_dir):
            os.makedirs(self.summary_dir)

        n_members = len(self.output_dirs)
        sample_files = glob.glob(os.path.join(self.output_dirs[0], 'timeseries_*.nc'))
        site_list = [os.path.basename(f).split('.')[0].split('_')[1] for f in sample_files]

        for site in site_list:
            output_file = os.path.join(self.summary_dir, 'timeseries_%s.nc'%site)
            comp_dict = dict(zlib=True, shuffle=True, complevel=5)
            with Dataset(output_file, 'w') as ofid:
                for i, dir_name in enumerate(tqdm.tqdm(self.output_dirs, desc='Summarizing observations for %s'%site)):
                    input_fname = os.path.join(dir_name, 'timeseries_%s.nc'%site)
                    with Dataset(input_fname, 'r') as ifid:
                        if i == 0:
                            # initialize the output file
                            ofid.createDimension('obs', len(ifid.dimensions['obs']))
                            ofid.createDimension('idate', 6)
                            ofid.createDimension('members', n_members)

                            for var_name in ['integer_times', 'decimal_times', 'observed', 'obs_error']:
                                var_in = ifid.variables[var_name]
                                v = ofid.createVariable(var_name, var_in.dtype, var_in.dimensions, **comp_dict)
                                v[:] = var_in[:]
                                for attr_name in var_in.ncattrs():
                                    setattr(v, attr_name, getattr(var_in, attr_name))

                            for var_name in ['modeled_apri', 'modeled_apos']:
                                var_in = ifid.variables[var_name]
                                v = ofid.createVariable(var_name, var_in.dtype, ('members',)+var_in.dimensions, **comp_dict)
                                for attr_name in var_in.ncattrs():
                                    setattr(v, attr_name, getattr(var_in, attr_name))

                        for var_name in ['modeled_apri', 'modeled_apos']:
                            ofid.variables[var_name][i] = ifid.variables[var_name][:]

                setattr(ofid, 'creation_date', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    def summarize_emissions(self):
        if not os.path.isdir(self.summary_dir):
            os.makedirs(self.summary_dir)
        output_file = os.path.join(self.summary_dir, 'optim_summary_spread.nc')

        var_list_to_write = ['x_final', 'j_final', 'prior_flux', 'poste_flux']
        n_members = len(self.output_dirs)

        comp_dict = dict(zlib=True, shuffle=True, complevel=5)
        with Dataset(output_file, 'w') as ofid:
            for i, dir_name in enumerate(tqdm.tqdm(self.output_dirs, desc='Summarizing emissions')):
                fname = os.path.join(dir_name, 'optim_summary.nc')
                with Dataset(fname, 'r') as ifid:
                    if i == 0:
                        # initialize the output file
                        ofid.createDimension("members", n_members)
                        ofid.createDimension("n_state", len(ifid.dimensions["n_state"]))
                        ofid.createDimension("n_region", len(ifid.dimensions['n_region']))
                        ofid.createDimension('n_month', len(ifid.dimensions['n_month']))

                        for var_name in var_list_to_write:
                            v_in = ifid.variables[var_name]
                            out_dims = ('members',) + v_in.dimensions
                            v_out = ofid.createVariable(var_name, v_in.dtype, out_dims, **comp_dict)
                            for attr_name in v_in.ncattrs():
                                setattr(v_out, attr_name, getattr(v_in, attr_name))

                    for var_name in var_list_to_write:
                        ofid.variables[var_name][i] = ifid.variables[var_name][:]

            setattr(ofid, 'creation_date', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

class Visualize_Fluxes(Visualize):

    def __init__(self, project):
        super(Visualize_Fluxes, self).__init__(project)

    def plot_region(self, region_names, plot_errs=False):
        if isinstance(region_names, str):
            region_names = [region_names]

        num_regions = len(region_names)
        ax_dict = dict.fromkeys(region_names)
        height = 3.0
        width = 4.0
        lpad = 0.6
        rpad = 0.1
        bpad = 0.5
        wd_pad = 0.4
        tpad = 0.3
        fig_width = lpad + num_regions*width + (num_regions-1)*wd_pad + rpad
        fig_height = bpad + height + tpad
        fig = plt.figure(figsize=(fig_width, fig_height))
        for i, region in enumerate(region_names):
            ax_dict[region] = plt.axes([(lpad+i*(width+wd_pad))/fig_width, bpad/fig_height, width/fig_width, height/fig_height])

        result_file = os.path.join(self.output_dir, 'optim_summary.nc')
        with Dataset(result_file, 'r') as fid:
            month_times = np.array([datetime(*d) for d in fid.variables['time_coordinate'][:]])
            n_month = len(month_times)
            region_names_in_file = [s.lower() for s in fid.variables['region_areas'].region_names.split(',')]
            region_areas = fid.variables['region_areas'][:] # m^2

            flux_conversion_factor = (12.01/44.01) * (1.0E-12 * 86400.0 * 365.25) # from Kg CO2/s to Pg C/year
            for iax,region in enumerate(region_names):
                region_index = region_names_in_file.index(region.lower())
                prior_flux = fid.variables['prior_flux'][region_index] * region_areas[region_index] * flux_conversion_factor # PgC/year
                poste_flux = fid.variables['poste_flux'][region_index] * region_areas[region_index] * flux_conversion_factor # PgC/year
                true_flux  = fid.variables['true_flux'][region_index]  * region_areas[region_index] * flux_conversion_factor # PgC/year

                if plot_errs and self.project in self.error_dirs:
                    error_file = os.path.join(self.output_root, self.error_dirs[self.project], 'optim_summary_spread.nc')
                    with Dataset(error_file, 'r') as e_fid:
                        prior_ensemble = e_fid.variables['prior_flux'][:,region_index]
                        poste_ensemble = e_fid.variables['poste_flux'][:,region_index]
                    prior_err = np.std(prior_ensemble, axis=0) * region_areas[region_index] * flux_conversion_factor
                    poste_err = np.std(poste_ensemble, axis=0) * region_areas[region_index] * flux_conversion_factor
                else:
                    plot_errs = False

                time_vals = np.arange(len(month_times))

                plot_ax = ax_dict[region]
                if plot_errs:
                    plot_ax.fill_between(time_vals, prior_flux-prior_err, prior_flux+prior_err,
                        zorder=0, fc=self.plot_styles['apri']['mfc'], ec=self.plot_styles['apri']['mec'], lw=0.5, alpha=0.6)
                    plot_ax.fill_between(time_vals, poste_flux-poste_err, poste_flux+poste_err,
                        zorder=0, fc=self.plot_styles['apos']['mfc'], ec=self.plot_styles['apos']['mec'], lw=0.5, alpha=0.6)
                plot_ax.plot(time_vals, true_flux, label='True', ls='-', lw=1, **self.plot_styles['obs'])
                plot_ax.plot(time_vals, prior_flux, label='Prior', ls='--', lw=1, **self.plot_styles['apri'])
                plot_ax.plot(time_vals, poste_flux, label='Poste', ls='-', lw=1, **self.plot_styles['apos'])

                # plot the totals
                indices_for_total = np.array([i for i,d in enumerate(month_times) if d.year == 2015])
                month_lengths = np.array([calendar.monthrange(2015,m)[1] for m in range(1,13)])
                tot_true = np.average(true_flux[indices_for_total], weights=month_lengths)
                tot_prior = np.average(prior_flux[indices_for_total], weights=month_lengths)
                tot_poste = np.average(poste_flux[indices_for_total], weights=month_lengths)

                plot_ax.plot(time_vals[-1]+1, tot_true, **self.plot_styles['obs'])
                plot_ax.plot(time_vals[-1]+1, tot_prior, **self.plot_styles['apri'])
                plot_ax.plot(time_vals[-1]+1, tot_poste, **self.plot_styles['apos'])

                leg = plot_ax.legend(loc='best', numpoints=1, **self.legend_props)
                leg.set_draggable(True)
                plt.setp(leg.texts, family=self.label_font_property['family'])

                xtick_locs = time_vals[::3]
                xtick_labels = [d.strftime('%b\n%Y') for d in month_times[::3]]
                # add the "Total" label
                xtick_locs = list(xtick_locs) + [time_vals[-1]+1]
                xtick_labels.append('2015\nTotal')

                plot_ax.set_xticks(xtick_locs)
                plot_ax.yaxis.set_major_locator(MaxNLocator(integer=True))

                plot_ax.set_xticklabels(xtick_labels, **self.tick_font_property)
                plt.setp(plot_ax.get_yticklabels(), **self.tick_font_property)

                plot_ax.grid(True, ls='--')
                plot_ax.axvspan(time_vals[-1]+0.5, time_vals[-1]+1.5, ec=None, fc='0.75')

                plot_ax.set_xlim(time_vals[0]-0.5, time_vals[-1]+1.5)

                fig.text((lpad+0.5*width+iax*(width+wd_pad))/fig_width, 1.0-0.5*tpad/fig_height, region.title(), ha='center', va='center', **self.label_font_property)

        fig.text(0.05/fig_width, (bpad+0.5*height)/fig_height, u'CO\u2082 flux (PgC/year)', ha='left', va='center', rotation=90, **self.label_font_property)

class Visualize_Obs(Visualize):

    def __init__(self, project):
        super(Visualize_Obs, self).__init__(project)

    def plot_site(self, site_codes, plot_errs=False):
        if isinstance(site_codes, str):
            site_codes = [site_codes]

        num_sites = len(site_codes)
        ax_dict = dict.fromkeys(site_codes)
        height = 3.0
        diff_height = 1.5
        width = 5.0
        lpad = 0.7
        rpad = 0.1
        bpad = 0.1
        ht_pad = 0.6
        wd_pad = 0.5
        tpad = 0.3
        fig_width = lpad + num_sites*width + (num_sites-1)*wd_pad + rpad
        fig_height = bpad + diff_height + ht_pad + height + tpad
        fig = plt.figure(figsize=(fig_width, fig_height))
        for i, site in enumerate(site_codes):
            ax_dict[site] = {
                'plot': plt.axes([(lpad+i*(width+wd_pad))/fig_width, (bpad+diff_height+ht_pad)/fig_height, width/fig_width, height/fig_height]),
                'diff': plt.axes([(lpad+i*(width+wd_pad))/fig_width, bpad/fig_height, width/fig_width, diff_height/fig_height]),
            }

        for iax,site in enumerate(site_codes):
            site_file = os.path.join(self.output_dir, 'timeseries_%s.nc'%site)
            with Dataset(site_file, 'r') as fid:
                times = fid.variables['decimal_times'][:]
                obs = fid.variables['observed'][:]
                model_apri = fid.variables['modeled_apri'][:]
                model_apos = fid.variables['modeled_apos'][:]

            if plot_errs and self.project in self.error_dirs:
                error_file = os.path.join(self.output_root, self.error_dirs[self.project], 'timeseries_%s.nc'%site)
                with Dataset(error_file, 'r') as fid:
                    apri_err = fid.variables['modeled_apri'][:].std(axis=0)
                    apos_err = fid.variables['modeled_apos'][:].std(axis=0)
            else:
                plot_errs = False

            time_range = times.max()-times.min()
            xmin = times.min() - 0.025*time_range
            xmax = times.max() + 0.025*time_range

            plot_ax = ax_dict[site]['plot']
            diff_ax = ax_dict[site]['diff']

            if plot_errs:
                plot_ax.fill_between(times, model_apri-apri_err, model_apri+apri_err,
                    zorder=0, fc=self.plot_styles['apri']['mfc'], ec=self.plot_styles['apri']['mec'], lw=0.25, alpha=0.6)
                plot_ax.fill_between(times, model_apos-apos_err, model_apos+apos_err,
                    zorder=0, fc=self.plot_styles['apos']['mfc'], ec=self.plot_styles['apos']['mec'], lw=0.25, alpha=0.6)
            plot_ax.plot(times, obs, label='Observed', ls='none', **self.plot_styles['obs'])
            plot_ax.plot(times, model_apri, label='Model (apri)', ls='none', **self.plot_styles['apri'])
            plot_ax.plot(times, model_apos, label='Model (apos)', ls='none', **self.plot_styles['apos'])

            leg = plot_ax.legend(loc='best', numpoints=2, **self.legend_props)
            leg.set_draggable(True)
            plt.setp(leg.texts, family=self.label_font_property['family'])

            diff_ax.plot(times, model_apri-obs, ls='none', **self.plot_styles['apri'])
            diff_ax.plot(times, model_apos-obs, ls='none', **self.plot_styles['apos'])

            plot_ax.set_xlim(xmin, xmax)
            diff_ax.set_xlim(xmin, xmax)
            plt.setp(plot_ax.get_yticklabels(), **self.tick_font_property)
            plt.setp(diff_ax.get_yticklabels(), **self.tick_font_property)
            plt.setp(plot_ax.get_xticklabels(), **self.tick_font_property)

            plot_ax.locator_params(axis='x', nbins=5)
            diff_ax.locator_params(axis='x', nbins=5)
            plot_ax.grid(True, ls='--')
            diff_ax.grid(True, ls='--')

            fig.text((lpad+0.5*width+iax*(width+wd_pad))/fig_width, (bpad+diff_height+0.05)/fig_height, 'Decimal year', ha='center', va='bottom', size=14)
            fig.text((lpad+0.5*width+iax*(width+wd_pad))/fig_width, 1.0-0.5*tpad/fig_height, site.upper(), ha='center', va='center', size=14)

        fig.text(0.05/fig_width, (bpad+diff_height+ht_pad+0.5*height)/fig_height, u'CO\u2082 (ppm)', ha='left', va='center', rotation=90, **self.label_font_property)
        fig.text(0.05/fig_width, (bpad+0.5*diff_height)/fig_height, u'Model\u2009\u2212\u2009Obs (ppm)', ha='left', va='center', rotation=90, **self.label_font_property)
