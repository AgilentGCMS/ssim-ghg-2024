from netCDF4 import Dataset
from matplotlib import pyplot as plt
import os, copy, calendar
from datetime import datetime, timedelta
from convert_jacobian import Paths
import numpy as np
from matplotlib.ticker import MaxNLocator

class Visualize(Paths):

    def __init__(self, project):
        super(Visualize, self).__init__()
        self.output_dir = os.path.join(self.output_root, project)
        self.plot_styles = {
            'obs': {'c': 'xkcd:puce', 'mec': 'xkcd:puce', 'mfc': 'xkcd:goldenrod', 'marker': 'o', 'ms': 5, 'mew': 0.5},
            'apri': {'c': 'xkcd:cerulean', 'mec': 'xkcd:cerulean', 'mfc': 'xkcd:powder blue', 'marker': 'd', 'ms': 5, 'mew': 0.5},
            'apos': {'c': 'xkcd:orange red', 'mec': 'xkcd:orange red', 'mfc': 'xkcd:baby pink', 'marker': 's', 'ms': 5, 'mew': 0.5},
            }
        self.legend_props = dict(fontsize=12, labelspacing=0.2, handlelength=1, handletextpad=0.25, borderpad=0.2)
        self.tick_font_property = dict(fontsize=12, family='Inconsolata')
        self.label_font_property = dict(fontsize=14, family='Inconsolata')

class Visualize_Fluxes(Visualize):

    def __init__(self, project):
        super(Visualize_Fluxes, self).__init__(project)

    def plot_region(self, region_names):
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

                time_vals = np.arange(len(month_times))

                plot_ax = ax_dict[region]
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

    def plot_site(self, site_codes):
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

            time_range = times.max()-times.min()
            xmin = times.min() - 0.025*time_range
            xmax = times.max() + 0.025*time_range

            plot_ax = ax_dict[site]['plot']
            diff_ax = ax_dict[site]['diff']

            plot_ax.plot(times, obs, label='Observed', ls='none', **self.plot_styles['obs'])
            plot_ax.plot(times, model_apri, label='Model (apri)', ls='none', **self.plot_styles['apri'])
            plot_ax.plot(times, model_apos, label='Model (apos)', ls='none', **self.plot_styles['apos'])

            leg = plot_ax.legend(loc='best', numpoints=2, **self.legend_props)
            leg.set_draggable(True)
            plt.setp(leg.texts, family=self.label_font_property['family'])

            diff_ax.plot(times, model_apri-obs, **self.plot_styles['apri'])
            diff_ax.plot(times, model_apos-obs, **self.plot_styles['apos'])

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
        fig.text(0.05/fig_width, (bpad+0.5*diff_height)/fig_height, 'Difference (ppm)', ha='left', va='center', rotation=90, **self.label_font_property)
