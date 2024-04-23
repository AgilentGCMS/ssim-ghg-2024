from netCDF4 import Dataset
from matplotlib import pyplot as plt
import os
from datetime import datetime, timedelta
from convert_jacobian import Paths
import numpy as np

class Visualize_Obs(Paths):

    def __init__(self, project):
        super(Visualize_Obs, self).__init__()
        self.output_dir = os.path.join(self.output_root, project)
        self.plot_styles = {
            'obs': {'c': 'xkcd:puce', 'mec': 'xkcd:puce', 'mfc': 'xkcd:goldenrod', 'marker': 'o', 'ls': 'none', 'ms': 5, 'mew': 0.5},
            'apri': {'c': 'xkcd:cerulean', 'mec': 'xkcd:cerulean', 'mfc': 'xkcd:powder blue', 'marker': 'd', 'ls': 'none', 'ms': 5, 'mew': 0.5},
            'apos': {'c': 'xkcd:orange red', 'mec': 'xkcd:orange red', 'mfc': 'xkcd:baby pink', 'marker': 's', 'ls': 'none', 'ms': 5, 'mew': 0.5},
            }

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
        tpad = 0.1
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

            plot_ax.plot(times, obs, label='Observed', **self.plot_styles['obs'])
            plot_ax.plot(times, model_apri, label='Model (apri)', **self.plot_styles['apri'])
            plot_ax.plot(times, model_apos, label='Model (apos)', **self.plot_styles['apos'])

            leg = plot_ax.legend(loc='best', numpoints=2, fontsize=14, labelspacing=0.2, handlelength=1)
            leg.set_draggable(True)

            diff_ax.plot(times, model_apri-obs, **self.plot_styles['apri'])
            diff_ax.plot(times, model_apos-obs, **self.plot_styles['apos'])

            plot_ax.set_xlim(xmin, xmax)
            diff_ax.set_xlim(xmin, xmax)
            plt.setp(plot_ax.get_yticklabels(), size=12)
            plt.setp(diff_ax.get_yticklabels(), size=12)
            plt.setp(plot_ax.get_xticklabels(), size=12)

            plot_ax.locator_params(axis='x', nbins=5)
            diff_ax.locator_params(axis='x', nbins=5)
            plot_ax.grid(True, ls='--')
            diff_ax.grid(True, ls='--')

            fig.text((lpad+0.5*width+iax*(width+wd_pad))/fig_width, (bpad+diff_height+0.05)/fig_height, 'Decimal year', ha='center', va='bottom', size=14)
            fig.text((lpad+(iax+1)*width+iax*wd_pad-0.1)/fig_width, 1.-(tpad+0.1)/fig_height, site.upper(), ha='right', va='top', size=14, bbox={'boxstyle':'round', 'facecolor':'wheat', 'alpha':0.5})
            
        fig.text(0.05/fig_width, (bpad+diff_height+ht_pad+0.5*height)/fig_height, u'CO\u2082 (ppm)', ha='left', va='center', size=14, rotation=90)
        fig.text(0.05/fig_width, (bpad+0.5*diff_height)/fig_height, 'Difference (ppm)', ha='left', va='center', size=14, rotation=90)
