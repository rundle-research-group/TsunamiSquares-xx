#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pylab
from matplotlib import pyplot as plt
import matplotlib.colors as mcolor
import matplotlib.animation as manimation
import matplotlib.colorbar as mcolorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.font_manager as mfont
import scipy as sp
import argparse
from geographiclib.geodesic import Geodesic as geo
from netCDF4 import Dataset
import os
import read_ETOPO1
import pandas as pd
import json
#import quakelib

# --------------------------------------------------------------------------------

class simAnalyzer:
    def __init__(self, sim_file_path):
        self.save_file_prefix = os.path.splitext(sim_file_path)[0]
        #sim_data = np.genfromtxt(sim_file, dtype=[('time','f8'),('lat','f8'),('lon','f8'), ('z','f8'), ('alt','f8')])
        self.sim_data = Dataset(sim_file_path, 'r', format='NETCDF4')
        
        # These arrays shouldn't be too big, so go ahead and load them into memory as numpy arrays    
        self.times = np.array(self.sim_data.variables['time'])
        self.lons = np.array(self.sim_data.variables['longitude'])
        self.lats = np.array(self.sim_data.variables['latitude'])
        
        #TODO: check for simulation that wraps around int date line.
        self.minlon, self.maxlon, self.meanlon = self.lons.min(), self.lons.max(), self.lons.mean()
        self.minlat, self.maxlat, self.meanlat = self.lats.min(), self.lats.max(), self.lats.mean()
        self.dlon = abs(self.lons[1]-self.lons[0])
        self.dlat = abs(self.lats[1]-self.lats[0])      
        
        height_ncVar = self.sim_data.variables['height']
        alt_ncVar = self.sim_data.variables['altitude']
        
        # Calculate locations where we simulated ocean runup
        ever_has_water = np.any(height_ncVar, axis=0)
        self.sim_inundation_array = np.ma.masked_where(alt_ncVar[0] < 0, ever_has_water)
        
        
    def make_grid_animation(self, FPS, DPI, zminmax=None, doBasemap=False):
        
        save_file = self.save_file_prefix+"_grid.mp4"
        
        # Keep the data from each time step in netCDF variable form, and slice into it as needed
        level_ncVar = self.sim_data.variables['level']
        height_ncVar = self.sim_data.variables['height']
        alt_ncVar = self.sim_data.variables['altitude']
        
        # Get ranges
        N_STEP = len(self.times)
        z_min =  np.inf
        z_max = -np.inf
        z_avs = []
        for i, levelstep in enumerate(level_ncVar):
            masked_data = np.ma.masked_where(height_ncVar[i] == 0.0000, levelstep)  
            z_min = min(masked_data.min(), z_min)
            z_max = max(masked_data.max(), z_max)
            z_avs.append(masked_data.mean())
        
        z_max = np.max(np.ma.masked_where(height_ncVar[0] == 0.0000, level_ncVar[0]))    
        
        print("min: {}, max: {}, av: {}".format(z_min, z_max, np.array(z_avs).mean()))
        if(zminmax != None): z_min,z_max = zminmax
    
        # Initialize movie writing stuff
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='TsunamiSquares', artist='Matplotlib', comment='Animation')
        writer = FFMpegWriter(fps=FPS, metadata=metadata, bitrate=1000)
    
        # Initialize the frame and axes
        fig = plt.figure()
        
        if not doBasemap:
            ax = fig.add_subplot(111)
            plt.xlim(self.minlon, self.maxlon)
            plt.ylim(self.minlat, self.maxlat)
            #ax.get_xaxis().get_major_formatter().set_useOffset(False)
            #ax.get_yaxis().get_major_formatter().set_useOffset(False)
        else:
            m = Basemap(projection='cyl', llcrnrlat=self.minlat, urcrnrlat=self.maxlat,
                        llcrnrlon=self.minlon, urcrnrlon=self.maxlon, lat_0=self.meanlat, lon_0=self.meanlon, resolution='h')
            m.drawmeridians(np.linspace(self.minlon, self.maxlon, num=5.0), labels=[0,0,0,1], linewidth=0)
            m.drawparallels(np.linspace(self.minlat, self.maxlat, num=5.0), labels=[1,0,0,0], linewidth=0)
            m.drawcoastlines(linewidth=0.5)
            m.ax = fig.add_subplot(111)
            ax = m.ax
        
        # Colorbar
        cmap = plt.get_cmap('Blues_r')
        landcolor = 'orange'#'black'#'#FFFFCC'
        cmap.set_bad(landcolor, 1.0)
        
        norm = mcolor.Normalize(vmin=z_min, vmax=z_max)
        divider = make_axes_locatable(ax)
        cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
        
        cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)
        framelabelfont = mfont.FontProperties(style='normal', variant='normal', size=14)
        plt.figtext(0.95, 0.7, r'water altitude $[m]$', rotation='vertical', fontproperties=framelabelfont)
    
        surface = None
        with writer.saving(fig, save_file, DPI):
            for index in range(int(N_STEP)):
                # Get the subset of data corresponding to current time
                this_level  = level_ncVar[index]
                this_height = height_ncVar[index]
                this_alt    = alt_ncVar[index]
                time        = self.times[index]
                
                # Masked array via conditional, don't color the land unless it has water on it
                masked_data = np.ma.masked_where(this_height == 0.0000, this_level)            
                
                print("step: {}  time: {}".format(index, time))
    
                # Plot the surface for this time step
                if surface is None:
                    ax.imshow(masked_data, cmap=cmap,origin='lower', norm=norm, extent=[self.minlon, self.maxlon, self.maxlat, self.minlat], interpolation='none')
                else:
                    surface.set_data(masked_data)
                    
                # Text box with the time
                plt.figtext(0.02, 0.5, 'Time: {:02d}:{:02d}'.format(int(time/60), int(time%60)), bbox={'facecolor':'yellow', 'pad':5})
                
                writer.grab_frame()
    
    
    def make_crosssection_animation(self, FPS, DPI):
        
        save_file = self.save_file_prefix+"_crosssection.mp4"
        
        #sim_data is expected to be a netcdf dataset 
        lons = np.array(self.sim_data.variables['longitude'])
        lats = np.array(self.sim_data.variables['latitude'])
        
        # But keep the data from each time step in netCDF variable form, and slice into it as needed
        level_ncVar = self.sim_data.variables['level']
        height_ncVar = self.sim_data.variables['height']
        alt_ncVar = self.sim_data.variables['altitude']
        
        # Get ranges
        N_STEP = len(self.times)
        z_min =  np.inf
        z_max = -np.inf
        z_avs = []
        for i, levelstep in enumerate(level_ncVar):
            masked_data = np.ma.masked_where(height_ncVar[i] == 0.0000, levelstep)  
            z_min = min(masked_data.min(), z_min)
            z_max = max(masked_data.max(), z_max)
            z_avs.append(masked_data.mean())
        
        print("min: {}, max: {}, av: {}".format(z_min, z_max, np.array(z_avs).mean()))
        
    
        # Initialize movie writing stuff
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='TsunamiSquares', artist='Matplotlib', comment='Animation')
        writer = FFMpegWriter(fps=FPS, metadata=metadata, bitrate=1000)
    
        # Initialize the frame and axes
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.xlim(self.minlon, self.maxlon)
        plt.ylim(z_min, z_max)
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        divider = make_axes_locatable(ax)        
        
        midlat_index = int(len(self.lats)/2)
            
        # Make array of distances from mean lat/lon
        dist_array = []
        geod = geo(6371000, 0)
        for lon in self.lons:
            distance = geod.Inverse(self.lats[midlat_index], lon, self.meanlat, self.meanlon)['s12']
            dist_array.append(abs(distance))
        dist_array = np.array(dist_array)
    
        sim_plot = ax.plot([], [], 'b')[0]
        analytic_plot = ax.plot([], [], 'r')[0]
        
        with writer.saving(fig, save_file, DPI):
            for index in range(int(N_STEP)):
                time   = self.times[index]
                print("step: {}  time: {}".format(index, time))
                
                analytic_Z = []
                for dist in dist_array:
                    analytic_Z.append(analyticGaussPile(dist, time, 10, 5000, 1000))
                
                # Plot the cross section for this time step
                sim_plot.set_data(self.lons, level_ncVar[index][midlat_index])
                analytic_plot.set_data(self.lons, analytic_Z)
                
                # Text box with the time
                plt.figtext(0.02, 0.5, 'Time: {:02f}:{:02f}'.format(int(time)/60, int(time)%60), bbox={'facecolor':'yellow', 'pad':5})
                    
                writer.grab_frame()


    def load_historical_runup_CSV(self, csv_file_path, year, month, day):
        inund_df = pd.read_csv(csv_file_path, sep='\t')

        inund_df = inund_df[inund_df.LATITUDE.notnull()]
        inund_df = inund_df[inund_df.LONGITUDE.notnull()]
        inund_df = inund_df[inund_df.HORIZONTAL_INUNDATION.notnull()]
        
        self.inund_df = inund_df[((inund_df.YEAR == year) & (inund_df.MONTH == month) & (inund_df.DAY == day) &
                                 (inund_df.LATITUDE >= self.minlat) & (inund_df.LATITUDE <= self.maxlat) &
                                 (inund_df.LONGITUDE >= self.minlon) & (inund_df.LONGITUDE <= self.maxlon) &
                                 (inund_df.HORIZONTAL_INUNDATION > 0))]
        
        # Now create an array to match the simulated inundation array showing observed inundations
        obs_histogram, xedge, yedge = np.histogram2d(inund_df.LONGITUDE, inund_df.LATITUDE, bins=[len(self.lons), len(self.lats)], 
                                                          range=[[self.minlon-self.dlon/2.0, self.maxlon+self.dlon/2.0], [self.minlat-self.dlat/2.0, self.maxlat+self.dlat/2.0]])
        obs_histogram = np.flipud(obs_histogram.T)
        alt_ncVar = self.sim_data.variables['altitude']
        self.obs_inundation_array = np.ma.masked_where(alt_ncVar[0] < 0, (obs_histogram>0))

        

    def compare_sim_and_obs_runup(self):
        
        all_inundation_results = np.zeros_like(self.sim_inundation_array.astype(int))
        all_inundation_results[np.logical_and(self.sim_inundation_array==0, self.obs_inundation_array==1)] = 1
        all_inundation_results[np.logical_and(self.sim_inundation_array==1, self.obs_inundation_array==0)] = 2
        all_inundation_results[np.logical_and(self.sim_inundation_array==1, self.obs_inundation_array==1)] = 3
        
        alt_ncVar = self.sim_data.variables['altitude']
        self.all_inundation_results = np.ma.masked_where(alt_ncVar[0]<0, all_inundation_results)
        
        
        plt.close(1)
        fig, ax = plt.subplots(num=1)
        m = Basemap(projection='cyl',llcrnrlat=self.minlat, urcrnrlat=self.maxlat,
                    llcrnrlon=self.minlon, urcrnrlon=self.maxlon, lat_0=self.meanlat, lon_0=self.meanlon, resolution='i')
        m.drawmeridians(np.linspace(self.minlon, self.maxlon, num=5.0), labels=[0,0,0,1], linewidth=0)
        m.drawparallels(np.linspace(self.minlat, self.maxlat, num=5.0), labels=[1,0,0,0], linewidth=0)
        #m.drawcoastlines(linewidth=0.5)

        cm = mcolor.LinearSegmentedColormap.from_list('custom_cmap', ['gray', 'maroon', 'blue', 'lime'], N=4)
        
        map_ax = m.imshow(self.all_inundation_results, origin='upper', extent=[self.minlon, self.maxlon, self.minlat, self.maxlat], interpolation='nearest', cmap=cm)
        
        cbar = fig.colorbar(map_ax, ticks=[3/8., 9/8., 15/8., 21/8.])
        cbar.ax.set_yticklabels(['Dry', 'Miss', 'False\nAlarm', 'Success'])
        
        plt.savefig(self.save_file_prefix+'_inundation.png',dpi=100)
        
        

def analyticGaussPileIntegrand(k, r, t, Dc, Rc, depth):
    dispersion = np.sqrt(9.80665*k*np.tanh(k*depth))
    #dispersion = k*np.sqrt(9.80665*depth)
    return Dc*Rc**2*k/2*np.cos(dispersion*t)*sp.special.jv(0, k*r)*np.exp(-(k*Rc/2)**2)


def analyticGaussPile(r, t, Dc, Rc, depth):
    return sp.integrate.quad(analyticGaussPileIntegrand, 0, 1e3, args=(r, t, Dc, Rc, depth), points=[0, 2e-3])[0]
#    k = np.linspace(0, 2e-3, 1e4)
#    sumd = np.sum(analyticGaussIntegrand(k, r, t, Dc, Rc, depth))
#    return np.diff(k)[0]*sumd


def plot_eq_displacements(disp_file):

    save_file = os.path.splitext(disp_file)[0] + "_disp_z.png"    
    
    # Read displacement data
    disp_data = np.genfromtxt(LLD_FILE, dtype=[('lat','f8'),('lon','f8'), ('z','f8')],skip_header=3)

    # Data ranges
    lon_min,lon_max = disp_data['lon'].min(),disp_data['lon'].max()
    lat_min,lat_max = disp_data['lat'].min(),disp_data['lat'].max()
    mean_lat = 0.5*(lat_min + lat_max)
    mean_lon = 0.5*(lon_min + lon_max)
    lon_range = lon_max - lon_min
    lat_range = lat_max - lat_min
    z_min,z_max = disp_data['z'].min(),disp_data['z'].max()
    z_lim = max(np.abs(z_min),np.abs(z_max))
    cmap = plt.get_cmap('seismic')

    LEVELS = np.concatenate((-1*np.linspace(0.01, z_lim, 6)[::-1], np.linspace(0.01, z_lim, 6)))    
    
    norm = mcolor.Normalize(vmin=-z_lim, vmax=z_lim)
    interp = 'cubic'
    landcolor = '#FFFFCC'
    framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14)

    # Initialize the frame and axes
    fig = plt.figure()
    
    m = Basemap(projection='cyl',llcrnrlat=lat_min, urcrnrlat=lat_max,
                llcrnrlon=lon_min, urcrnrlon=lon_max, lat_0=mean_lat, lon_0=mean_lon, resolution='h')
    m.ax = fig.add_subplot(111)
    
    m.drawmeridians(np.linspace(lon_min,lon_max,num=5.0),labels=[0,0,0,1], linewidth=0)
    m.drawparallels(np.linspace(lat_min,lat_max,num=5.0),labels=[1,0,0,0], linewidth=0)
    m.drawcoastlines(linewidth=0.5)
    m.fillcontinents(color=landcolor, zorder=0)

    # Colorbar
    divider = make_axes_locatable(m.ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
    plt.figtext(0.96, 0.7, r'displacement $[m]$', rotation='vertical', fontproperties=framelabelfont)
    cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)

    # Reshape into matrices
    Ncols = len(np.unique(disp_data['lon']))
    Nrows = len(np.unique(disp_data['lat']))
    
    X = disp_data['lon'].reshape(Nrows, Ncols)
    Y = disp_data['lat'].reshape(Nrows, Ncols)
    Z = disp_data['z'].reshape(Nrows, Ncols)
    
    # Masked array via conditional, don't color the land unless it has water on it
    zero_below = int(len(LEVELS)/2)-1
    zero_above = zero_below+1
    masked_data = np.ma.masked_where(np.logical_and(np.array(Z <= LEVELS[zero_above]),np.array(Z >= LEVELS[zero_below])), Z)
    
    # Set masked pixels to the land color
    cmap.set_bad(landcolor, 0.0)  # set alpha=0.0 for transparent
    
    # Plot the contours
    m.contourf(X, Y, masked_data, LEVELS, cmap=cmap, norm=norm, extend='both', zorder=1)

    plt.savefig(save_file,dpi=100)
    print("Saved to "+save_file)


def plot_eq_disps_horiz(disp_file):
    # Read displacement data
    disp_data = np.genfromtxt(disp_file, dtype=[('lon','f8'), ('lat','f8'), ('z','f8'), ('eU','f8'), ('nV','f8')])
    
    save_file_prefix = os.path.splitext(disp_file)[0]+"_disp"    
    
    # Data ranges
    lon_min,lon_max = disp_data['lon'].min(),disp_data['lon'].max()
    lat_min,lat_max = disp_data['lat'].min(),disp_data['lat'].max()
    mean_lat = 0.5*(lat_min + lat_max)
    mean_lon = 0.5*(lon_min + lon_max)
    lon_range = lon_max - lon_min
    lat_range = lat_max - lat_min
    
    # Reshape into matrices
    Ncols = len(np.unique(disp_data['lon']))
    Nrows = len(np.unique(disp_data['lat']))
    X = disp_data['lon'].reshape(Nrows, Ncols)
    Y = disp_data['lat'].reshape(Nrows, Ncols)
    Z = disp_data['z'].reshape(Nrows, Ncols)
    eU = disp_data['eU'].reshape(Nrows, Ncols)
    nV = disp_data['nV'].reshape(Nrows, Ncols)
    
    cmap = plt.get_cmap('seismic')
    
    z_min,z_max = disp_data['z'].min(),disp_data['z'].max()
    z_lim = max(np.abs(z_min),np.abs(z_max))
    normz = mcolor.Normalize(vmin=-z_lim, vmax=z_lim)
    #LEVELSz = np.concatenate((-1*np.logspace(-3, np.log10(z_lim), 6)[::-1], np.logspace(-3, np.log10(z_lim), 6)))
    LEVELSz = np.concatenate((-1*np.linspace(0.01, z_lim, 6)[::-1], np.linspace(0.01, z_lim, 6)))
    
    e_min,e_max = disp_data['eU'].min(),disp_data['eU'].max()
    e_lim = max(np.abs(e_min),np.abs(e_max))
    norme = mcolor.Normalize(vmin=-e_lim, vmax=e_lim)
    #LEVELSe = np.concatenate((-1*np.logspace(-3, np.log10(e_lim), 6)[::-1], np.logspace(-3, np.log10(e_lim), 6)))
    LEVELSe = np.concatenate((-1*np.linspace(0.01, e_lim, 6)[::-1], np.linspace(0.01, e_lim, 6)))
    
    n_min,n_max = disp_data['nV'].min(),disp_data['nV'].max()
    n_lim = max(np.abs(n_min),np.abs(n_max))
    normn = mcolor.Normalize(vmin=-n_lim, vmax=n_lim)
    #LEVELSn = np.concatenate((-1*np.logspace(-3, np.log10(n_lim), 6)[::-1], np.logspace(-3, np.log10(n_lim), 6)))
    LEVELSn = np.concatenate((-1*np.linspace(0.01, n_lim, 6)[::-1], np.linspace(0.01, n_lim, 6)))
    
    interp = 'cubic'
    landcolor = '#FFFFCC'
    framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14)

    # Initialize the frame and axes
    fig = plt.figure(0)    
    m = Basemap(projection='cyl',llcrnrlat=lat_min, urcrnrlat=lat_max,
                llcrnrlon=lon_min, urcrnrlon=lon_max, lat_0=mean_lat, lon_0=mean_lon, resolution='h')
    m.ax = fig.add_subplot(111)
    m.drawmeridians(np.linspace(lon_min,lon_max,num=5.0),labels=[0,0,0,1], linewidth=0)
    m.drawparallels(np.linspace(lat_min,lat_max,num=5.0),labels=[1,0,0,0], linewidth=0)
    m.drawcoastlines(linewidth=0.5)
    m.fillcontinents(color=landcolor, zorder=0)

    # Colorbar
    divider = make_axes_locatable(m.ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
    plt.figtext(0.96, 0.7, r'Vertical disp $[m]$', rotation='vertical', fontproperties=framelabelfont)
    cbz = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=normz)
    
    # Masked array via conditional, don't color the land unless it has water on it
    zero_below = int(len(LEVELSz)/2)-1
    zero_above = zero_below+1
    masked_data = np.ma.masked_where(np.logical_and(np.array(Z <= LEVELSz[zero_above]),np.array(Z >= LEVELSz[zero_below])), Z)
    
    # Set masked pixels to the land color
    cmap.set_bad(landcolor, 0.0)  # set alpha=0.0 for transparent
    
    # Plot the contours
    m.contourf(X, Y, masked_data, LEVELSz, cmap=cmap, norm=normz, extend='both', zorder=1)

    plt.savefig(save_file_prefix+'_z.png',dpi=100)
    
    # Initialize the frame and axes
    fig = plt.figure(1)    
    m = Basemap(projection='cyl',llcrnrlat=lat_min, urcrnrlat=lat_max,
                llcrnrlon=lon_min, urcrnrlon=lon_max, lat_0=mean_lat, lon_0=mean_lon, resolution='h')
    m.ax = fig.add_subplot(111)
    m.drawmeridians(np.linspace(lon_min,lon_max,num=5.0),labels=[0,0,0,1], linewidth=0)
    m.drawparallels(np.linspace(lat_min,lat_max,num=5.0),labels=[1,0,0,0], linewidth=0)
    m.drawcoastlines(linewidth=0.5)
    m.fillcontinents(color=landcolor, zorder=0)

    # Colorbar
    divider = make_axes_locatable(m.ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
    plt.figtext(0.96, 0.7, r'East disp $[m]$', rotation='vertical', fontproperties=framelabelfont)
    cbe = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norme)
    
    # Masked array via conditional, don't color the land unless it has water on it
    zero_below = int(len(LEVELSe)/2)-1
    zero_above = zero_below+1
    masked_data = np.ma.masked_where(np.logical_and(np.array(eU <= LEVELSe[zero_above]),np.array(eU >= LEVELSe[zero_below])), eU)
    
    # Set masked pixels to the land color
    cmap.set_bad(landcolor, 0.0)  # set alpha=0.0 for transparent
    
    # Plot the contours
    m.contourf(X, Y, masked_data, LEVELSe, cmap=cmap, norm=norme, extend='both', zorder=1)

    plt.savefig(save_file_prefix+'_e.png',dpi=100)

    # Initialize the frame and axes
    fig = plt.figure(2)    
    m = Basemap(projection='cyl',llcrnrlat=lat_min, urcrnrlat=lat_max,
                llcrnrlon=lon_min, urcrnrlon=lon_max, lat_0=mean_lat, lon_0=mean_lon, resolution='h')
    m.ax = fig.add_subplot(111)
    m.drawmeridians(np.linspace(lon_min,lon_max,num=5.0),labels=[0,0,0,1], linewidth=0)
    m.drawparallels(np.linspace(lat_min,lat_max,num=5.0),labels=[1,0,0,0], linewidth=0)
    m.drawcoastlines(linewidth=0.5)
    m.fillcontinents(color=landcolor, zorder=0)

    # Colorbar
    divider = make_axes_locatable(m.ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
    plt.figtext(0.96, 0.7, r'North disp $[m]$', rotation='vertical', fontproperties=framelabelfont)
    cbn = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=normn)
    
    # Masked array via conditional, don't color the land unless it has water on it
    zero_below = int(len(LEVELSn)/2)-1
    zero_above = zero_below+1
    masked_data = np.ma.masked_where(np.logical_and(np.array(nV <= LEVELSn[zero_above]),np.array(nV >= LEVELSn[zero_below])), nV)
    
    # Set masked pixels to the land color
    cmap.set_bad(landcolor, 0.0)  # set alpha=0.0 for transparent
    
    # Plot the contours
    m.contourf(X, Y, masked_data, LEVELSn, cmap=cmap, norm=normn, extend='both', zorder=1)

    plt.savefig(save_file_prefix+'_n.png',dpi=100)    
    
    print("Saved to "+save_file)


def bathy_topo_map(LLD_FILE):
    
    save_file = os.path.splitext(LLD_FILE)[0]+'_bathymap.png'    
    
    # Read bathymetry/topography data
    data = np.genfromtxt(LLD_FILE, dtype=[('lat','f8'),('lon','f8'), ('z','f8')],skip_header=3)

    # Data ranges
    lon_min,lon_max = data['lon'].min(),data['lon'].max()
    lat_min,lat_max = data['lat'].min(),data['lat'].max()
    mean_lat = 0.5*(lat_min + lat_max)
    mean_lon = 0.5*(lon_min + lon_max)
    lon_range = lon_max - lon_min
    lat_range = lat_max - lat_min
    cmap = plt.get_cmap('terrain')
    interp = 'none'
    framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14)
    
    # Reshape into matrices
    Ncols = len(np.unique(data['lon']))
    Nrows = len(np.unique(data['lat']))
    
    X = data['lon'].reshape(Nrows, Ncols)
    Y = data['lat'].reshape(Nrows, Ncols)
    Z = data['z'].reshape(Nrows, Ncols)
        
    # catch any nan values
    masked_data = np.ma.masked_invalid(Z)
    cmap.set_bad('red')
    
    # Color limits
    z_min,z_max = masked_data.min(),masked_data.max()
    z_lim = max(np.abs(z_min),np.abs(z_max))
    norm = mcolor.Normalize(vmin=-z_lim, vmax=z_lim)

    # Initialize the frame and axes
    fig = plt.figure()
    
    m = Basemap(projection='cyl',llcrnrlat=lat_min, urcrnrlat=lat_max,
                llcrnrlon=lon_min, urcrnrlon=lon_max, lat_0=mean_lat, lon_0=mean_lon, resolution='h')
    m.ax = fig.add_subplot(111)
    
    m.drawmeridians(np.linspace(lon_min,lon_max,num=5.0),labels=[0,0,0,1], linewidth=0)
    m.drawparallels(np.linspace(lat_min,lat_max,num=5.0),labels=[1,0,0,0], linewidth=0)
    m.drawcoastlines(linewidth=0.5)

    # Colorbar
    divider = make_axes_locatable(m.ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
    plt.figtext(0.96, 0.7, r'elevation $[m]$', rotation='vertical', fontproperties=framelabelfont)
    cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)    
        
    # Plot the contours
    #m.contourf(X, Y, masked_data, 100, cmap=cmap, norm=norm, extend='both', zorder=1)
    m.ax.imshow(masked_data,cmap=cmap,origin='lower',norm=norm,extent=[lon_min,lon_max,lat_max,lat_min],interpolation=interp)

    plt.savefig(save_file,dpi=100)
    print("Saved to "+save_file)
    

# --------------------------------------------------------------------------------
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()    
    subparsers = parser.add_subparsers(title='mode', description='valid modes of usage', dest='mode')     
    
    # Arguments for generating bathymetry LLD files
    parser_gen = subparsers.add_parser('generate_bathy', help='generate interpolated bathymetry subset from NOAA\'s ETOPO1 dataset')
    parser_gen.add_argument('--info_file', required=True, help='json containing regional and earthquake information')  
    parser_gen.add_argument('--etopo1_file', required=False, default="~/Tsunami/ETOPO1/ETOPO1_Ice_g_gmt4.grd",
                            help="NOAA ETOPO1 combined topography and bathymetry file path")
    parser_gen.add_argument('--resolution', type=int, required=False, default=1,
                            help="Resolution interpolation multiplier for NOAA topography")
    parser_gen.add_argument('--text', action="store_true", required=False,
                            help="Store bathymetry as a text lld file")
    
    # Arguments for plotting bathymetry 
    parser_plot_bathy = subparsers.add_parser('plot_bathy', help='Plot previously-generated bathymetry LLD file')
    parser_plot_bathy.add_argument('--lld_file', required=True,
            help="Path of bathymetry file")    
    
    # Arguments for generating EQ surface displacement fields
    parser_field_eval = subparsers.add_parser('eq_field_eval', help='Generate surface displacement for a VQ fault model')
    parser_field_eval.add_argument('--info_file', required=True, help='json containing regional and earthquake information')  
    parser_field_eval.add_argument('--lld_file', required=True,
            help="Path of bathymetry file")
    parser_field_eval.add_argument('--slip_from', choices=['vq_sim', 'uniform', 'slipmap'], required=True,
                                  help="Sources of displacements")
    
    # Arguments for plotting eq displacement fields
    parser_field_plot = subparsers.add_parser('eq_field_plot', help='Plot a surface displacement field')
    parser_field_plot.add_argument('--field_file', required=True,
            help="Path of surface displacement file") 
    parser_field_plot.add_argument('--plot_horizontal', action="store_true", required=False,
            help="Whether to plot the horizontal displacements in addition to vertical") 
    
    # Arguments for plotting simulation results  
    parser_animate = subparsers.add_parser('animate', help='Make animations from simulation output')        
    parser_animate.add_argument('--type', choices=['grid', 'xsection'], required=True,
            help="Type of animation to make; birds-eye grid or cross-sectional verification against analytic solution")
    parser_animate.add_argument('--sim_file', required=True,
            help="Name of simulation file to analyze.")
    parser_animate.add_argument('--zminmax', type=float, nargs=2, required=False,
            help="Bounds for water height color bar")
    parser_animate.add_argument('--use_basemap', action="store_true", required=False,
            help="Whether to plot a basemap coastal outline over the grid animation")
    parser_animate.add_argument('--fps', type=int, required=False, default=20,
            help="Frames per second for animations")
    parser_animate.add_argument('--dpi', type=int, required=False, default=100,
            help="Bounds for water height color bar")
    
    # Arguments for verifying simualtion against observed historical runups
    parser_verify = subparsers.add_parser('verify', help='Verify simulation results against NOAA historical runup data')
    parser_verify.add_argument('--sim_file', required=True,
            help="Name of simulation file to analyze.")
    parser_verify.add_argument('--obs_file', required=True,
            help="Observed historical tsunami runup CSV file")
    parser_verify.add_argument('--ymd', nargs=3, type=int, required=True,
            help="Year, Month, and Day of historical runup information needed")
    
    
#    args = parser.parse_args(['verify', '--sim_file', 'outputs/Tohoku/TohokuSmall_x2_contRealisticX1_output_000-1600.nc', 
#                              '--obs_file', '../historical_runups/tsrunup.csv', '--ymd', '2011', '3', '11'])
    args = parser.parse_args()
    
    
    if args.mode == 'generate_bathy':
        
        with open(args.info_file, "r") as open_info_file:
            region_info = json.load(open_info_file)        
        
        # ====== PARSE ETOPO1 FILE, SAVE SUBSET =====
        ETOPO1_FILE = args.etopo1_file
        FACTOR  = args.resolution
        SAVE_NAME = os.path.join(os.path.split(args.info_file)[0], 'bathymetry', region_info['name']+'_x'+str(FACTOR)+'_lld.txt')
        
        save_dir = os.path.join(os.path.split(args.info_file)[0], 'bathymetry')
        if args.text:
            save_name = region_info['name']+'_x'+str(FACTOR)+'_lld.txt'
        else:
            save_name = region_info['name']+'_x'+str(FACTOR)+'_lld.nc'
        
        if not os.path.isdir(save_dir):
            os.makedirs(save_dir)
        
        SAVE_PATH = os.path.join(save_dir, save_name)
        
        MIN_LAT = region_info['lat_bounds'][0]
        MAX_LAT = region_info['lat_bounds'][1]
        MIN_LON = region_info['lon_bounds'][0]
        MAX_LON = region_info['lon_bounds'][1]

        # --- write grid ------
        # TODO: transition from txt files to netCDF files containing bathymetry data
        lats, lons, bathy = read_ETOPO1.grab_ETOPO1_subset_interpolated(ETOPO1_FILE, min_lat=MIN_LAT, max_lat=MAX_LAT, min_lon=MIN_LON, max_lon=MAX_LON, factor=FACTOR)
        if args.text:
            read_ETOPO1.write_grid(SAVE_PATH, lats, lons, bathy)
        else:
            read_ETOPO1.write_grid_netCDF(SAVE_PATH, lats, lons, bathy)
        
    
    if args.mode == 'plot_bathy':
        bathy_topo_map(args.lld_file)
    
    
    if args.mode == 'eq_field_eval':
        
        with open(args.info_file, "r") as open_info_file:
            region_info = json.load(open_info_file)          
        
        VQ_DIR = "~/VirtQuake/"
        
        system_command = "python "+VQ_DIR+"vq/PyVQ/pyvq/pyvq.py --field_eval --netCDF --horizontal --model_file {} --lld_file {}".format(region_info['model_file'], args.lld_file)
        
        if args.slip_from == 'vq_sim':
            if not region_info['event_file'] or not region_info['event_id']:
                raise BaseException("Must specify both an event file and an event ID in info file")
                
            system_command += " --event_file {} --event_id {}".format(region_info['event_file'], region_info['event_id'])
            
        elif args.slip_from == 'uniform':
            system_command += " --uniform_slip 10"      
        
        elif args.slip_from == 'slipmap':
            if not region_info['slip_map']:
                raise BaseException("Must specify slipmap file")
            system_command += " --slipmap_file {}".format(region_info['slip_map'])   
           
        os.system(system_command)
       
       
    if args.mode == 'eq_field_plot':
        
        if not args.plot_horizontal:
            plot_eq_displacements(args.field_file)
        else:
            if not os.path.splitext(args.field_file)[1] == 'xyuen':
                raise BaseException("Must have .xyuen file format for horizontal field plotting")
            plot_eq_disps_horiz(args.field_file)
    
    
    if args.mode in ['animate', 'verify']:
        
        this_sim = simAnalyzer(args.sim_file)

        if args.mode == 'verify':
            # Load historical inundation data and compare it to simulated results
            this_sim.load_historical_runup_CSV(args.obs_file, args.ymd[0], args.ymd[1], args.ymd[2])
            this_sim.compare_sim_and_obs_runup()
        
        if args.mode == 'animate':
            if args.type == 'grid':
                #zminmax = (-1,1)#(-sim_data['z'].std(), sim_data['z'].std())
                # Makes animation
                this_sim.make_grid_animation(args.fps, args.dpi, zminmax = args.zminmax, doBasemap = args.use_basemap)
            
            if args.type == 'xsection':
                this_sim.make_crosssection_animation(args.fps, args.dpi)

                
        
        
        
        
        
