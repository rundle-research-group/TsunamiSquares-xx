#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.colors as mcolor
import matplotlib.animation as manimation
import matplotlib.colorbar as mcolorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap
import matplotlib.font_manager as mfont
# -------
#import quakelib
from os import system
import read_ETOPO1

"""
Tsunami Squares output files are in the following columns:

TIME (secs) |  LON (dec. degrees)  |  LAT (dec. deg.)  |  WATER HEIGHT (meters)  |  BATHY. DEPTH (meters)
  
In the code below: X = LON, Y = LAT, Z = WATER HEIGHT, ALT = BATHYMETRY DEPTH

Water height and bathymetry depth are given in units of meters relative to global mean sea level.
If a square sits on a beach at 5m above sea level and has 4m of water on it, then Z=9m and ALT=5m.

"""

VQ_DIR = "~/VirtQuake/"

# --------------------------------------------------------------------------------
def make_animation(sim_data, FPS, DPI, T_MIN, T_MAX, T_STEP, N_STEP):
    # Get ranges
    print "min, max, av, std: ", sim_data['z'].min(), sim_data['z'].max(), sim_data['z'].mean(), sim_data['z'].std()
    lon_min,lon_max = sim_data['lon'].min(),sim_data['lon'].max()
    lat_min,lat_max = sim_data['lat'].min(),sim_data['lat'].max()
    z_min,z_max = sim_data['z'].min(), sim_data['z'].max()#-sim_data['z'].std(), sim_data['z'].std()#-.1, .1#-2, 2#
    cmap = plt.get_cmap('Blues_r')
    norm = mcolor.Normalize(vmin=z_min, vmax=z_max)
    interp = 'none'
    
    # Split the data up into arrays for each time step
    split_data = np.split(sim_data, np.unique(sim_data['time']).shape[0])

    # Initialize movie writing stuff
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='TsunamiSquares', artist='Matplotlib',
            comment='Animation')
    writer = FFMpegWriter(fps=FPS, metadata=metadata, bitrate=1000)

    # Initialize the frame and axes
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlim(lon_min, lon_max)
    plt.ylim(lat_min, lat_max)
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    ax.get_yaxis().get_major_formatter().set_useOffset(False)
    
    # I don't know why, but the y-axis is backwards
    ax.invert_yaxis()
    
    divider = make_axes_locatable(ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
    cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)
    # Increment the time from T_MIN
    TIME = T_MIN

    first_step = sim_data[ sim_data['time'] == T_MIN ]
    Ncols = len(np.unique(first_step['lon']))

    surface = None
    with writer.saving(fig, save_file, DPI):
        for index in range(int(N_STEP)):
            # Get the subset of data corresponding to current time
            this_step = split_data[index]
            time = this_step['time'][0]
            
            print "step: "+str(index)+"  time: "+str(time)+" num_points: "+str(len(this_step))
            assert len(this_step) > 0

            X = this_step['lon'].reshape(-1, Ncols)
            Y = this_step['lat'].reshape(-1, Ncols)
            Z = this_step['z'].reshape(-1, Ncols)
            ALT = this_step['alt'].reshape(-1, Ncols)

            # Plot the surface for this time step
            if surface is None:
                surface = ax.imshow(Z,cmap=cmap,origin='upper',norm=norm,extent=[lon_min,lon_max,lat_max,lat_min],interpolation=interp)
            else:
                surface.set_data(Z)
                
            # Text box with the time
            plt.figtext(0.02, 0.5, 'Time: {:02d}:{:02d}'.format(int(time)/60, int(time)%60), bbox={'facecolor':'yellow', 'pad':5})
            
                
            writer.grab_frame()
            
            TIME +=T_STEP


# --------------------------------------------------------------------------------
def make_map_animation(sim_data, FPS, DPI, T_MIN, T_MAX, T_STEP, N_STEP, save_file):
    # Get ranges
    lon_min,lon_max = sim_data['lon'].min(),sim_data['lon'].max()
    lat_min,lat_max = sim_data['lat'].min(),sim_data['lat'].max()
    mean_lat = 0.5*(lat_min + lat_max)
    mean_lon = 0.5*(lon_min + lon_max)
    lon_range = lon_max - lon_min
    lat_range = lat_max - lat_min
    z_min,z_max = -sim_data['z'].std(),sim_data['z'].std()#sim_data['z'].min(),sim_data['z'].max()
    cmap = plt.get_cmap('Blues_r')
    norm = mcolor.Normalize(vmin=z_min, vmax=-z_min)
    interp = 'none'
    landcolor = '#FFFFCC'
    framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14)
    
    # Split the data up into arrays for each time step
    split_data = np.split(sim_data, np.unique(sim_data['time']).shape[0])
    
    # Initialize movie writing stuff
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='TsunamiSquares', artist='Matplotlib',
                    comment='Basemap Animation')
    writer = FFMpegWriter(fps=FPS, metadata=metadata)
    
    # Initialize the frame and axes
    fig = plt.figure()

    m = Basemap(projection='cyl',llcrnrlat=lat_min, urcrnrlat=lat_max,
            llcrnrlon=lon_min, urcrnrlon=lon_max, lat_0=mean_lat, lon_0=mean_lon, resolution='h')
    m.ax = fig.add_subplot(111)

    m.drawmeridians(np.linspace(lon_min,lon_max,num=5.0),labels=[0,0,0,1], linewidth=0)
    m.drawparallels(np.linspace(lat_min,lat_max,num=5.0),labels=[1,0,0,0], linewidth=0)
    m.drawcoastlines(linewidth=0.5)
    #m.drawcountries()
    #m.drawstates()
    #m.fillcontinents(color=landcolor)
    #m.shadedrelief()
    
    # Colorbar
    divider = make_axes_locatable(m.ax)
    cbar_ax = divider.append_axes("right", size="5%",pad=0.05)
    plt.figtext(0.95, 0.7, r'water altitude $[m]$', rotation='vertical', fontproperties=framelabelfont)
    cb = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norm)
    
    # Increment the time from T_MIN
    TIME = T_MIN
    
    first_step = sim_data[ sim_data['time'] == T_MIN ]
    Ncols = len(np.unique(first_step['lon']))
    
    surface = None
    with writer.saving(fig, save_file, DPI):
        for index in range(int(N_STEP)):
            # Get the subset of data corresponding to current time
            this_step = split_data[index]
            time = this_step['time'][0]
                    
            print "step: "+str(index)+"  time: "+str(time)+" num_points: "+str(len(this_step))
            assert len(this_step) > 0
                
            X = this_step['lon'].reshape(-1, Ncols)
            Y = this_step['lat'].reshape(-1, Ncols)
            Z = this_step['z'].reshape(-1, Ncols)
            ALT = this_step['alt'].reshape(-1, Ncols)
            
            # Masked array via conditional, don't color the land unless it has water on it
            masked_data = np.ma.masked_where(np.logical_and(np.array(Z == 0.0),np.array(ALT >= 0.0)), Z)
            
            # Set masked pixels to the land color
            cmap.set_bad(landcolor, 1.0)  # set alpha=0.0 for transparent
            
            # Plot the surface for this time step
            if surface is None:
                surface = m.ax.imshow(masked_data,cmap=cmap,origin='lower',norm=norm,extent=[lon_min,lon_max,lat_max,lat_min],interpolation=interp)
            else:
                surface.set_data(masked_data)
                
            # Text box with the time
            plt.figtext(0.129, 0.82, 'Time: {:02d}:{:02d}'.format(int(time)/60, int(time)%60), bbox={'facecolor':'yellow', 'pad':5})
                
            writer.grab_frame()
        
            TIME +=T_STEP


# =============================================================
def plot_eq_displacements(LLD_FILE, LEVELS, save_file):
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

# =============================================================
def plot_eq_disps_horiz(disp_file, save_file):
    # Read displacement data
    disp_data = np.genfromtxt(disp_file, dtype=[('lon','f8'), ('lat','f8'), ('z','f8'), ('eU','f8'), ('nV','f8')])

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

    plt.savefig(save_file+'_z.png',dpi=100)
    
    
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
    cbz = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=norme)
    
    # Masked array via conditional, don't color the land unless it has water on it
    zero_below = int(len(LEVELSe)/2)-1
    zero_above = zero_below+1
    masked_data = np.ma.masked_where(np.logical_and(np.array(eU <= LEVELSe[zero_above]),np.array(eU >= LEVELSe[zero_below])), eU)
    
    # Set masked pixels to the land color
    cmap.set_bad(landcolor, 0.0)  # set alpha=0.0 for transparent
    
    # Plot the contours
    m.contourf(X, Y, masked_data, LEVELSe, cmap=cmap, norm=norme, extend='both', zorder=1)

    plt.savefig(save_file+'_e.png',dpi=100)


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
    cbz = mcolorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=normn)
    
    # Masked array via conditional, don't color the land unless it has water on it
    zero_below = int(len(LEVELSn)/2)-1
    zero_above = zero_below+1
    masked_data = np.ma.masked_where(np.logical_and(np.array(nV <= LEVELSn[zero_above]),np.array(nV >= LEVELSn[zero_below])), nV)
    
    # Set masked pixels to the land color
    cmap.set_bad(landcolor, 0.0)  # set alpha=0.0 for transparent
    
    # Plot the contours
    m.contourf(X, Y, masked_data, LEVELSn, cmap=cmap, norm=normn, extend='both', zorder=1)

    plt.savefig(save_file+'_n.png',dpi=100)    
    
    
    print("Saved to "+save_file)


# =============================================================
def bathy_topo_map(LLD_FILE, save_file):
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
    
#    MODE = "generate"
    MODE = "animate"
#    MODE = "eq_field_plot"
#    MODE = "eq_field_plot_horiz"
#    MODE = "plot_bathy"
#    MODE = "gen_bathyfile_interp"
#    MODE = "eq_field_eval"
    
    if MODE == "generate": #read bethymetry file
        # ====== PARSE ETOPO1 FILE, SAVE SUBSET, EVALUATE EVENT FIELD AT THE LAT/LON, SAVE =====
        ETOPO1_FILE = "ETOPO1/ETOPO1_Ice_g_gmt4.grd"
        # =================================
#        SAVE_NAME = "bathymetry/Channel_Islands_largest_subset_lld.txt"
#        MODEL     = "~/VirtQuake/UCERF3_reindexed/Model/UCERF3_VQmeshed_from_EQSIM_ReIndexed_AseismicCut_0-11_taper_drops0.9.h5"
#        EVENTS    = "~/VirtQuake/UCERF3_reindexed/VQ_runs/events_UCERF3_ReIndexed_AseismicCut_0-11_taper_drops0-9_50kyr_dyn0-2_greenLimits.h5"
#        EVID      = 39951
#        # Larger subset
#        MIN_LAT = 33.75
#        MAX_LAT = 34.3
#        MIN_LON = -120.2
#        MAX_LON = -119.2
        # =================================
        SAVE_NAME = "bathymetry/Tohoku_lld.txt"
        MODEL     = "~/VirtQuake/Tohoku/simple_Tohoku_50000m_drops0.txt"
        EVENTS    = "~/VirtQuake/Tohoku/events_Tohoku_100kyr_drops0_dyn0-2_greenLimits.h5"
        EVID      = 1196
        MIN_LAT = 33
        MAX_LAT = 44
        MIN_LON = 138.4
        MAX_LON = 145
        # =================================
        # --- write grid ------
        lats,lons,bathy = read_ETOPO1.grab_ETOPO1_subset(ETOPO1_FILE,min_lat=MIN_LAT,max_lat=MAX_LAT,min_lon=MIN_LON,max_lon=MAX_LON)
        read_ETOPO1.write_grid(SAVE_NAME,lats,lons,bathy)
        # ---- compute field and write it ------
        system("python "+VQ_DIR+"vq/PyVQ/pyvq/pyvq.py --field_eval --netCDF  --horizontal --model_file {} --event_file {} --event_id {} --lld_file {} ".format(MODEL, EVENTS, EVID, SAVE_NAME))
    
    if MODE == "animate":
        sim_file = "tsunami_output.txt"
        save_file = sim_file.split(".")[0]+".mp4"
        sim_data = np.genfromtxt(sim_file, dtype=[('time','f8'),('lat','f8'),('lon','f8'), ('z','f8'), ('alt','f8')])
        FPS = 10 #FRAMES PER SECOND
        DPI = 100
        T_MAX,T_MIN = sim_data['time'].max(),sim_data['time'].min()
        T_STEP = np.unique(sim_data['time'])[1] - np.unique(sim_data['time'])[0]
        assert T_STEP > 0
        N_STEP = float(T_MAX-T_MIN)/T_STEP
        # Makes animation on a Basemap plot, currently misbehaving
#        make_map_animation(sim_data, FPS, DPI, T_MIN, T_MAX, T_STEP, N_STEP, save_file)
        # Makes animation without any background Basemap
        make_animation(sim_data, FPS, DPI, T_MIN, T_MAX, T_STEP, N_STEP)

    if MODE == "eq_field_plot":
        Levels = [-.3, -.2, -.1, -.05, -.008, .008, .05, .1, .2, .3]
        plot_eq_displacements("bathymetry/Channel_Islands_largest_subset_lld_dispField_event39951.txt",Levels, "outputs/disp_map.png")
    
    if MODE == "eq_field_plot_horiz":
        plot_eq_disps_horiz("bathymetry/Tohoku_lld_dispField_realisticSlips.xyuen", "outputs/realisticSlips_disp_map")
        
    if MODE == "plot_bathy":
        #Levels = [-3, -.2, -.1, -.05, -.008, .008, .05, .1, .2, .3]
        #bathy_topo_map("local/Channel_Islands.txt",Levels, "bathy_map.png")
        bathy_topo_map("local/Channel_Islands_interp_larger.txt", "bathy_map_interp_larger_imshow.png")
        
    if MODE == "gen_bathyfile_interp":
        # ====== PARSE ETOPO1 FILE, SAVE SUBSET, EVALUATE EVENT FIELD AT THE LAT/LON, SAVE =====
        ETOPO1_FILE = "ETOPO1/ETOPO1_Ice_g_gmt4.grd"
        # =================================
#        SAVE_NAME = "bathymetry/Channel_Islands_largest_subset_lld.txt"
#        MODEL     = "~/VirtQuake/UCERF3_reindexed/Model/UCERF3_VQmeshed_from_EQSIM_ReIndexed_AseismicCut_0-11_taper_drops0.9.h5"
#        EVENTS    = "~/VirtQuake/UCERF3_reindexed/VQ_runs/events_UCERF3_ReIndexed_AseismicCut_0-11_taper_drops0-9_50kyr_dyn0-2_greenLimits.h5"
#        EVID      = 39951
#        # Larger subset
#        MIN_LAT = 33.75
#        MAX_LAT = 34.3
#        MIN_LON = -120.2
#        MAX_LON = -119.2
        # =================================
        SAVE_NAME = "bathymetry/Tohoku_lld.txt"
        MODEL     = "~/VirtQuake/Tohoku/simple_Tohoku_50000m_drops0.txt"
        EVENTS    = "~/VirtQuake/Tohoku/events_Tohoku_100kyr_drops0_dyn0-2_greenLimits.h5"
        EVID      = 1196
        MIN_LAT = 33
        MAX_LAT = 44
        MIN_LON = 138.4
        MAX_LON = 145
        # =================================
        FACTOR  = 3
        # --- write grid ------
        lats, lons, bathy = read_ETOPO1.grab_ETOPO1_subset_interpolated(ETOPO1_FILE,min_lat=MIN_LAT,max_lat=MAX_LAT,min_lon=MIN_LON,max_lon=MAX_LON, factor=FACTOR)
        read_ETOPO1.write_grid(SAVE_NAME,lats,lons,bathy)
        # ---- compute field and write it ------
        system("python "+VQ_DIR+"vq/PyVQ/pyvq/pyvq.py --field_eval --netCDF  --horizontal --model_file {} --event_file {} --event_id {} --lld_file {} ".format(MODEL, EVENTS, EVID, SAVE_NAME))


    if MODE == "eq_field_eval":
#        LLD_NAME = "bathymetry/Channel_Islands_largest_subset_lld.txt"
#        MODEL     = "~/VirtQuake/UCERF3_reindexed/Model/UCERF3_VQmeshed_from_EQSIM_ReIndexed_AseismicCut_0-11_taper_drops0.9.h5"
#        EVENTS    = "~/VirtQuake/UCERF3_reindexed/VQ_runs/events_UCERF3_ReIndexed_AseismicCut_0-11_taper_drops0-9_50kyr_dyn0-2_greenLimits.h5"
#        EVID      = 39951
        # =================================
        SAVE_NAME = "bathymetry/Tohoku_lld.txt"
        MODEL     = "~/VirtQuake/Tohoku/simple_Tohoku_50000m_drops0.txt"
        EVENTS    = "~/VirtQuake/Tohoku/events_Tohoku_100kyr_drops0_dyn0-2_greenLimits.h5"
        EVID      = 1196
        # ---- compute field and write it ------
#        system("python "+VQ_DIR+"vq/PyVQ/pyvq/pyvq.py --field_eval --horizontal --model_file {} --event_file {} --event_id {} --lld_file {} ".format(MODEL, EVENTS, EVID, SAVE_NAME))
        system("python "+VQ_DIR+"vq/PyVQ/pyvq/pyvq.py --field_eval --netCDF --horizontal --model_file {} --event_id {} --lld_file {} ".format(MODEL, EVID, SAVE_NAME))



