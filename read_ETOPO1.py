#!/usr/bin/env python

import sys
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

def get_subset_indices(min_lat, max_lat, min_lon, max_lon, lats, lons):
    # These are the indices of the closest lat/lon values to
    # (min_lat, max_lat, min_lon, max_lon)
    indices = []
    indices.append(int((np.abs(np.array(lats)-min_lat)).argmin()))
    indices.append(int((np.abs(np.array(lats)-max_lat)).argmin()))
    indices.append(int((np.abs(np.array(lons)-min_lon)).argmin()))
    indices.append(int((np.abs(np.array(lons)-max_lon)).argmin())) 
    # return [min_lat_index, max_lat_index, min_lon_index, max_lon_index]
    return indices


def grab_ETOPO1_subset(file_name, min_lat, max_lat, min_lon, max_lon):
    
    ETOPO1 = Dataset(file_name, 'r')
    lons = ETOPO1.variables["x"][:]
    lats = ETOPO1.variables["y"][:]
    
    # Grab indices for max/min lat/lon bounds
    minLat, maxLat, minLon, maxLon = get_subset_indices(min_lat, max_lat, min_lon, max_lon, lats, lons)
    bathy = ETOPO1.variables["z"][minLat:maxLat,minLon:maxLon]
    lons,lats = np.meshgrid(lons[minLon:maxLon],lats[minLat:maxLat])
    print("== Selected {} points ({}x{}) from {}".format(bathy.size,bathy.shape[1],bathy.shape[0],file_name))   
    print("---- Lats: {} to {},   Lons: {} to {}".format(min_lat, max_lat, min_lon, max_lon))
            
    return lats,lons,bathy
    
    
def write_grid(out_file_name, lats, lons, bathy):
    outfile = open(out_file_name, 'w')

    dlon = abs(np.unique(lons)[1] - np.unique(lons)[0])
    dlat = abs(np.unique(lats)[1] - np.unique(lats)[0])
    
    outfile.write("# N_lats\n")
    outfile.write("# N_lons\n")
    outfile.write("# dlat\n")
    outfile.write("# dlon\n")
    outfile.write("{} {} {} {}\n".format(lons.shape[0], lons.shape[1], dlat, dlon))
    outfile.write("##################################\n")
    # Write vertices from top left (Northwest) to bottom right (Southeast)
    for i in list(reversed(range(lons.shape[0]))):
        for j in range(lons.shape[1]):
            outfile.write("{}\t{}\t{}\n".format(lats[i][j],lons[i][j],bathy[i][j]))
    outfile.close()
    print("output written to {}".format(out_file_name))
    
    
    
def grab_ETOPO1_subset_interpolated(file_name, min_lat, max_lat, min_lon, max_lon, factor=3):
    debug = False

    ETOPO1 = Dataset(file_name, 'r')
    
    #adjusting lons to be able to cross the dateline
    lons = ETOPO1.variables["x"][:]
    lons2 = np.copy(lons)
    if min_lon>0 and max_lon<0:
        for i in range(len(lons)):
            if lons[i]>0:
                lons[i]=lons[i]-180
            elif lons[i]<0:
                lons[i]=lons[i]+180
            elif lons[i]==180 or lons[i]==-180:
                lons[i]=0
            else:
                lons[i]=180

    lats = ETOPO1.variables["y"][:]

    # Extend the bounds to add a buffer, ensures that interpolation has enough data on the boundaries
    # TODO: better algorithm here, the min/max position of 0.98 vs 1.02 changes in different regions
    # TODO: need to handle crossing the int date line
    min_lat_buff = min_lat-3
    max_lat_buff = max_lat+3
    min_lon_buff = min_lon-3
    max_lon_buff = max_lon+3
    
    # Grab indices for larger buffered max/min lat/lon bounds
    minLat_buff, maxLat_buff, minLon_buff, maxLon_buff = get_subset_indices(min_lat_buff, max_lat_buff, min_lon_buff, max_lon_buff, lats, lons2)
    minLat_buff2, maxLat_buff2, minLon_buff2, maxLon_buff2 = get_subset_indices(min_lat_buff, max_lat_buff, min_lon_buff, max_lon_buff, lats, lons)
    lons_buff,lats_buff = np.meshgrid(lons2[minLon_buff2:maxLon_buff2],lats[minLat_buff2:maxLat_buff2])
    if min_lon>0 and max_lon<0:
        left_bathy = ETOPO1.variables["z"][minLat_buff:maxLat_buff,minLon_buff:-1]
        right_bathy = ETOPO1.variables["z"][minLat_buff:maxLat_buff,:maxLon_buff]
        bathy_buff = np.zeros((len(left_bathy),len(left_bathy[0])+len(right_bathy[0])))
        for i in range(len(left_bathy)):
            bathy_buff[i] = np.append(left_bathy[i],right_bathy[i])
    else:
        bathy_buff = ETOPO1.variables["z"][minLat_buff:maxLat_buff,minLon_buff:maxLon_buff]

    #Grab indices for requested lat/lon bounds
    minLat, maxLat, minLon, maxLon = get_subset_indices(min_lat, max_lat, min_lon, max_lon, lats, lons)
    #lons,lats = np.meshgrid(lons[minLon:maxLon],lats[minLat:maxLat])
    lons = lons[minLon:maxLon]
    lats = lats[minLat:maxLat]
    bathy = ETOPO1.variables["z"][minLat:maxLat,minLon:maxLon]
    
    print("== Selected {} points ({}x{}) from {}".format(bathy.size, bathy.shape[1], bathy.shape[0], file_name))   
    print("---- Lats: {} to {},   Lons: {} to {}".format(min_lat, max_lat, min_lon, max_lon))
    
    # Flatten everything into 1D arrays
    bathy_flat = bathy_buff.flatten()
    lons_buff_flat = lons_buff.flatten()
    lats_buff_flat = lats_buff.flatten()
    print(len(lons_buff_flat),len(lats_buff_flat))
    points = list(zip(lons_buff_flat,lats_buff_flat))
    
    # Create the grid for interpolation
    grid_lon_vals = np.linspace(min_lon, max_lon, num=int(len(lons)*factor))
    grid_lat_vals = np.linspace(min_lat, max_lat, num=int(len(lats)*factor))
    if min_lon>0 and max_lon<0:
        grid_lon_vals = np.linspace(min_lon-180, max_lon+180, num=int(len(lons)*factor))
        grid_lon_vals2 = np.copy(grid_lon_vals) #actual lons array (not transformed)
        for i in range(grid_lon_vals2.size):
            if grid_lon_vals2[i]>0:
                grid_lon_vals2[i]-=180
            else:
                grid_lon_vals2[i]+=180
        grid_lons2, grid_lats2 = np.meshgrid(grid_lon_vals2, grid_lat_vals)
    grid_lons, grid_lats = np.meshgrid(grid_lon_vals, grid_lat_vals)
    print("grid_lons shape = ", grid_lons.shape)

    print("points length = ", len(points))
    print(len(bathy_flat))
    print(grid_lons)
    print(grid_lats)
    
    if debug:
        print("buffered")
        print(minLat_buff, maxLat_buff, minLon_buff, maxLon_buff)
        print(min_lat_buff, max_lat_buff, min_lon_buff, max_lon_buff)
        print(lons_buff)
        print(" ")
        print("original")
        print(minLat, maxLat, minLon, maxLon)
        print(min_lat, max_lat, min_lon, max_lon)

    # Interpolate
    if factor>=1:
        bathy_interp = griddata(points, bathy_flat, (grid_lons, grid_lats), method='cubic')
    elif factor<1 and factor>=0.1:
        bathy_interp = griddata(points, bathy_flat, (grid_lons, grid_lats), method='nearest')
        #bathy_interp = np.zeros((int((len(bathy_buff)-360)*factor),int((len(bathy_buff[0])-360)*factor)))
        #print(len(bathy_interp),len(bathy_interp[0]))
    else:
        bathycut = bathy_buff[(len(bathy_buff)-len(bathy))/2:(len(bathy_buff)-len(bathy))/2+len(bathy),(len(bathy_buff[0])-len(bathy[0]))/2:(len(bathy_buff[0])-len(bathy[0]))/2+len(bathy[0])]
        print(len(bathycut))
        print(len(bathycut[0]))
        bathy_interp = np.zeros((len(grid_lons),len(grid_lons[0])))
        print(len(bathycut)*factor)
        print(len(bathycut[0])*factor)
        invfac = int(1.0/factor)
        invfac2 = int(1.0/factor)
        icount = -1
        jcount = -1
        for i in range(len(bathycut)):
            if invfac==i:
                invfac = int((icount+3)*1.0/factor)
                icount+=1
                jcount=-1
                for j in range(len(bathycut[0])):
                    if invfac2==j:
                        invfac2 = int((jcount+3)*1.0/factor)
                        jcount+=1
                        bathy_interp[icount,jcount] = bathycut[i,j]
                    if j==len(bathycut[0])-1:
                        invfac2 = 0
        
    print(bathy_interp)
    print("SQUARES: "+str(len(grid_lons)*len(grid_lons[0])))
    
    print("** Interpolated {} points ({}x{})".format(bathy_interp.size,bathy_interp.shape[1],bathy_interp.shape[0]))
    return grid_lats, grid_lons, bathy_interp
    return 0,0,0

    
def write_grid_netCDF(out_file_name, lats, lons, bathy):
    
    out_dataset = Dataset(out_file_name, 'w', format='NETCDF4')  
    
    lons = np.unique(lons)
    lats = np.unique(lats)
    
    out_dataset.createDimension('latitude', np.size(lats))
    out_dataset.createDimension('longitude', np.size(lons))
    
    lats_data   = out_dataset.createVariable('latitude',  'f4', ('latitude',))
    lons_data   = out_dataset.createVariable('longitude', 'f4', ('longitude',))
    altitude_data = out_dataset.createVariable('altitude',    'f4', ('latitude','longitude'))
    
    lats_data[:]       = lats
    lons_data[:]       = lons
    altitude_data[:,:] = bathy
    
    out_dataset.close()
    print("output written to {}".format(out_file_name))
    

    
    
    
    
