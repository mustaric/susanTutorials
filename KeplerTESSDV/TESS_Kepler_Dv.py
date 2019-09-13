#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 16:26:55 2019

@author: smullally
"""
from astroquery.mast import Observations
from astroquery.mast import Catalogs
from astropy.io import fits
from astropy import table
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
#%%
#Goals: Obtain all time series data from Kepler and TESS on a given target 
#using astroquery.

star_name = "L98-59"

#Use Astroquery to find all dvt products
#Start by getting all time series observations at exactly the coordinates of our star.

observations = Observations.query_object(star_name,radius="0 deg")
obs_wanted = observations['dataproduct_type'] == 'timeseries'
print(observations[obs_wanted]['obs_collection','project','obs_id'])

#Get all Products and determine which have DV files.
#Then plot the time series for each that has DVT files. 


data_products = Observations.get_product_list(observations[obs_wanted])

files_wanted = (data_products["productSubGroupDescription"] == "DVT") | \
                (data_products["productSubGroupDescription"] == "DVM") | \
                (data_products["productSubGroupDescription"] == "DVS") | \
                (data_products["productSubGroupDescription"] == "DVR")
    
    
#files_wanted = list( 
#        map( lambda x: x[-8:] == 'dvt.fits', data_products['productFilename'] ))

print(data_products["productFilename"][files_wanted])
manifest = Observations.download_products(data_products[files_wanted])


#%%


def parse_manifest(manifest):
    """
    Parse manifest and add back columns that are useful for TESS DV exploration.
    """
    results = deepcopy(manifest)
    filenames = []
    sector_range = []
    exts = []
    for i,f in enumerate(manifest['Local Path']):
        file_parts = np.array(np.unique(f.split(sep='-')))
        sectors = list( map ( lambda x: x[0:2] == 's0', file_parts))
        s1 = file_parts[sectors][0]
        try:
            s2 = file_parts[sectors][1]
        except:
            s2 = s1
        sector_range.append("%s-%s" % (s1,s2))
        path_parts = np.array(f.split(sep='/'))
        filenames.append(path_parts[-1])
        exts.append(path_parts[-1][-8:])

    results.add_column(table.Column(name="filename",data=filenames))
    results.add_column(table.Column(name="sectors",data=sector_range))
    results.add_column(table.Column(name="fileType", data=exts))
    results.add_column(table.Column(name="index", data=np.arange(0,len(manifest))))
    
    return results
#
#Examine the Manifest that is returned
#

results = parse_manifest(manifest)
print(results['index','sectors','fileType','filename'][results['sectors'] == "s0001-s0013"])

     
#%%
#Plot the light curve from the DVT file for each TCE in the file
want = (results['sectors'] == "s0001-s0013") & (results['fileType'] == "dvt.fits")
dvt_filename = manifest[want]['Local Path'][0]

fits.info(dvt_filename)

#%%
#Notice that there are 3 TCEs found
#Plot the entire detrended light curve sent into the plant hunting algorithm.

data = fits.getdata(dvt_filename,1)
time = data['TIME']
relflux = data['LC_DETREND']

plt.figure(figsize=(10,4))
plt.plot (time, relflux, 'b.')
plt.ylim(1.2* np.nanpercentile(relflux, .5) , 1.2 * np.nanpercentile(relflux,99.5))
plt.title('Data Validation Detrended Light Curve for %s' % (star_name))

#%%
#
#Plot the folded light curves. First write a function 

def plot_folded(phase, data, model, ext, period):
    isort = phase.argsort()
    
    plt.plot(phase[isort], data[isort], '.', ms=.5)
    plt.plot(phase[isort], model[isort], '-', lw=1, label="TCE %i" % ext)
    plt.xlabel('Phase (Period = %5.2f days)' % period)
    plt.ylim(1.5* np.nanpercentile(data, .5) , 1.4 * np.nanpercentile(data,99.5))
    plt.legend(loc="lower right")

#%%
plt.figure(figsize=(10,12))

nTCEs = fits.getheader(dvt_filename)['NEXTEND'] - 2

for ext in range(1,nTCEs+1):
    data = fits.getdata(dvt_filename, ext)
    head = fits.getheader(dvt_filename, ext)
    period = head['TPERIOD']
    phase = data['PHASE']
    flux = data['LC_INIT']
    model = data['MODEL_INIT']
    plt.subplot(3,1,ext)
    plot_folded(phase,flux,model, ext, period)
    

    
    
