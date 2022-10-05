#%%
import numpy as np
from astrodendro import Dendrogram, ppv_catalog
from astropy import units as u
from astropy.table import Table, Column
from astropy.io import fits
from scimes import SpectralCloudstering
from astropy.wcs import WCS
import os
import random
import shutil
import math
# %%
def cluster(file, sigma, voxels):
    '''
    Identify molecular clouds using Scimes from a given 3D data cube. The output file 'clusters.fits' recording the masks of the identified clouds, i.e., the voxels belonging to a cloud are marked with the same dendrogram ID.
    '''
    dendrofile = 'dendro.fits'
    clusterfile = 'clusters.fits'
    catalogfile = 'dendro_catalog.cat'
    selected_catalogfile = 'clusters.cat'
    sigma = sigma
    ppb = voxels
    min_delt = 3.0*sigma
    min_v = 2.0*sigma
    path = 'ppb'+str(ppb)
    cubefile = file

    print('Make dendrogram from the full cube')
    
    data = fits.getdata(cubefile)
    hd = fits.getheader(cubefile)
    d = Dendrogram.compute(data, min_value=min_v, min_delta=min_delt, min_npix=ppb, verbose = 1)
    d.save_to(dendrofile)


    print("Generate a catalog of dendrogram structures")
    metadata = {}
    metadata['data_unit'] = u.Jy
    metadata['spatial_scale'] = 30. * u.arcsec
    metadata['velocity_scale'] = hd['CDELT3'] * u.meter / u.second
    w = WCS(naxis=3)
    w.wcs.crpix = [hd['CRPIX1'], hd['CRPIX2'], hd['CRPIX3']]
    w.wcs.cdelt = [hd['CDELT1'], hd['CDELT2'], hd['CDELT3']]
    w.wcs.crval = [hd['CRVAL1'], hd['CRVAL2'], hd['CRVAL3']]
    w.wcs.ctype = [hd['CTYPE1'], hd['CTYPE2'], hd['CTYPE3']]
    w.wcs.set_pv([(hd['NAXIS1'], hd['NAXIS2'], hd['NAXIS3'])])
    metadata['wcs'] = w


    cat = ppv_catalog(d, metadata)   
    cat.write(catalogfile, format = 'ascii', overwrite=True)

    print("Running SCIMES")
    dclust = SpectralCloudstering(d, cat, hd, criteria = ['volume'])
    # dclust.asgncube(header=hd)
    clust_hdu = dclust.clusters_asgn
    # clust_hdu = dclust.asgn 
    clust_hdu.writeto(clusterfile, overwrite=True)

    clusts = dclust.clusters
    cat_clusters = cat[clusts]
    cat_selected_clusters = cat_clusters   
    cat_selected_clusters.write(selected_catalogfile, format='ascii', overwrite=True)


    if os.path.exists(path): shutil.rmtree(path)
    os.makedirs(path)
    record = open(path+'/parameters.txt', 'w')
    record.write('MIN_NPIX: '+str(ppb)+'\n')
    record.write('MIN_DELT: '+str(min_delt/sigma)+' sigma\n')
    record.write('MIN_VAL: '+str(min_v/sigma)+' sigma\n')
    record.write('Nclusts: '+str(len(clusts))+'\n')
    record.write('\n')
    record.close()
    shutil.move(dendrofile, path)
    shutil.move(clusterfile, path)
    shutil.move(catalogfile, path)
    shutil.move(selected_catalogfile, path)


if __name__ == '__main__':
    cubefile = 'G105L_resample_clip.fits'
    sigma = 0.17
    voxels = 50
    cluster(cubefile, sigma, voxels)