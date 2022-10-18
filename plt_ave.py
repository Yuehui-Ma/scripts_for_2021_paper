import os
from astropy.io import fits
import numpy as np
from matplotlib import rcParams
from matplotlib import pyplot as plt
rcParams["savefig.dpi"] = 300
rcParams["figure.dpi"] = 300
rcParams["font.size"] = 15
rcParams['ytick.direction'] = 'in'
rcParams['xtick.direction'] = 'in'

def plt_ave(path, region):
    file1 = path + region + '/' + region+'U.fits'
    file2 = path + region + '/' + region+'L.fits'
    figname = region+'_ave.pdf'
    data1 = fits.getdata(file1)
    hd1 = fits.getheader(file1)
    ch1 = np.linspace(0, hd1['NAXIS3']-1, num=hd1['NAXIS3']) 
    v1 = ((ch1 - hd1['CRPIX3'] + 1) * hd1['CDELT3'] + hd1['CRVAL3'])/1000.
    spec1 = np.zeros_like(ch1)
    for i in range(hd1['NAXIS3']-1):
        spec1[i] = np.mean(data1[i, :, :])

    data2 = fits.getdata(file2)
    hd2 = fits.getheader(file2)
    ch2 = np.linspace(0, hd2['NAXIS3']-1, num=hd2['NAXIS3']) 
    v2 = ((ch2 - hd2['CRPIX3'] + 1) * hd2['CDELT3'] + hd2['CRVAL3'])/1000.
    spec2 = np.zeros_like(ch2)
    for i in range(hd2['NAXIS3']-1):
        spec2[i] = np.mean(data2[i, :, :])
    fig = plt.figure(figsize = (8, 6))
    ax = fig.add_subplot()
    # ax.scatter(v1, spec1, color = 'red', label = r'$^{12}CO$')
    ax.plot(v1, spec1, ds = 'steps', color = 'blue', label = r'$^{12}CO$')
    ax.plot(v2, spec2*5, ds = 'steps', color = 'red', label=r'$^{13}CO\times5$')
    ax.set_xlabel(r'$\rm{Velocity\ (km\ s^{-1})}$')
    ax.set_ylabel(r'$\rm{T_{mb}}$')
    ax.legend(loc='upper left')
    plt.minorticks_on()
    fig.savefig(figname, bbox_inches = 'tight')
    plt.close()


if __name__=='__main__':
    path = './clouds/'
    clouds = os.listdir(path)
    for cld in clouds:
        plt_ave(path, cld)
