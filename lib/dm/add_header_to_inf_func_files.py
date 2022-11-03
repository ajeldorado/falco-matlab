"""Add new header values to the influence function FITS files to function
after prop_dm.py in PROPER changed how the header values are used."""
from astropy.io import fits
import os

HERE = os.path.abspath(os.path.dirname(__file__))

fnList = [
    'influence_BMC_2kDM_400micron_res10.fits',
    'influence_BMC_2kDM_400micron_res20.fits',
    'influence_BMC_kiloDM_300micron_res10_spline.fits',
    'influence_dm5v2.fits',
    'influence_dm_FEA_10milsV2_20160330.fits',
]

for index, fn in enumerate(fnList):
    
    print(fn)
    fnFull = os.path.join(HERE, fn)
    inf_func = fits.getdata(fnFull)
    hdul = fits.open(fnFull)  # open a FITS file
    hdr = hdul[0].header  # the primary HDU header
    P2PD_M = hdr['P2PDX_M']
    C2CD_M = hdr['C2CDX_M']

    # Output
    hdu = fits.PrimaryHDU(inf_func)
    hdu.header['P2PD_M'] = (P2PD_M, 'pixel2pix distance in meters')
    hdu.header['P2PDX_M'] = P2PD_M
    hdu.header['P2PDY_M'] = P2PD_M
    hdu.header['C2CD_M'] = (C2CD_M, 'center2cen dist of actuators in meters')
    hdu.header['C2CDX_M'] = C2CD_M
    hdu.header['C2CDY_M'] = C2CD_M

    # fnOut = os.path.join(HERE, 'test_file.fits')
    fnOut = fnFull
    hdu.writeto(fnOut, overwrite=True)
    