from astropy.io import fits
import numpy as np

# Use FITS file with radial component of magnetic field
fits_file = '/home/andre/cuboulder/abdullah-project/cgem.pdfi_output.101834.20120712_133000_TAI.Brll0.fits'

with fits.open(fits_file) as data:
    crpix1 = data[1].header['crpix1'] 
    crpix2 = data[1].header['crpix2'] 
    cdelt1 = data[1].header['cdelt1'] 
    cdelt2 = data[1].header['cdelt2'] 
    crval1 = data[1].header['crval1'] 
    crval2 = data[1].header['crval2'] 
    naxis1 = data[1].header['naxis1'] 
    naxis2 = data[1].header['naxis2'] 

x1b = np.empty(naxis1)
x1a = np.empty(naxis1+1)
for i in range(len(x1b)):
    x1b[i] = crval1 + cdelt1 * (i+1-crpix1)
    x1a[i] = crval1 + cdelt1 * (i+0.5-crpix1)
x1a[naxis1] = crval1 + cdelt1 * (naxis1+0.5-crpix1)

x2b = np.empty(naxis2)
x2a = np.empty(naxis2+1)
for j in range(len(x2b)):
    x2b[j] = crval2 + cdelt2 * (j+1-crpix2)
    x2a[j] = crval2 + cdelt2 * (j+0.5-crpix2)
x2a[naxis2] = crval2 + cdelt2 * (naxis2+0.5-crpix2)
# check 
print('difference between x1a and minlon value', x1a[0]-data[1].header['MINLON']*180./np.pi)
print('difference between x1a and maxlon value', x1a[-1]-data[1].header['MAXLON']*180./np.pi)
print('difference between x2a and maxcolat value', 90.-data[1].header['Maxcolat']*180./np.pi - x2a[0])
print('difference between x2a and mincolat value', 90.-data[1].header['Mincolat']*180./np.pi - x2a[-1])

# if Blat used
# x2a = np.empty(naxis2)
# x2b = np.empty(naxis2-1)
# for j in range(len(x2b)):
#     x2a[j] = crval2 + cdelt2 * (j+1-crpix2)
#     x2b[j] = crval2 + cdelt2 * (j+1.5-crpix2)
# x2a[naxis2-1] = crval2 + cdelt2 * (naxis2-crpix2)

# print('difference between x1b and minlon value', x1b[0]-data[1].header['MINLON']*180./np.pi)
# print('difference between x1b and maxlon value', x1b[-1]-data[1].header['MAXLON']*180./np.pi)

# print('difference between x2b and maxolat value', 90.-data[1].header['Maxcolat']*180./np.pi - x2b[0])
# print('difference between x2b and mincolat value', 90.-data[1].header['Mincolat']*180./np.pi - x2b[-1])

# print(data.info())
# print(data[0].header)
# print(data[1].header)