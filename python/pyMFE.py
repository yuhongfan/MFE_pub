    # *****************************************************************************
    # NAME: pyMFE
    # DESCRIPTION: Module contains functions to read data from MFE output binary .dat files,
    # and provide vars in CGS units, and compute electric current, and other functions
    # List of functions:
    #   readphysparams (Read physparams.dat file)
    #   readgrid (Read grid.dat file)
    #   readdata (Read output data file (dat files contain data with ghost zones))
    #   getdata (Read MFE data arrays in CGS units inside domain (no ghost zones))
    #   coords_from_fits (Compute coordinate arrays from JSOC FITS file of PDFI Bradial) 
    # VERSION: 1.2 (24 Jul 2023) Last update: added getting coodrinate arrays from JSOC FITS file
    # USAGE: import pyMFE as mfe
    # AUTHOR: Dr. Andrey Afanasyev, LASP CU Boulder
    # COPYRIGHT:
    # BUGS REPORT: andrei.afanasev@colorado.edu
    # *****************************************************************************

import numpy as np
from astropy.io import fits

def readphysparams(filename):

    # DESCRIPTION: Read physparams.dat file
    # USAGE:
    # g_const, rgas, kboltz, mproton, pi, gamma,\
    #         unit_rho, unit_len, unit_b, unit_temp, unit_v, unit_time,\
    #             msol, r_c, temp_c, rho_c, muconst = mfe.readphysparams(physparamsfile)

    print('reading ', filename)

    file = open(filename, 'rb')
    dat = np.fromfile(file, dtype='float64', count=-1)
    
    g_const, rgas, kboltz, mproton, pi, gamma,\
        unit_rho, unit_len, unit_b, unit_temp, unit_v, unit_time,\
            msol, r_c, temp_c, rho_c, muconst = dat
    
    file.close()
    
    return g_const, rgas, kboltz, mproton, pi, gamma,\
        unit_rho, unit_len, unit_b, unit_temp, unit_v, unit_time,\
            msol, r_c, temp_c, rho_c, muconst


def readgrid(filename):

    # DESCRIPTION: Read grid.dat file
    # USAGE:
    # inn, jn, kn, x1a, x1b, x2a, x2b, x3a, x3b, dx1a, dx1b, dx2a, dx2b, dx3a, dx3b,\
    #     g2a, g2b, g31a, g31b, g32a, g32b, dg2bd1, dg2ad1, dg31bd1, dg31ad1, dg32bd2, dg32ad2,\
    #         dvl1a, dvl1b, dvl2a, dvl2b, dvl3a, dvl3b, dx1ai, dx1bi, dx2ai, dx2bi, dx3ai, dx3bi,\
    #             g2ai, g2bi, g31ai, g31bi, g32ai, g32bi, dvl1ai, dvl1bi, dvl2ai, dvl2bi, dvl3ai, dvl3bi = \
    #                 mfe.readgrid(gridfile)

    print('reading ', filename)

    file = open(filename, 'rb')
    dat = np.fromfile(file, dtype=np.dtype('i4'), count=3)
    inn, jn, kn = dat
    del dat

    ist = 3-1
    js = 3-1
    ks = 3-1
    ie = ist+(inn-6)
    je = js+(jn-6)
    ke = ks+(kn-6)
    # ism1 = ist-1
    # jsm1 = js-1
    # ksm1 = ks-1
    ism2 = ist-2
    jsm2 = js-2
    ksm2 = ks-2
    # iep1 = ie+1
    # jep1 = je+1
    # kep1 = ke+1
    iep2 = ie+2
    jep2 = je+2
    kep2 = ke+2
    iep3 = ie+3
    jep3 = je+3
    kep3 = ke+3

    x1asize = iep3-ism2+1
    x1bsize = iep2-ism2+1
    x2asize = jep3-jsm2+1
    x2bsize = jep2-jsm2+1
    x3asize = kep3-ksm2+1
    x3bsize = kep2-ksm2+1

    dx1asize = iep2-ism2+1
    dx1bsize = iep3-ism2+1
    dx2asize = jep2-jsm2+1
    dx2bsize = jep3-jsm2+1
    dx3asize = kep2-ksm2+1
    dx3bsize = kep3-ksm2+1

    g2asize = iep3-ism2+1
    g2bsize = iep2-ism2+1
    g31asize = iep3-ism2+1
    g31bsize = iep2-ism2+1
    g32asize = jep3-jsm2+1
    g32bsize = jep2-jsm2+1

    dg2bd1size = iep3-ism2+1
    dg2ad1size = iep2-ism2+1
    dg31bd1size = iep3-ism2+1
    dg31ad1size = iep2-ism2+1
    dg32bd2size = jep3-ism2+1
    dg32ad2size = jep2-ism2+1

    dvl1asize = iep2-ism2+1
    dvl1bsize = iep3-ism2+1
    dvl2asize = jep2-jsm2+1
    dvl2bsize = jep3-jsm2+1
    dvl3asize = kep2-ksm2+1
    dvl3bsize = kep3-ksm2+1

    dx1aisize = iep2-ism2+1
    dx1bisize = iep3-ism2+1
    dx2aisize = jep2-jsm2+1
    dx2bisize = jep3-jsm2+1
    dx3aisize = kep2-ksm2+1
    dx3bisize = kep3-ksm2+1

    g2aisize = iep3-ism2+1
    g2bisize = iep2-ism2+1
    g31aisize = iep3-ism2+1
    g31bisize = iep2-ism2+1
    g32aisize = jep3-jsm2+1
    g32bisize = jep2-jsm2+1

    dvl1aisize = iep2-ism2+1
    dvl1bisize = iep3-ism2+1
    dvl2aisize = jep2-jsm2+1
    dvl2bisize = jep3-jsm2+1
    dvl3aisize = kep2-ksm2+1
    dvl3bisize = kep3-ksm2+1

    dt = np.dtype(
        str(x1asize)+'f8,' + str(x1bsize)+'f8,' + str(x2asize)+'f8,' + 
        str(x2bsize)+'f8,' + str(x3asize)+'f8,' + str(x3bsize)+'f8,' +

        str(dx1asize)+'f8,' + str(dx1bsize)+'f8,' + str(dx2asize)+'f8,' +
        str(dx2bsize)+'f8,' + str(dx3asize)+'f8,' + str(dx3bsize)+'f8,' +

        str(g2asize)+'f8,' + str(g2bsize)+'f8,' + str(g31asize)+'f8,' +
        str(g31bsize)+'f8,' + str(g32asize)+'f8,' + str(g32bsize)+'f8,' +
        
        str(dg2bd1size)+'f8,' + str(dg2ad1size)+'f8,' + str(dg31bd1size)+'f8,' +
        str(dg31ad1size)+'f8,' + str(dg32bd2size)+'f8,' + str(dg32ad2size)+'f8,' +
        
        str(dvl1asize)+'f8,' + str(dvl1bsize)+'f8,' + str(dvl2asize)+'f8,' +
        str(dvl2bsize)+'f8,' + str(dvl3asize)+'f8,' + str(dvl3bsize)+'f8,' +
        
        str(dx1aisize)+'f8,' + str(dx1bisize)+'f8,' + str(dx2aisize)+'f8,' +
        str(dx2bisize)+'f8,' + str(dx3aisize)+'f8,' + str(dx3bisize)+'f8,' +
        
        str(g2aisize)+'f8,' + str(g2bisize)+'f8,' + str(g31aisize)+'f8,' +
        str(g31bisize)+'f8,' + str(g32aisize)+'f8,' + str(g32bisize)+'f8,' +
        
        str(dvl1aisize)+'f8,' + str(dvl1bisize)+'f8,' + str(dvl2aisize)+'f8,' +
        str(dvl2bisize)+'f8,' + str(dvl3aisize)+'f8,' + str(dvl3bisize)+'f8'
        )

    dat = np.fromfile(file, dtype=dt, count=1)

    x1a, x1b, x2a, x2b, x3a, x3b, dx1a, dx1b, dx2a, dx2b, dx3a, dx3b,\
        g2a, g2b, g31a, g31b, g32a, g32b, dg2bd1, dg2ad1, dg31bd1, dg31ad1, dg32bd2, dg32ad2,\
            dvl1a, dvl1b, dvl2a, dvl2b, dvl3a, dvl3b, dx1ai, dx1bi, dx2ai, dx2bi, dx3ai, dx3bi,\
                g2ai, g2bi, g31ai, g31bi, g32ai, g32bi, dvl1ai, dvl1bi, dvl2ai, dvl2bi, dvl3ai, dvl3bi = \
                    [dat[0][i] for i in range(48)]

    file.close()

    return inn, jn, kn, x1a, x1b, x2a, x2b, x3a, x3b, dx1a, dx1b, dx2a, dx2b, dx3a, dx3b,\
        g2a, g2b, g31a, g31b, g32a, g32b, dg2bd1, dg2ad1, dg31bd1, dg31ad1, dg32bd2, dg32ad2,\
            dvl1a, dvl1b, dvl2a, dvl2b, dvl3a, dvl3b, dx1ai, dx1bi, dx2ai, dx2bi, dx3ai, dx3bi,\
                g2ai, g2bi, g31ai, g31bi, g32ai, g32bi, dvl1ai, dvl1bi, dvl2ai, dvl2bi, dvl3ai, dvl3bi


def readdata(filename):

    # DESCRIPTION: Read output data file (dat files contain data with ghost zones)
    # USAGE: data, time, inn, jn, kn = mfe.readdata(datafile)

    print('reading ', filename)

    file = open(filename, 'rb')
    dat = np.fromfile(file, dtype=np.dtype('f8,i4,i4,i4'), count=1)
    time, inn, jn, kn = dat[0]
    del dat

    ist = 3-1
    js = 3-1
    ks = 3-1
    ie = ist+(inn-6)
    je = js+(jn-6)
    ke = ks+(kn-6)
    # ism1 = ist-1
    # jsm1 = js-1
    # ksm1 = ks-1
    ism2 = ist-2
    jsm2 = js-2
    ksm2 = ks-2
    # iep1 = ie+1
    # jep1 = je+1
    # kep1 = ke+1
    iep2 = ie+2
    jep2 = je+2
    kep2 = ke+2
    # iep3 = ie+3
    # jep3 = je+3
    # kep3 = ke+3

    datasize = (iep2-ism2+1,jep2-jsm2+1,kep2-ksm2+1)
    sizeflat = datasize[0]*datasize[1]*datasize[2]
    dat = np.fromfile(file, dtype=np.dtype(str(sizeflat)+'f8'), count=1)
    data = dat[0]
    data = data.reshape(datasize[::-1]).T

    file.close()

    return data, time, inn, jn, kn


def getdata(var, it):

    # DESCRIPTION: Read MFE data arrays in CGS units inside domain (no ghost zones)
    # USAGE: data_cgs, x1_cgs, x2_cgs, x3_cgs, time_cgs = mfe.getdata(varname, snapshot_number)

    physparamsfile = './physparams.dat'
    _, _, kboltz, mproton, _, gamma,\
            unit_rho, unit_len, unit_b, _, unit_v, unit_time,\
                _, _, _, _, muconst = readphysparams(physparamsfile)
    rstar = kboltz / muconst / mproton

    gridfile = './grid.dat'
    inn, jn, kn, x1a, x1b, x2a, x2b, x3a, x3b, _, dx1b, _, dx2b, _, dx3b,\
        _, g2b, _, g31b, _, g32b, _, _, _, _, _, _,\
            _, _, _, _, _, _, _, _, _, _, _, _,\
                g2ai, g2bi, g31ai, g31bi, g32ai, g32bi, _, _, _, _, _, _ = \
                    readgrid(gridfile)

    ist = 3-1
    js = 3-1
    ks = 3-1
    ie = ist+(inn-6)
    je = js+(jn-6)
    ke = ks+(kn-6)
    # ism1 = ist-1
    # jsm1 = js-1
    # ksm1 = ks-1
    # ism2 = ist-2
    # jsm2 = js-2
    # ksm2 = ks-2
    iep1 = ie+1
    jep1 = je+1
    kep1 = ke+1
    # iep2 = ie+2
    # jep2 = je+2
    # kep2 = ke+2
    # iep3 = ie+3
    # jep3 = je+3
    # kep3 = ke+3

    ddatafile    = './d_'  + str(it).zfill(4) +'.dat'
    edatafile    = './e_'  + str(it).zfill(4) +'.dat'
    v1datafile   = './v1_' + str(it).zfill(4) +'.dat'
    v2datafile   = './v2_' + str(it).zfill(4) +'.dat'
    v3datafile   = './v3_' + str(it).zfill(4) +'.dat'
    b1datafile   = './b1_' + str(it).zfill(4) +'.dat'
    b2datafile   = './b2_' + str(it).zfill(4) +'.dat'
    b3datafile   = './b3_' + str(it).zfill(4) +'.dat'

    if var == 'd':
        d, time, _, _, _ = readdata(ddatafile)
        data_cgs = d[ist:ie+1, js:je+1, ks:ke+1] * unit_rho
        x1_cgs = x1b[ist:ie+1] * unit_len
        x2_cgs = x2b[js:je+1] * unit_len
        x3_cgs = x3b[ks:ke+1] * unit_len

    elif var == 'pres':
        e, time, _, _, _ = readdata(edatafile)
        data_cgs = e[ist:ie+1, js:je+1, ks:ke+1] * (gamma-1.0) * unit_rho * unit_v * unit_v
        x1_cgs = x1b[ist:ie+1] * unit_len
        x2_cgs = x2b[js:je+1] * unit_len
        x3_cgs = x3b[ks:ke+1] * unit_len

    elif var == 'temp':
        d, time, _, _, _ = readdata(ddatafile)
        e, _, _, _, _ = readdata(edatafile)
        rho_cgs = d * unit_rho
        pres_cgs = e * (gamma-1.0) * unit_rho * unit_v * unit_v
        data_cgs = (pres_cgs / rstar / rho_cgs)[ist:ie+1, js:je+1, ks:ke+1]
        x1_cgs = x1b[ist:ie+1] * unit_len
        x2_cgs = x2b[js:je+1] * unit_len
        x3_cgs = x3b[ks:ke+1] * unit_len
    
    elif var == 'v1':
        v1, time, _, _, _ = readdata(v1datafile)
        data_cgs = v1[ist:iep1+1, js:je+1, ks:ke+1] * unit_v
        x1_cgs = x1a[ist:iep1+1] * unit_len
        x2_cgs = x2b[js:je+1] * unit_len
        x3_cgs = x3b[ks:ke+1] * unit_len

    elif var == 'v2':
        v2, time, _, _, _ = readdata(v2datafile)
        data_cgs = v2[ist:ie+1, js:jep1+1, ks:ke+1] * unit_v
        x1_cgs = x1b[ist:ie+1] * unit_len
        x2_cgs = x2a[js:jep1+1] * unit_len
        x3_cgs = x3b[ks:ke+1] * unit_len

    elif var == 'v3':
        v3, time, _, _, _ = readdata(v3datafile)
        data_cgs = v3[ist:ie+1, js:je+1, ks:kep1+1] * unit_v
        x1_cgs = x1b[ist:ie+1] * unit_len
        x2_cgs = x2b[js:je+1] * unit_len
        x3_cgs = x3a[ks:kep1+1] * unit_len

    elif var == 'b1':
        b1, time, _, _, _ = readdata(b1datafile)
        data_cgs = b1[ist:iep1+1, js:je+1, ks:ke+1] * unit_b
        x1_cgs = x1a[ist:iep1+1] * unit_len
        x2_cgs = x2b[js:je+1] * unit_len
        x3_cgs = x3b[ks:ke+1] * unit_len
    
    elif var == 'b2':
        b2, time, _, _, _ = readdata(b2datafile)
        data_cgs = b2[ist:ie+1, js:jep1+1, ks:ke+1] * unit_b
        x1_cgs = x1b[ist:ie+1] * unit_len
        x2_cgs = x2a[js:jep1+1] * unit_len
        x3_cgs = x3b[ks:ke+1] * unit_len

    elif var == 'b3':
        b3, time, _, _, _ = readdata(b3datafile)
        data_cgs = b3[ist:ie+1, js:je+1, ks:kep1+1] * unit_b
        x1_cgs = x1b[ist:ie+1] * unit_len
        x2_cgs = x2b[js:je+1] * unit_len
        x3_cgs = x3a[ks:kep1+1] * unit_len

    elif var == 'j1':
        b2, time, _, _, _ = readdata(b2datafile)
        b3, _, _, _, _ = readdata(b3datafile)
        x1_cgs = x1b[ist:ie+1] * unit_len
        x2_cgs = x2a[js:jep1+1] * unit_len
        x3_cgs = x3a[ks:kep1+1] * unit_len
        # j1
        tmp = (b3.transpose(0,2,1) * g32b).transpose(0,2,1)
        tmp = (tmp[ist:ie+1, js:jep1+1, ks:kep1+1] - tmp[ist:ie+1, js-1:jep1, ks:kep1+1])
        tmp = (tmp.transpose(0,2,1) / dx2b[js:jep1+1] * g32ai[js:jep1+1]).transpose(0,2,1)
        tmp = (tmp.transpose(2,1,0) * g2bi[ist:ie+1]).transpose(2,1,0)
        j1 = tmp
        tmp = (b2[ist:ie+1, js:jep1+1, ks:kep1+1] - b2[ist:ie+1, js:jep1+1, ks-1:kep1]) / dx3b[ks:kep1+1]
        tmp = (tmp.transpose(2,1,0) * g31bi[ist:ie+1]).transpose(2,1,0) 
        tmp = (tmp.transpose(0,2,1) * g32ai[js:jep1+1]).transpose(0,2,1)
        j1 = j1 - tmp
        data_cgs = j1 * unit_b / unit_len

    elif var == 'j2':
        b1, time, _, _, _ = readdata(b1datafile)
        b3, _, _, _, _ = readdata(b3datafile)
        x1_cgs = x1a[ist:iep1+1] * unit_len
        x2_cgs = x2b[js:je+1] * unit_len
        x3_cgs = x3a[ks:kep1+1] * unit_len
        # j2
        tmp = (b3.transpose(2,1,0) * g31b).transpose(2,1,0)
        tmp = tmp[ist-1:iep1, js:je+1, ks:kep1+1] - tmp[ist:iep1+1, js:je+1, ks:kep1+1]
        tmp = (tmp.transpose(2,1,0) / dx1b[ist:iep1+1] * g31ai[ist:iep1+1]).transpose(2,1,0)
        j2 = tmp
        tmp = (b1[ist:iep1+1, js:je+1, ks:kep1+1] - b1[ist:iep1+1, js:je+1, ks-1:kep1]) / dx3b[ks:kep1+1]
        tmp = (tmp.transpose(0,2,1) * g32bi[js:je+1]).transpose(0,2,1)
        tmp = (tmp.transpose(2,1,0) * g31ai[ist:iep1+1]).transpose(2,1,0)
        j2 = j2 + tmp
        data_cgs = j2 * unit_b / unit_len

    elif var == 'j3':
        b1, time, _, _, _ = readdata(b1datafile)
        b2, _, _, _, _ = readdata(b2datafile)
        x1_cgs = x1a[ist:iep1+1] * unit_len
        x2_cgs = x2a[js:je+1] * unit_len
        x3_cgs = x3b[ks:ke+1] * unit_len
        # j3
        tmp = (b2.transpose(2,1,0) * g2b).transpose(2,1,0)
        tmp = tmp[ist:iep1+1, js:jep1+1, ks:ke+1] - tmp[ist-1:iep1, js:jep1+1, ks:ke+1]
        tmp = (tmp.transpose(2,1,0) / dx1b[ist:iep1+1] * g2ai[ist:iep1+1]).transpose(2,1,0)
        j3 = tmp
        tmp = b1[ist:iep1+1, js-1:je+1, ks:ke+1] - b1[ist:iep1+1, js:jep1+1, ks:ke+1] 
        tmp = (tmp.transpose(0,2,1) / dx2b[js:jep1+1]).transpose(0,2,1)
        tmp = (tmp.transpose(2,1,0) * g2ai[ist:iep1+1]).transpose(2,1,0)
        j3 = j3 + tmp
        data_cgs = j3 * unit_b / unit_len

    time_cgs = time * unit_time

    return data_cgs, x1_cgs, x2_cgs, x3_cgs, time_cgs


def coords_from_fits(radial_Bfield_fits_file):

    # DESCRIPTION: Compute coordinate arrays for MFE code using the JSOC FITS file of RADIAL component of
    #   PDFI magnetic field
    # USAGE: x1b, x1a, x2b, x2a = coords_from_fits(radial_Bfield_fits_file)

    with fits.open(radial_Bfield_fits_file) as data:
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

    return x1b, x1a, x2b, x2a