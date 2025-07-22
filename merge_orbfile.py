#!/usr/bin/env python

import os, sys
import datetime
from astropy.time import Time, TimeDelta
import astropy.io.fits as pyfits

import argparse


def get_orbfilelist(mjd_start, mjd_stop, obsdir='obs') :
    orbfilelist = []
    for mjd in range(mjd_start, mjd_stop+1) :
        pardir = obsdir+os.sep+"MJD%d000"%(int(mjd/1000))
        if not os.path.isdir(pardir) : continue
        daydir = pardir+os.sep+"MJD%d"%(mjd)
        if not os.path.isdir(daydir) : continue
        obsfname = daydir+os.sep+"auxil"+os.sep+"mx_mjd%d.orb.gz"%(mjd)
        if os.path.isfile(obsfname) :
            orbfilelist.append(obsfname)
    return orbfilelist
    

def create_orbfile(orbfilelist, outfname) :

    if os.path.isfile(outfname) :
        print("Warning: outfile %s has been already created"%(outfname))
        return

    listfname = "temp_orbfile.list"
    outfile = open(listfname, "wt")
    for orbfname in orbfilelist :
        outfile.write("%s[orbit][#row<#naxis2]\n"%(orbfname))
    outfile.close()
    cmd = "ftmerge @%s %s"%(listfname, outfname)
    print(cmd)
    os.system(cmd)

    if not os.path.isfile(outfname) :
        print("Error: outfile %s has not been created"%(outfname))
        return


    ### update header  ####
    hdul = pyfits.open(outfname, mode='update')

    hdu = hdul["ORBIT"]
    #print(hdu.data)
    vtime = hdu.data.field("TIME")
    tstart = vtime[0]
    tstop  = vtime[-1]

    print("Data tstart = ", tstart)
    print("Data tstop  = ", tstop)


    tstart_hdr = hdu.header["TSTART"]
    tstop_hdr  = hdu.header["TSTOP"]
    print("HDU Header tstart = ", tstart_hdr)
    print("HDU Header tstop  = ", tstop_hdr)
    ### update
    hdu.header["TSTART"] = tstart
    hdu.header["TSTOP"]  = tstop


    mission_time0 = Time('2000-01-01T00:00:00', format='isot', scale='tt')
    date_obs = mission_time0 + TimeDelta(tstart, format='sec')
    date_end = mission_time0 + TimeDelta(tstop,  format='sec')
    hdu.header["DATE-OBS"] = date_obs.isot
    hdu.header["DATE-END"]  = date_end.isot
    print("HDU Header DATE-OBS = ", date_obs.isot)
    print("HDU Header DATE-END  = ", date_end.isot)

    ### update file creation date
    timenow = Time.now()
    for hdu in hdul :
        hdu.header["DATE"]  = timenow.isot
    
    hdul.flush()
    hdul.close()

    return
    

if __name__=="__main__" :
  
    #mjd_start = 55000
    #mjd_stop  = int(Time.now().mjd)
    date_from_default = '2009-08-10T00:00:00'
    date_to_default   =  Time.now().fits

    parser = argparse.ArgumentParser(description='merge orbit files in a MAXI obs data directry')

    # Optional argument with a short and long form
    parser.add_argument("-o", "--output_file",
                        type=str, default="orbit.fits",
                        help="Output filename (default: orbit.fits)")

    parser.add_argument("-f", "--date_from",
                        type=str, default=date_from_default,
                        help="Start date YYYY-mm-dd or mjd (integer)")

    parser.add_argument("-t", "--date_to",
                        type=str, default=date_to_default,
                        help="Start date YYYY-mm-dd or mjd (integer)")
    
    args = parser.parse_args()

    
    #print(args)
    outfname = args.output_file
    
    date_from = args.date_from 
    try  :
        astime_from = Time(date_from, format='fits')
        mjd_from = int(astime_from.mjd)
    except ValueError :
        try :
            mjd_from = int(date_from)
        except :
            mjd_from = ValueError
        
    date_to = args.date_to
    try  :
        astime_to = Time(date_to, format='fits')
        mjd_to = int(astime_to.mjd)
    except ValueError :
        try :
            mjd_to = int(date_to)
        except :
            mjd_to = ValueError
        
    if ValueError in [mjd_from, mjd_to] :
        print("Error: date_from=%s or date_to=%s is illegal")
        sys.exit()
            
    print("output filename = %s, mjd_from=%d, mjd_to=%d"%(outfname, mjd_from, mjd_to))
    if os.path.isfile(outfname) :
        print("Error: output file = %s already exists"%(outfname))              
        sys.exit()

    ### run merge
    orbfilelist =  get_orbfilelist(mjd_from, mjd_to, obsdir='obs') 
    print("Merging orbit files")
    for fname in orbfilelist: 
        print(fname)
    if len(orbfilelist)>0 :
        create_orbfile(orbfilelist, outfname)
    

