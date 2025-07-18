#!/usr/bin/env python

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import astropy.io.fits as pyfits

lcfname = "products/crab_g_low.lc"
hdul = pyfits.open(lcfname)

hdu = hdul[1]
mjdref = hdu.header["MJDREFI"]+hdu.header["MJDREFF"]
elo1 = hdu.header["PI_LOW"]*0.05
ehi1 = hdu.header["PI_HI"]*0.05
vtc1 = hdu.data.field("TIME")
vts1 = hdu.data.field("START")
vte1 = hdu.data.field("STOP")
vrate1 = hdu.data.field("RATE")
vrerr1 = hdu.data.field("ERROR")

vmjd1 = vtc1/86400. + mjdref
xle1 = (vtc1-vts1)/86400.
xhe1 = (vte1-vtc1)/86400.
xle1 *= (xle1>0.)
xhe1 *= (xhe1>0.)

hdu = hdul[2]
mjdref = hdu.header["MJDREFI"]+hdu.header["MJDREFF"]
elo2 = hdu.header["PI_LOW"]*0.05
ehi2 = hdu.header["PI_HI"]*0.05
vtc2 = hdu.data.field("TIME")
vts2 = hdu.data.field("START")
vte2 = hdu.data.field("STOP")
vrate2 = hdu.data.field("RATE")
vrerr2 = hdu.data.field("ERROR")

vmjd2 = vtc2/86400. + mjdref
xle2 = (vtc2-vts2)/86400.
xhe2 = (vte2-vtc2)/86400.
xle2 *= (xle2>0.)
xhe2 *= (xhe2>0.)


fig, axs = plt.subplots(2, 1, sharex=True, figsize=(6,4))
fig.subplots_adjust(hspace=0.01, left=0.13, right=0.95, bottom=0.12, top=0.95)

axs[0].errorbar(vmjd1, vrate1, xerr=(xle1, xhe1), yerr=vrerr1, fmt='.', ms=3, lw=1, color='C0')
ymax = vrate1.max()*1.25
axs[0].set_ylim(0, ymax)
axs[0].set_ylabel("%.1lf-%.1lf keV\n Counts cm$^{-1}$ s$^{-1}$"%(elo1, ehi1))
axs[0].yaxis.set_label_coords(-0.08, 0.5)


axs[1].errorbar(vmjd2, vrate2, xerr=(xle2, xhe2), yerr=vrerr2, fmt='.', ms=3, lw=1, color='C1')
ymax = vrate2.max()*1.25
axs[1].set_ylim(0, ymax)
axs[1].set_ylabel("%.1lf-%.1lf keV\n Counts cm$^{-1}$ s$^{-1}$"%(elo2, ehi2))
axs[1].yaxis.set_label_coords(-0.08, 0.5)

axs[1].set_xlabel("MJD")
axs[1].xaxis.set_major_formatter(ticker.FormatStrFormatter('%d'))

plt.ion()
plt.show()
