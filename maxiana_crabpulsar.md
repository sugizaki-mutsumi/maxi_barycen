# MAXI GSC pulsar timing analysis 

## References
- MAXI data analysis https://darts.isas.jaxa.jp/missions/maxi/analysis/dataanalysis.html
- NICER timing analysis https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/bary_corr_phase_resolved_spectra/


## Test environment 
- Heasoft Ver.6.34
- 

## Step by step procedure

### 1. Download data 
```
% mxdownload_wget --coordinates 83.633083,22.0145 -date_from 2010-01-01 -date_to 2010-01-10 -uri https://darts.isas.jaxa.jp/pub/maxi/mxdata/
```

### 2. Run standard analysis script mxproduct 
Use `mxproduct` in HEASOFT 
```
% mxproduct 83.633083 22.014500 2010-01-01 2010-01-10 object=crab
```

#### 2.1 Check Image
```
% ds9 products/crab_g_low.evt
... 
```
<img src="ds9.png" width="60%">
<!-- ![da9 image](ds9.png =100x100) -->

#### 2.2 Check light curve
Using python script [plot_lc.py](plot_lc.py), 
``` 
% python -i plot_lc.py
>>>
```
<img src="lc.png" width="60%">

#### 2.3 Check spectrum
```
% xspec
XSPEC12>data products/crab_g_low_src.pi 
XSPEC12>setplot ene
XSPEC12>ignore 1-**
XSPEC12>notice 2.0-20.
XSPEC12>model tbabs*power
XSPEC12>statistic cstat
XSPEC12>fit
XSPEC12>setplot rebin 3 10
XSPEC12>plot lda del
```
<img src="spec.png" width="60%">


### 3. Pulsar timing analysis 
#### 3.0 Create orbit file
Merge orbit files separated for each day into one file using a python script [`merge_orbfile.py`](merge_orbfile.py)
```
% merge_orbfile.py
```
#### 3.1 Barycentric time correction 
Sort event with time 
```
% ftsort products/crab_g_low.evt crab_g_low_tsort.evt TIME
```

Apply barycentric time correction with `barycen`
```
% barycen crab_g_low_tsort.evt crab_g_low_barycen.evt orbit.fits 83.633083 22.0145 orbext="ORBIT" orbform=COMPONENTS orbcol="X,Y,Z,vX,vY,vZ"
```

#### 3.2 Epoch-folding search and Folded pulse profile
Get Crab pulsar ephemeris

```
% efsearch cfile1="crab_g_low_barycen.evt" window="-" sepoch="15180 0.011470" dper=3.363725432368840e-02 dpdot=4.202665622939618e-13 nphase=16 nbint=INDEF dres=1e-10 nper=128
```
<img src="fes_nphase16.png" width="60%">


```
% efold nser=1 cfile1="crab_g_low_barycen.evt" window="-" sepoch="15180 0.011470" dper=3.363725432368840e-02 dpdot=4.202665622939618e-13 nphase=16 nbint=INDEF nintfm=INDEF
```
<img src="efold1.png" width="60%">


### 3.3 Phase resolved analysis

