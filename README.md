# rf-ccp
Common Conversion Point Stacking for receiver function (RF) analysis. Input data must be in SAC format. The header must contain BAZ, STLA and STLO, as well as an extra field where the ray parameter is stored.

Example parameter file:

    lat0,lon0 (initial latitude and longitude)
    lat1,lon1 (final latitude and longitude)
    t0 (initial time of data, in s)
    dx,dz (lateral and vertical resolution of profile, or cell size, in km)
    zmin,zmax (initial and final depth of profile, in km)
    hw (half width, or lateral sampling, in km)
    outfile (output file name, with full path)
    model (earth model file, with full path, in format Z,Vp,Vs and without header)
    p (name of ray parameter variable in SAC header. Must be in s/km)
    zv (exponential term for depth-amplitude scaling. As amplitude decreases with depth, it is sometimes useful to apply a scale factor. The amplitude (A) is recalculated as A=A*(z**zv), where z is depth. Leave as 0 if you do not want depth scaling at all)
    v (exponential term for phase weighting. 0 for linear stacking, i.e. no phase weighting)
    a (frequency of the RF)
    0 (stacking of first Fresnel zone) or 1 (stacking of ray paths)
    
Required libraries

    sacio
    fftw3
    
This program migrates Ps waves along the first fresnel zone and outputs a text file:
X (distance in km), Z (depth in km), A (amplitude), LAT, LON, Z (depth in deg), D (n. of stacked RFs in cell, or density)

To run the program:

    ccp [parameter file] [RF list file - with full paths]
    
RFs in the list must include the full path.
