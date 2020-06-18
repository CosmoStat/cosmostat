;+
;  NAME: 
;        MRS_LATPOWSPEC
;
;  PURPOSE: 
;        Estimate the power spectrum per latitude from an Healpix map         
;
; CALLING: 
;        mrs_latpowspec, Map, field_deg=field_deg
;
; INPUTS: 
;	       Map --- Healpix map
;	
; OUTPUTS: 
;       IDL structure cl_g: 
;             l --- modes
;             mcl --- mean power spectrum
;             vcl --- variance on the mean power spectrum
;             tcl --- power spectrum per patch
;             mlcl --- power spectrum per latitude 
;             vlcl --- variance on the power spectrum per latitude
;             nlat --- number of latitudes
;             lat --- array of latitudes
;             nlati --- index of patches per latitude
;  KEYWORDS: 
;       field_deg --- size of the patches in degrees (default is 10°)        
;       wind --- window function ('rect', 'hanning', 'hamming', 
;                                'bartlett','cosine')
; HISTORY:
;	Written: Sandrine Pires 2010.
;-
;-------------------------------------------------------------------------------
function mrs_latpowspec, Map, field_deg=field_deg, wind=wind, proj=proj

nside = gnside(Map)
if not keyword_set(wind) then wind = 'hanning'
if not keyword_set(proj) then proj = 'gnome'
if not keyword_set(field_deg) then field_deg = 10.
psize_arcmin = pixel_size(nside)
nx0 = fix(field_deg*60./psize_arcmin)
cube_g = mrs_exact_split(Map, nx=nx0, wind=wind, proj=proj)
fft_g = mrs_cube2ftoe(cube_g)
cl_g = mrs_cftoe2latspec(fft_g, wind=wind)

return, cl_g
end

