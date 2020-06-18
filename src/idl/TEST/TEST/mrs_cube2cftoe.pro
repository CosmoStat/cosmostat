;+
; NAME:
;       mrs_cube2cftoe
;
; PURPOSE:
;       Estimate the fourier transform of a cube of data in a non regular grid
;       using Fourier sToeplitz
;
; CALLING:
;       cftoe = mrs_cube2cftoe(cube)
;
; INPUTS:
;       Cube --- structure of cube obtained with mrs_exact_split.pro
;    
; OUTPUTS:
;      cftoe ---  IDL structures with the following fields: 
;          cftoe : Fourier transform of the patches
;          nside : healpix map resolution
;          nx : nb of pixel per patch
;          nmaps : nb of patches
;          lat :  patch latitude
;          lon : patch longitude
;
; EXTERNAL CALLS:
;       mk_fourierToe
;
; EXAMPLE:
;          
; HISTORY:
;	Written by Sandrine Pires, Nov. 2009
;------------------------------------------------------------------------------
function mrs_cube2cftoe, Cube

psize_arcmin = pixel_size(Cube.nside) 
nmaps = Cube.nmaps
nfreq = fix(Cube.nx/4.)
NameCatIn= gettmpfilename()
NameImagOut = gettmpfilename()

cftoe = complexarr(nfreq*2.+1, nfreq*2.+1, nmaps)

for l = 0l, nmaps-1 do begin
  print, 'map =', l
  writefits, NameCatIn+'.fits', Cube.map(*,0:Cube.npix_map[l]-1,l)
  com='mk_fourierToe -v -N'+strtrim(string(nfreq),1)+' '+ NameCatIn+'.fits' +' '+ NameImagOut
  spawn, com  
  tfoe_re=readfits(NameImagOut+'_re.fits')
  tfoe_im=readfits(NameImagOut+'_im.fits')
  cftoe(*, *, l) = complex(tfoe_re, tfoe_im)
endfor

Cubeftoe={cftoe:cftoe, nside:Cube.nside, nx:Cube.nx, nmaps:nmaps, lat:Cube.lat, lon:Cube.lon}
delete, NameCatIn+'.fits'
delete, NameImagOut+'_re.fits'
delete, NameImagOut+'_im.fits'
return, Cubeftoe
end
