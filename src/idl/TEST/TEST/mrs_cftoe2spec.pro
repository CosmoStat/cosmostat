;+
; NAME:
;       mrs_cftoe2spec
;
; PURPOSE:
;       Estimate a mean cl from N patches of Fourier Toe
;
; CALLING:
;
;       Cl = mrs_cftoe2spec(CFToe)
;
; INPUTS:
;       CFToe -- cube of Fourier Toeplitz
;    
; OUTPUTS:
;      Cl ---  IDL structures with the following fields: 
;          l : modes
;          mcl : mean cl
;          vcl : variance on cl
;          mlcl : mean cl per latitude
;          vlcl : variance on cl per latitude
;          nlat : nb de latitude
;          latitude : latitude values
;          nlati : latitude index 
;
; KEYWORD:
;      whanning --- correct from a hanning window
;
; EXTERNAL CALLS:
;    fft_spec.cc
;
; EXAMPLE:
;          
; HISTORY:
;	Written by Sandrine Pires, Nov. 2009
;------------------------------------------------------------------------------
function mrs_cftoe2spec, CFToe, whanning=whanning

a = systime(1)

Name = gettmpfilename()
print, Name
psize_arcmin = pixel_size(CFToe.nside) 
nmaps = CFToe.nmaps
lat = CFToe.lat

ilat = 0
lat_i = -1
for l=0l, nmaps-1 do begin
  if lat[l] ne lat_i then begin
    ilat = ilat + 1
    lat_i = lat[l]
  endif
endfor
nlat = ilat
print, 'nlat=', nlat

ilat = 0
lat_i = -1
nlati=fltarr(nlat+1)
latitude=fltarr(nlat)
for l=0l, nmaps-1 do begin
  if lat[l] ne lat_i then begin
    lat_i = lat[l]
    latitude[ilat]= lat[l]
    nlati(ilat) = l 
    ilat = ilat + 1
  endif
endfor
nlati(nlat)=nmaps
nfreq = CFToe.nx/4.
field_rad = (psize_arcmin/60.)*float(CFToe.nx)*!dtor
Nh = CFToe.nx
w = hanning(Nh, Nh, alpha=0.5)
if keyword_set(whanning) then weight = Nh^2/total(w^2) else weight = 1.

for l=0l, nmaps-1 do begin
  print, 'l= ', l
  writefits, Name+'_re.fits', re(CFtoe.cftoe(*,*,l))
  writefits, Name+'_im.fits', im(CFToe.cftoe(*,*,l))
  com = 'fft_spec '+Name +' '+Name+'_spec.fits'
  spawn, com
  spec=readfits(Name+'_spec.fits')
  if l eq 0 then begin
    sz = size(spec)
    Nsample = sz[1]
    som = fltarr(Nsample)
    x= fltarr(Nsample)
    x= spec(*, 0)
    tab_spec = dblarr(Nsample, nmaps)
  endif
  tab_spec(*, l) = spec(*,1)*weight
endfor

mean_spec = dblarr(Nsample)
var_spec = dblarr(Nsample)
for i=0l, Nsample-1 do begin
  mean_spec(i) = mean(tab_spec(i, *))
  var_spec(i) = sigma(tab_spec(i, *))
endfor

mean_spec_lat = dblarr(Nsample, nlat)
var_spec_lat = dblarr(Nsample, nlat)
for j= 0, nlat-1 do begin
  for i=0l, Nsample-1 do begin
    mean_spec_lat(i,j) = mean(tab_spec(i, nlati(j):nlati(j+1)-1))
    var_spec_lat(i, j) = sigma(tab_spec(i, nlati(j):nlati(j+1)-1))
  endfor
endfor



;normalisation 
l = x * (2.*!pi*float(nfreq*2.+1))/float(field_rad)
mcl = mean_spec*(field_rad)^2/float(nfreq*2.+1)^2
vcl = var_spec*(field_rad)^2/float(nfreq*2.+1)^2
mlcl = mean_spec_lat*(field_rad)^2/float(nfreq*2.+1)^2
vlcl = var_spec_lat*(field_rad)^2/float(nfreq*2.+1)^2

tcl = tab_spec*(field_rad)^2/float(nfreq*2.+1)^2

;pow = {x:x, l:l, mean_spec:mean_spec, mcl:mcl, var_spec:var_spec, vcl:vcl, tab_spec:tab_spec, tcl:tcl,mlcl:mlcl,vlcl:vlcl, nlat:nlat, latitude:latitude, nlati:nlati}
pow = {l:l, mcl:mcl, vcl:vcl, tcl:tcl, mlcl:mlcl, vlcl:vlcl, nlat:nlat, latitude:latitude, nlati:nlati}

delete, Name+'_re.fits'
delete, Name+'_im.fits'
delete, 'spec_'+Name+'.fits'

print, 'time =', systime(1)-a
return, pow
end
