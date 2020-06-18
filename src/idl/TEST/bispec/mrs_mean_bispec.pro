;+
; NAME:
;    mrs_mean_bispec
;
; PURPOSE:
;    Compute the mean bispectrum
;
; CALLING:
;    mrs_mean_bispec, Imag, Nmap, iso=iso, tab_bispec, mean_bispec, var_bispec, x  
;  
; INPUTS:
;    Imag --- Healpix Map (default is nested)  
;    Nmap --- Imag is splitted in maps Nmap x Nmap
;    
; OUTPUTS:
;    tab_bispec --- array with all the bispectrum
;    mean_bispec --- mean bispectrum (complex)
;    var_bispec --- std of the bispectrum (complex)
;    x --- x-axis (frequencies)
;
; OPTIONS:
;    iso --- angle in degrees (defaul is 120 : equilateral) 
;         
; HISTORY:
;    Written: Sandrine Pires in September 2009 
;-------------------------------------------------------------------
pro mrs_mean_bispec, Map, Nmap, iso=iso, tab_bispec, mean_bispec, var_bispec, x, Name=Name

if not keyword_set(iso) then iso = 120
Cube = mrs_split(Map, nx=Nmap)
sz = size(Cube.map)
n3 = sz[3]
print, 'Nb maps = ', n3 
n1 = fix(Nmap*3./4.)

tab_bispec = complexarr(n1, n3)
mean_bispec = complexarr(n1)
var_bispec = complexarr(n1)
for i=0, n3-1 do begin
  print, 'Computing the',i+1, 'eme bispectre'
  ;bispec = polfft_bispec_eq(Cube.map(*,*,i), x=x)
  bispec = polfft_bispec_iso(Cube.map(*,*,i), x=x, iso=iso, Name=Name)
  tab_bispec(*, i) = bispec
endfor

for i = 0, n1-1 do begin 
  mean_bispec(i) = mean(tab_bispec(i, *))
  var_bispec(i) = sigma(tab_bispec(i, *))
endfor

end


