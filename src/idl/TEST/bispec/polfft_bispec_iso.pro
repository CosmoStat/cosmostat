;+
; NAME: 
;       POLFFT_BISPEC_ISO
;
; PURPOSE: 
;       Calculate a particular configuration of the bispectrum 
;       
; CALLING:
;     	bispec= polfft_bispec_iso(imag, iso = 30) 
;
; INPUTS:
;	imag --- input image
;
; KEYWORD:
;       iso : theta of isocele triangle (in degrees)
;     	PolarFFT\_CutOffParam : cut-off parameter in the polar FFT approximation
;	(the larger it is the better approximation is)
;   	equi : estimate only the equilateral configuration
;   	iso : estimate the isocele configuration with an angle equal to iso
;       Name : Name output file
;
; OUTPUT KEYWORDS:
;	x ---  x-axis (frequencies)
;
; OUTPUTS:
;       bispec --- bispectrum
;
; EXTERNAL CALLS:
;       cea_polar_fft (C++ program) 
;
;-----------------------------------------------------------------------
function polfft_bispec_iso, imag, x=x, iso = iso, PolarFFT_CutOffParam=PolarFFT_CutOffParam, Name=Name

b = systime(1)
if not keyword_set(PolarFFT_CutOffParam) then PolarFFT_CutOffParam = 8

;run polar fft
print, 'THE TIME TO PROCESS AN EXACT POWER SPECTRUM IS MUCH MORE LONGER'
polar_fft, imag, real, ima, PolarFFT_CutOffParam=PolarFFT_CutOffParam, Name=Name

;set dimensions
sz =size(imag)
n1 = sz[1]
n2 = sz[2]
if (n1 ne n2) then print, "map must be square !"
sz =size(real)
R = sz[1]
T = sz[2]

;reform real and ima
Nr = (R/2.)
Nt = T*2
re_imag = fltarr(Nr, Nt)
im_imag = fltarr(Nr, Nt)
re_imag(*, 0: T-1) = real(n1*3./4.: n1*3./2-1, *) 
im_imag(*, 0: T-1) = ima(n1*3./4.: n1*3./2-1, *) 
re_imag(*, T: 2*T-1) = reverse(real(1:n1*3./4., *)) 
im_imag(*, T: 2*T-1) = reverse(ima(1:n1*3./4., *))

;bispectrum calculation
re_bispec = fltarr(Nr)
im_bispec = fltarr(Nr)
som = fltarr(Nr)
x = fltarr(Nr)

dtheta = 2.*!dpi/double(Nt)
iso_rad = iso*!pi/180.
t_iso = (iso*double(Nt)/360.)

for i = 1, Nr-1 do begin
  l3 = double(i) 
  l1 = round(sqrt((l3^2/2.)/(1.+cos(iso_rad))))
  l2 = l1
  if l3 lt l1 or l3 gt 2*l1 then continue ;triangle not defined
  t13 =round((atan((l2*sin(iso_rad))/(l1+l2*cos(iso_rad))))/dtheta)+Nt/double(2)
  for k = 0, Nt-1 do begin
    t13k = t13+k
    t12k = t_iso+k
    if t13k gt Nt-1 then t13k = t13k-Nt
    if t12k gt Nt-1 then t12k = t12k-Nt
    a1 = re_imag(l1, k)
    b1 = im_imag(l1, k)
    a2 = re_imag(l2, t12k)
    b2 = im_imag(l2, t12k)
    a3 = re_imag(l3, t13k)
    b3 = im_imag(l3, t13k)
    ;B(k1)B(k2)B(k3)* = B(k1)B(k2)B(-k3) = (a1+ib1)(a2+ib2)(a3+ib3)
    re_bispec(l3) = re_bispec(l3) +a1*a2*a3 -a1*b2*b3 -b1*a2*b3 -b1*b2*a3
    im_bispec(l3) = im_bispec(l3) -b1*b2*b3 +b1*a2*a3 +a1*b2*a3 +a1*a2*b3
    som(l3) = som(l3) + 1.
    x(l3) = l3
  endfor
endfor

bispec=complex(re_bispec, im_bispec)
ind = where(som ne 0, count)
if count ne 0 then bispec(ind) = bispec(ind)/(double(som(ind))*n1^2)
;x=indgen(Nr)/ float(R)
x=x/float(R)
print, 'time =', systime(1)-b
return, bispec
end
