; NAME			
;    POWER_TO_CORREL
;
; PURPOSE
;    find the link between correl_2pts and my_isospec
;voici ce que j'ai retrouvé dans mes archives  :
; un code en fortran :
;   http://casa.colorado.edu/~ajsh/FFTLog/#download
;  un code IDL en pièce jointe.
; Je ne sais plus lequel j'ai utilisé....
; Une référence en pièce jointe.
;  Szapudi, et al, FAST EDGE-CORRECTED MEASUREMENT OF THE TWO-POINT 
;     CORRELATION FUNCTION AND THE POWER SPECTRUM
;  The Astrophysical Journal, 631:L1–L4, 2005 September 20
;
; Sandrine
;
; Input
;
;
; Output
;    correl --- power spectrum
;
; History
; 	  S. Pires Feb. 2007
;;-------------------------------------------------------------------------
pro power_to_correl, theta, l, Ptot, scale=scale

NN = 128
k=readfits('/Users/sandrinepires/UNIX/Image/a_kapint.fits')
kk= k(0:NN-1,0:NN-1)
;kk=congrid(k, NN, NN)


my_isospec, kk, powspec, x=x, /exact
l = x*2.*!pi
sz = size(l)
N=sz[1]
print, 'N=', N
Ptot = fltarr(N)

n_theta=(3.*NN)/4 ;number of regular bins in theta
print, 'n_theta=', n_theta
theta=fltarr(n_theta+1)
pas = (4.*!pi)/(3.*NN)
print, 'pas=', pas
for i = 0, n_theta-1 do begin
 theta(n_theta-(i+1)) = 1/(pas*(i+1))
endfor
theta(n_theta)= NN/2. 
n_radius = (3.*float(NN))/4.
;n_radius = (3.*float(NN))/2. ;number of radius in power spectrum when r = [0, Rmax] in integration
for j = 0, N-1 do begin
  for i = 0, N-1 do begin
    ;Ptot(j) = Ptot(j) + (1./n_radius) * powspec(i)*(beselj(l(i)*theta(j), 0))*l(i)
    Ptot(j) = Ptot(j) + (1/(2.*n_radius)) * powspec(i)*(beselj(l(i)*theta(j), 0))*l(i)
  endfor
endfor


;plot, theta*(120./128.),  Ptot, color=20
oplot, theta,  Ptot, color=20

end
