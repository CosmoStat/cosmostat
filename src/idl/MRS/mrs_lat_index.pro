;+
; NAME:
;        mrs_lat_index
;
; PURPOSE:
;	 Return the pixel indices of all pixels in a latitude band. The default number of band is 25.
;
; CALLING:
;      LatInd = mrs_lat_index(Nside, Lat, NbrLat=NbrLat, Mask=Mask,  TabBLat=TabBLat)
;       
; INPUTS:
;		Nside -- int:  Healpix nside number 
;		Lat -- scalar: Latitude band number. 
;
; INPUT KEYWORDS:
;		NbrLat -- int: Number  of latitude bands.   
;
; OUTPUT KEYWORDS: 
;     Mask -- Healpix image:   Mask[LatInd] 1 and other pixels are set to zero.
;	  TabBLat -- 1D float array [0: NbrLat-1]: B latitude for each latitude band
;    
; EXAMPLE:
;     Compute the standard deviation in each latitude band:
;     for l=0,24 do begin  LatInd = mrs_lat_index(Nside, l,  TabBLat=TabBLat) & print, 'Lat ', TabBLat[l], ', Sigma = ', sigma( Ima[LatInd] ) & end
;
; HISTORY:
;	Written: Jean-Luc Starck & Jerome Bobin, April 2011
;-


function mrs_lat_index, Nside, Lat, NbrLat=NbrLat, Mask=Mask, MinLat=MinLat, MaxLat=MaxLat, MeanLat=MeanLat, freqband= freqband, theta=theta, phi=phi, TabBLat=TabBLat

if N_PARAMS() LT 2   then begin 
        print, 'CALLING SEQUENCE:  LatInd = mrs_lat_index(Nside, Lat, NbrLat=NbrLat, Mask=Mask,  TabBLat=TabBLat)'
        ind=-1
        goto, DONE
        end

if not keyword_set(NbrLat) then NbrLat = 25
NBand = NbrLat
LBand = 180./NBand

freqband = 0.*dblarr(2,Nband)
freqband[0,*] = indgen(Nband)*LBand  ;--- Lowest latitude
freqband[1,*] = indgen(Nband)*LBand + LBand  ;--- Highest latitude
freq_mean = 0.5*total(freqband,1)
ang2lb, freq_mean / 180. * !PI,  freq_mean*0, l, TabBLat

 resband = 0.*dblarr(NBand)
 
ipix = linspace(0,double(nside)^2.*12.-1,double(nside)^2.*12.)
ipix = long(ipix)
; help, ipix, nside
pix2ang_nest, nside, ipix,theta,phi
ang2lb, theta, phi, l, b

ll = Lat
 rad_max = freqband[1,ll]/180.*!PI
 rad_min = freqband[0,ll]/180.*!PI
MinLat =  freqband[0,ll]
MaxLat = freqband[1,ll]
MeanLat = freq_mean[ll]
Ns = long(Nside)
 mask = intarr( Ns^2*12) + 1 
        
 ind_max = where(theta gt rad_max, count)
 if count gt 0 then mask[ind_max] = 0.
ind_min = where(theta lt rad_min,count)
if count gt 0  then mask[ind_min] = 0.
               
ind = where(mask ne 0, count)
if count gt 0 then theta= theta[ind] else theta = 0
if count gt 0 then phi = phi[ind] else phi = 0

DONE:
return, ind
end

;==================================================================================
 
