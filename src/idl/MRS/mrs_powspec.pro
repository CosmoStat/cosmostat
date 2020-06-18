;+
; NAME:
;        mrs_powspec
;
; PURPOSE:
;   Computes the power spectrum of a map,  using the HEALPix representation (NESTED data
;   representation by default) or the GLESP representation. 
;   If the keyword log is set, it is the log-power spectrum which  is returned. 
;
;  If the global variable DEF_NORM_POWSPEC is equal to 1 or if the keyword /Normalisation is set, 
;  then a normalization is performed, so that a Gaussian randomn noise with variance equal to 1 
;  has a power spectrum equal to 1. 
;
; CALLING:
;     P = mrs_powspec( Imag, plot=plot, lplot=lplot, log=log, IndL=IndL, PowSpecIma=PowSpecIma, StdPS=StdPS, Normalisation=Normalisation, nonorm=nonorm, NormVal=NormVal, alm=alm, lmax=lmax )
;
; INPUTS:
;     Imag -- IDL array of healpix map or GLESP structure: Input image to be transformed 
;    
; OUTPUTS:
;     P -- 1D IDL fltarr: Power Spectrum. P[k] = Mean(  POWSPECTRUM[*,l]  )
;
; INPUT KEYWORDS:
;		plot  -- int: if set, the power spectrum is plotted.
;		lplot -- int: if set, the power spectrum multiplied by l(l+1) is plotted.
;		log   -- int: if set, the log power spectrum is calculated instead of the power  spectrum
;		Lmax  -- int: Number of spherical harmonics computed in the decomposition and size of the computed spectrum (Lmax+1)
;					  (HEALPIX==> default is 3*nside, should be between 2*nside and 4*nside)
;					  (GLESP==> default is: min([Imag.nx/2, Imag.np/4]) )
;      Normalisation: if set, a l2 normalization if perform, so a Gaussian randomn noise with variance equal to 1 will have a power spectrum equal to 1. 
;	   nonorm: If set, no normalisation is performed on the alm computed.
;
; OUTPUT KEYWORDS:
;		IndL -- 1D IDL array: array contains the l(l+1) values.
;		PowSpecIma -- 2D IDL fltarr: Power spectrum of the input data
;		NormVal -- float: Normalisation value applied on the ALM (unless keyword nonorm set).
;		StdPS -- 1D IDL float array: contains the estimated standard deviation on the power spectrum computed.
;		alm -- IDL structure: result of the alm transform of the input image (see mrs_almtrans) with options lmax, /tab, /psp, norm=?
;
; EXTERNAL CALLS:
;       mrs_almtrans.pro
;
; EXAMPLE:
;       Compute the power spectrum of an image. 
;               P = mrs_powspec( Imag ) 
;         
; HISTORY:
;	Written: Jean-Luc Starck, 2005
;	December, 2005 File creation
;--------------------------------------------------------------------------------------------------------
 
;=====================================================================
 
function mrs_powspec, Imag, plot=plot, lplot=lplot, log=log, IndL=IndL, PowSpecIma=PowSpecIma, StdPS=StdPS, nonorm=nonorm, Normalisation=Normalisation, NormVal=NormVal, alm=alm,lmax=lmax, median=median
COMMON MR1ENV

if N_PARAMS() LT 1  then begin 
        print, 'CALLING SEQUENCE: p = mrs_powspec(Imag, plot=plot, lplot=lplot, log=log, IndL=IndL, PowSpecIma=PowSpecIma, StdPS=StdPS, Normalisation=Normalisation, nonorm=nonorm, NormVal=NormVal, alm=alm, lmax=lmax)'
	Ret=-1
        goto, DONE
        end

L2PowSpecNorm=0	
if DEF_NORM_POWSPEC EQ 1 then L2PowSpecNorm=1
if keyword_set(nonorm) then L2PowSpecNorm=0 $
else if keyword_set(Normalisation) then L2PowSpecNorm=1
mrs_almtrans, Imag, ALM, lmax=lmax, ring=ring, /tab, /psp, norm=L2PowSpecNorm

NormVal = ALM.NormVal
NbrAlm = ALM.LMAX+1
TabM = ALM.TABNBRM
PowSpecIma = ALM.ALM[*,*,0]
StdPS = dblarr(NbrAlm)
Ret = dblarr(NbrAlm)
 
 
for i=0, NbrAlm-1 do  begin
   allalm = 2*TabM[i]-1
   if allalm GT 1 then begin
      P = dblarr(allalm)
      P[0:TabM[i]-2] = reverse( PowSpecIma[i, 1:TabM[i]-1])
      P[TabM[i]-1:*] = PowSpecIma[i, 0:TabM[i]-1]
      if not keyword_set(Median) then Ret[i] = mean( P) $
      else  Ret[i] = median( P ) ;  
      StdPS[i] = sigma(P)
   end else begin
       Ret[i] = PowSpecIma[i, 0]
       StdPS[i] = 0.
   end
end

if keyword_set(log) then Ret = alog(Ret)
if keyword_set(plot) then plot, Ret

IndL = dindgen(NbrAlm)
IndL1 = IndL + 1.
IndLp = IndL*IndL1
IndL=IndLp
if keyword_set(log) then Tit='Log Power Spectrum' else Tit='Power Spectrum'

if keyword_set(plot) then plot, Ret,  xtitle='l', ytitle='P(k)', title=Tit
if keyword_set(lplot) then plot, Ret*IndLp / (2.*!PI), ytitle='l(l+1) C(l) / 2PI', xtitle='l', title=Tit

DONE:
return, Ret
END

pro zz
; save, filename='sim.xdr', c, n50, cn50, s, sigmanoise
; restore, filename='sim.xdr'
zc = mrs_powspec(c)
zn50 = mrs_powspec(n50)
zcn50 = mrs_powspec(cn50)
L = findgen(769)
l1= l*(l+1) / (2.*!PI)
plot, xtitle='l', ytitle='Cl l(l+1)/PI', zc*l1, thick=2, line=2
oplot, (zcn50 - zn50)*l1
SigmaNoise = sigma(n50)
s1 = sigma(c)
zc1 = mrs_powspec(c/s1)
plot, zc*l1/s1^2., line=2
oplot, zc1*l1

cmblet_filter, cn50, F1,  Spec1d=S1D, SigmaNoise=SigmaNoise, DPS=DPS, FPS=FPS, trans=trans, Itrans=Itrans, pf=pf, WienerWindow=WienerWindow, pm=pm
pf1 =  mrs_powspec(f1)
mrs_wiener, cn50, SigmaNoise^2., WRec,  Spec1D=Sw1D, WienerFilter=WienerFilter
rw = c - WRec
print, sigma(rw)
print, sigma(c-wrec), sigma(c-f1)  
pf1 =  mrs_powspec(f1)
pw1 = mrs_powspec(WRec)
end
