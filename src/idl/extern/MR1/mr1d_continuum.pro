;+
; NAME:
;        MR1D_CONTINUUM
;
; PURPOSE:
;	Estimate the contimuum of a 1D signal.
;
;
; CALLING:
;
;   output=MR1D_CONTINUUM(Signal,Nscale,Sigma=Sigma, niter=niter
;                         median=median,mirror=mirror)
;
; INPUTS:
;     Signal -- 1D IDL array: Input spectrum
;    
;     Nscale -- scalar: number of scales
;
; OUTPUTS:
;     output -- 1D: continuum
;
; INPUT KEYWORDS:
;     niter -- scalar: number of iterations. Default is 5.
;     median -- scalar: if set use multiresolution median transform
;     mirror -- scalar: if set, use mirror borders.
;
; OUTPUT KEYWORDS:
;    Sigma -- float: noise estimation
;
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;-
 
function MR1D_CONTINUUM, Signal, Nscale, Sigma=Sigma, median=median, Niter=Niter, mirror=mirror

continuum=-1
if N_PARAMS() LT 2 then begin 
        print, 'CALLING SEQUENCE: MR1D_CONTINUUM, Signal, Nscale, Sigma=Sigma'
        goto, DONE
        end


TabSmooth = [1.,0.521788,0.352123,0.244073,0.167545, 0.116704, 0.0816195, 0.0816195, 0.0816195]
TabSmoothMed = [1.,0.531919,0.399170,0.311026,0.239450, 0.178777, 0.139950, 0.139950, 0.139950]
if keyword_set(median) then TabSmooth=TabSmoothMed

if keyword_set(median) then mr1d_pyrmed, Signal, Wavelet, Nscale, /interp $
ELSE  BEGIN
    if keyword_set(mirror) then mr1d_atrou, Signal, Wavelet, Nscale, /mirror $
    else mr1d_atrou, Signal, Wavelet, Nscale
END

if not keyword_set(niter) then niter=5

continuum = Wavelet[*,Nscale-1]

Ns=Nscale
for i=0, Niter-1 do BEGIN
  sig = Signal - continuum

;  if keyword_set(median) then mr1d_pyrmed, sig, Wavelet, Ns, /interp $
  if keyword_set(median) then mr1d_pavemed, sig, Wavelet, Ns  $
  ELSE BEGIN
    if keyword_set(mirror) then  mr1d_atrou, sig, Wavelet, Ns, /mirror $
    else mr1d_atrou, sig, Wavelet, Ns
END

  continuum = continuum + Wavelet[*,Ns-1]
END

Sigma = sigma_clip(Signal-continuum)*TabSmooth[Nscale-1]


DONE:
return, continuum
END
