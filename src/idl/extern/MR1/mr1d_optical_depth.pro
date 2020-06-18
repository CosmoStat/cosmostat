;+
; NAME:
;        MR1D_OPTICAL_DEPTH
;
; PURPOSE:
;	Estimate the optical depth of a spectrum.
;
; CALLING:
;
;      output = MR1D_OPTICAL_DEPTH(Signal, Nscale, Sigma=Sigma)
;
; INPUTS:
;     Signal -- 1D IDL array: Signal we want the continuum
;    
;     Nscale -- scalar: number of scales
;
; OUTPUTS:
;     output -- 1D: continuum
;
; OUTPUT KEYWORDS:
;    Sigma -- float: noise estimation
;
; HISTORY:
;	Written: Jean-Luc Starck 1995.
;	December, 1995 File creation
;-
 
function MR1D_OPTICAL_DEPTH, Signal, Nscale, Sigma=Sigma

tau=-1
if N_PARAMS() LT 2 then begin 
        print, 'CALLING SEQUENCE: MR1D_OPTICAL_DEPTH, Signal, Nscale, Sigma=Sigma'
        goto, DONE
        end

Log_Signal = alog(Signal)
continuum = mr1d_continuum(Log_Signal, Nscale, sigma=sigma)

tau = continuum - Log_Signal
DONE:
return, tau
END
