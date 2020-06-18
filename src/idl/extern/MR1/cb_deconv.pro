;+
; NAME: 
;       CB_DECONV
;
; PURPOSE: 
;       Deconvolve of an image by using wavelet and curvelet
;
; CALLING:
;
;       CB_DECONV, Imag, Psf, Result, Opt=Opt
;
; INPUTS:
;       Imag: image to deconvolve
;       Psf:  Point Spread Function
;
; KEYWORDS:
;       See binary cb_deconv 
;
; OUTPUTS:
;           Result: result of the deconvolution
;
; EXTERNAL CALLS:
;           cb_deconv (C++ program)
;
; EXAMPLE:
;           deconvolve an image with all default options  
;                CB_DECONV, Imag, Psf, Result
;
;          same example, but impose the number of iterations to be 30
;                CB_DECONV, Imag, Psf, Result, OPT='-i 30 -e 0'
;
; HISTORY:
;	Written: Jean-Luc Starck 2005.
;-

PRO cb_deconv, imag, Psf, result, OPT=OPT

if N_PARAMS() LT 3 then begin 
        spawn, 'mr_deconv'
        print, 'CALL SEQUENCE: cb_deconv Imag Psf OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

Nl = (size(imag))[2]
Nc = (size(imag))[1]

NameImag = 'xx_imag.fits'
NamePsf = 'xx_psf.fits'
NameResult = 'xx_result.fits'

writefits, NameImag, imag
writefits, NamePsf, Psf

com = 'cb_deconv ' + OPT + ' ' + NameImag + ' ' + NamePsf + ' ' + NameResult
spawn, com
Result = readfits(NameResult, /silent)

delete, NameImag
delete, NameResult
delete, NamePsf
DONE:
end
