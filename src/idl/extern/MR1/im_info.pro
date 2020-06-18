;+
; NAME:
;        IM_INFO
;
; PURPOSE:
;	Compute statictics of an image.
;
; CALLING:
;
;      IM_INFO, ImagIn, Opt=Opt
;       
;
; INPUTS:
;     Imag -- 2D IDL array: image we want transform
;    
;
; OUTPUTS:
;
; KEYWORDS:
;      Opt -- string: string which contains the differents options. Options are:
;
;        where options = 
;          [-e]
;             Entropy and Fisher Information calculation.
;             Histogram bin size value is 1. .
;
;          [-E HistoBinSize]
;             Entropy and Fisher Information calculation.. 
;
;          [-M]
;            Skewness, Kurtosis, HC and HC^+  calculation.  
;
;
; EXTERNAL CALLS:
;       im_info (C++ program)
;
; EXAMPLE:
;
;       Compute  information about an image
;              IM_INFO, Data, OPT='-M'
;
; HISTORY:
;	Written: Jean-Luc Starck
;	December, 2005 File creation
;-

pro IM_INFO,  ImagIn, OPT=Opt, Survival=Survival
ImagIn = float(ImagIn)
if N_PARAMS() LT 1 then begin 
        print, 'CALLING SEQUENCE: IM_INFO, ImagIn, OPT=Opt'
        goto, DONE
        end

vsize = size(ImagIn)
if vsize(0) NE 2 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: IM_INFO, ImagIn, OPT=Opt'
         goto, DONE
        end
 
if not keyword_set(Opt) then Opt = ' '
filename = 'xx_temp.fits'
p = strpos(filename, '.fits')
if p LT 0 then filename = filename + '.fits'

NameImag='xx_imag.fits'
writefits, NameImag, ImagIn

com =  'im_info -S ' + Opt + ' ' + NameImag 
spawn, com
Survival = readfits('xx_survival.fits')

delete, NameImag
delete, 'xx_survival.fits'
DONE:

END
