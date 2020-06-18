;+
; NAME: 
;       CUR_FILTER
;
; PURPOSE: 
;       filter  an image by using a multiresolution transform.
;       
; CALLING:
;
;       CUR_FILTER, Imag, Result, Opt
;
; INPUTS:
;       Imag: image to filter
;
; KEYWORDS:
;      Opt: string which contains the differents options. Options are:
;
; OUTPUTS:
;           Result: result of the filtering
;
; EXTERNAL CALLS:
;           cur_filter (C++ program)
;
; EXAMPLE:
;           filter an image with all default options 
;                CUR_FILTER, Imag, Result
;
;          same example, but impose the number of scales  to be 3, and 
;          a thresolding at 5 sigma
;                CUR_FILTER, Imag,  Result, OPT='-n 3 -s 5'
;-

PRO CUR_FILTER, imag, result, OPT=OPT

if N_PARAMS() LT 2 then begin 
        spawn, 'CUR_FILTER'
        print, 'CALL SEQUENCE: CUR_FILTER Imag Imag_Out OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

Nl = (size(imag))[2]
Nc = (size(imag))[1]

NameImag = 'xx_imag.fits'
NameResult = 'xx_result.fits'

writefits, NameImag, imag

com = 'cur_filter ' + OPT + ' ' + NameImag + ' ' + NameResult
spawn, com
Result = readfits(NameResult, /silent)

delete, NameImag
delete, NameResult
DONE:
end
