;+
; NAME:
;        IM_INTERPOLATE
;
; PURPOSE:
;	Interpolatation of 2D data
;
; CALLING:
;
;      IM_INTERPOLATE, ImagIn, ListIn, ListOut 
;       
;
; INPUTS:
;     ImagIn -- 2D IDL array: image we want to interpolate
;    
;     ListIn -- 2D IDL array:: Position to interpolate 
;                     [ListIn(0,i), ListIn(1,i)] Coordinates of the ith position 
;
; OUTPUTS:
;     ListOut -- 1D: Interpolation result
;
; KEYWORDS:
;
; EXTERNAL CALLS:
;       im_interpolate (C++ program)
;
; EXAMPLE:
;
; HISTORY:
;	Written: Jean-Luc Starck 2005.
;	March, 2005 File creation
;-

pro IM_INTERPOLATE, ImagIn, List, Result 

ImagIn = float(ImagIn)
if N_PARAMS() LT 3 then begin 
        print, 'CALLING SEQUENCE: IM_INTERPOLATE, ImagIn, List, Result'
        goto, DONE
        end

vsize = size(ImagIn)
if vsize(0) NE 2 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: IM_INTERPOLATE, ImagIn, List, Result'
         goto, DONE
        end
 
 
NameImag='xx_imag.fits'
NameList='xx_list.fits'
NameRes = 'xx_temp.fits'

writefits, NameImag, ImagIn
writefits, NameList, List

com = 'im_interpolate' + ' ' +  NameImag + ' ' + NameList + ' ' + NameRes 
spawn, com
Result = readfits(NameRes)

delete, NameImag
delete, NameRes
delete, NameList

DONE:

END
