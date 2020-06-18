;+
; NAME: 
;       CUR_STAT
;
; PURPOSE: 
;       Computes the standard deviation, the skewness and the curtosis
;       at each resolution level of the curvelet transform of the 2D input data.
;       
; CALLING:
;
;       CUR_STAT, Image, Result, Opt=Opt
;
; INPUTS:
;       Image -- IDL 2D array: Image to analyze
;
; KEYWORDS:
;
;      Opt: string which contains the differents options. Options are:
;
;           [-g sigma]
;               Gaussian noise
;                  sigma = noise standard deviation 
;               by default, the noise is gaussian, and the standard 
;               devaition is automatically estimated. 
;
;           [-n number_of_scales]
;               number of scales used in the multiresolution transform
;               default is 4.
;               default is 6 in case of poisson noise with few events 
;
;           [-N number_of_scales]
;               Number of scales used in the ridgelet transform.
;               default is automatically calculated. 
;
;           [-b BlockSize]
;               Block Size. Default is 16.
;
;           [-O]
;               Use overlapping block. Default is no. 
;
;           [-v]
;               Verbose. Default is no.
;
; OUTPUTS:
;           Result: IDL array = cotnains the results
;                     Result(i,0) = standard deviation of the band i
;                     Result(i,1) = skewness of the band i
;                     Result(i,2) = curtosis of the band i
;
; EXTERNAL CALLS:
;           cur_stat (C++ program)
;
; EXAMPLE:
;           caculate the statistics of an image using 5 scales in the
;           the curvelet transform. 
;                cur_STAT, Image, Result, opt='-n5'
;
;
; HISTORY:
;	Written: Jean-Luc Starck 2002.
;	January, 2002 File creation
;-

PRO cur_stat, image, result, OPT=OPT, NameRes=NameRes, nodel=nodel

if N_PARAMS() LT 2 then begin 
        spawn, 'cur_stat'
        print, 'CALL SEQUENCE: cur_stat, image, result, OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

NameSig = 'xx_imag.fits'
if not keyword_set(NameRes) then NameRes = 'xx_result.fits' 
NameResult = NameRes 

writefits, NameSig, Image

com = 'cur_stat ' + OPT + ' ' + NameSig + ' ' + NameResult
spawn, com
Result = readfits(NameResult, /silent)

if not keyword_set(nodel) then begin
delete, NameSig
delete, NameResult
end

DONE:
end
