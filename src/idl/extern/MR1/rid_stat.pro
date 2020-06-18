;+
; NAME: 
;       RID_STAT
;
; PURPOSE: 
;       Calculate statistical information about the curvelet transform
;       of an image.
;
; CALLING:
;       RID_Info, Imag, TabStat, OPT=Opt, nodel=nodel, NameRes=NameRes
;
; INPUT:
;       Imag: image 
;
; KEYWORDS:
;
;      NameRes: Outfile file name for the statistic table
;
;      nodel: if set, the created files (by the program rid_stat)
;             are not deleted.
;
;      Opt: string which contains the differents options. Options are:
;
;
;         [-n number_of_scales]
;             number of scales used in the multiresolution transform
;              default is 5
;
;         [-N number_of_scales]
;             Number of scales used in the ridgelet transform.
;             default is automatically calculated. 
;
;         [-b BlockSize]
;             Block Size.
;             default is 16. 
;
;         [-O]
;             Use overlapping block. Default is no. 
;
;	     
;         [-v] 
;              Verbose
;
; OUTPUTS:
;           TabStat: fltarr(NbrBand, NbrStat) 
;                    TabStat(j,0) = standard deviation of the jth band
;                    TabStat(j,1) = skewness of the jth band
;                    TabStat(j,2) = kurtosis of the jth band
;                    TabStat(j,3) = minimum of the jth band
;                    TabStat(j,4) = maximum of the jth band

; EXTERNAL CALLS:
;           rid_stat (C++ program)
;
; EXAMPLE:
;           calculate the statistics with all default options 
;                rid_stat, Imag, TabStat
;
; HISTORY:
;	Written: Jean-Luc Starck 2002.
;-

PRO rid_stat, imag, result, OPT=OPT, nodel=nodel, NameRes=NameRes

if N_PARAMS() LT 2 then begin 
        spawn, 'rid_stat'
        print, 'CALL SEQUENCE: rid_stat, Imag, Imag_Out, OPT=Opt, nodel=nodel'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

Nl = (size(imag))[2]
Nc = (size(imag))[1]

NameImag = 'xx_imag.fits'
if keyword_set(NameRes) then NameResult = NameRes $
else NameResult = 'xx_result.fits' 

writefits, NameImag, imag 
com = 'rid_stat ' + OPT + ' ' + NameImag + ' ' + NameResult
spawn, com

Result = readfits(NameResult, /silent)
if not keyword_set(nodel) then begin
delete, NameImag
if not keyword_set(NameRes) then delete, NameResult
end

DONE:
end
