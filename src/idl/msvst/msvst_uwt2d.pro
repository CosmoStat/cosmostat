;+
; NAME: 
;       MSVST_UWT2D
;
; PURPOSE: 
;       filter  an image contaminated by a Poisson noise by using the MSVST method and the bi-orthogonal undecimated wavelet  transform (i.e. 7/9 filter bank).
;       
; CALLING:
;
;       msvst_uwt2d, Imag, Result, , Opt=Opt
;
; INPUTS:
;       Imag: image to filter
;
; KEYWORDS:
;      Opt: string which contains the differents options. Options are:
;    [-M value]
;         M = 0 : MS-VST (default)
;         M = 1 : Anscombe
;         M = 2 : Haar-Fisz
;    [-m] Use FDR thresholding
;    [-t value] max. cycle translation for Haar-Fisz (default = 0)
;    [-E value] two sided Pr or FDR (default: Pr = 3.5*sigma or FDR = 0.1)
;    [-s value] Equivalent Gauss-detection level (default = 3.5*sigma)
;    [-c value] manual FDR setting, 0/1 = independence/dependence of coef. (default = 1)
;    [-n value] number of scales (default = 5)
;    [-F value] first detection scale (default = 1)
;    [-K] ignore the last approximation band (used with iteration)
;    [-p] detect only the positive coefficients (used with iteration)
;    [-i value] number of iteration (default = 10)
;    [-v] verbose mode
;
; OUTPUTS:
;           Result: result of the Poisson denoising
;
; EXTERNAL CALLS:
;           msvst_uwt2d (C++ program)
;
; EXAMPLE:
;           filter an image with all default options 
;                msvst_uwt2d, Imag, Result
;
; HISTORY:
;	Written: Jean-Luc Starck 2013.
;	November, 2013 File creation
;-

PRO msvst_uwt2d, imag, result, OPT=OPT
COMMON MR1ENV
 
if N_PARAMS() LT 2 then begin 
        spawn, 'msvst_uwt2d'
        print, 'CALL SEQUENCE: msvst_uwt2d, Imag, Imag_Out OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

NameImag = gettmpfilename() 
NameResult = gettmpfilename() 
writefits, NameImag, imag

com = BIN_ISAPCXX + '/msvst_uwt2d ' + OPT + ' ' + NameImag + ' ' + NameResult
spawn, com
Result = readfits(NameResult, /silent)

delete, NameImag
delete, NameResult
DONE:
end
