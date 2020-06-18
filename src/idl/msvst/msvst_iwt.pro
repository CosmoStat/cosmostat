;+
; NAME: 
;       MSVST_IWT
;
; PURPOSE: 
;       filter  a signal or an image contaminated by a Poisson noise by using the MSVST method and the starlet transform (i.e. isotropic wavelet transform).
;       
; CALLING:
;
;       msvst_iwt, Data, Result, , Opt=Opt
;
; INPUTS:
;       Data: Data to filter
;
; KEYWORDS:
;      Opt: string which contains the differents options. Options are:
;    [-T] use g=Id-h*h as iteration band-pass filter
;    [-M value]
;         M = 0 : p-value threshold (default)
;         M = 1 : FDR threshold
;    [-E value] two sided Pr or FDR (default: Pr = 3.5*sigma or FDR = 0.1)
;    [-s value] Equivalent Gauss-detection level (default = 3.5*sigma)
;    [-c value] manual FDR setting, 0/1 = independence/dependence of coef. (default = 1)
;    [-n value] number of scales (default = 5)
;    [-F value] first detection scale (default = 1)
;    [-K] ignore the last approximation band (used with iteration)
;    [-p] detect only the positive coefficients (used with iteration)
;    [-I value] iteration modes for M = 0,1
;                      I = 0 : Direct iterative mode
;                      I = 1 : L1-regularized iterative mode (default)
;    [-i value] number of iteration (default = 10)
;    [-v] verbose mode
;
; OUTPUTS:
;           Result: result of the Poisson denoising
;
; EXTERNAL CALLS:
;           msvst_iwt2d_coupled (C++ program)
;
; EXAMPLE:
;           filter an image with all default options 
;                msvst_iwt, Imag, Result
;
; HISTORY:
;	Written: Jean-Luc Starck 2013.
;	November, 2013 File creation
;-

PRO msvst_iwt, imag, result, OPT=OPT
COMMON MR1ENV
 
if N_PARAMS() LT 2 then begin 
        spawn, 'msvst_iwt2d_coupled'
        print, 'CALL SEQUENCE: msvst_iwt, Imag, Imag_Out, OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

NameImag = gettmpfilename() 
NameResult = gettmpfilename() 
writefits, NameImag, imag

com = BIN_ISAPCXX + '/msvst_iwt2d_coupled  -T ' + OPT + ' ' + NameImag + ' ' + NameResult
; com = 'msvst_iwt2d_coupled ' + OPT + ' ' + NameImag + ' ' + NameResult

spawn, com
Result = readfits(NameResult, /silent)

delete, NameImag
delete, NameResult
DONE:
end
