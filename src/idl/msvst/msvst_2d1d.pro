;+
; NAME: 
;       MSVST_2D1D
;
; PURPOSE: 
;       filter  a data cube (two spatial dimensions + 1 temporal or energy dimension) using the 2D1D wavelet transform and 
;       the multiscale variance stabilization. If a model is given, the routine estimates the difference between the data and the model.
;       
; CALLING:
;
;       msvst_2d1d, Data, Result, Model=Model, Opt=Opt
;
; INPUTS:
;       Data: Data to filter
;
; KEYWORDS:
;    Model -- IDL array: same size as input array. Contains a model of solution.
;
;      Opt: string which contains the differents options. Options are:
;    [-T] use g=Id-h*h as iteration band-pass filter
;    [-M value]
;         M = 0 : p-value threshold (default)
;         M = 1 : FDR threshold
;    [-E value] two sided Pr or FDR (default: Pr = 3.5*sigma or FDR = 0.1)
;    [-s value] Equivalent Gauss-detection level (default = 3.5*sigma)
;    [-c value] manual FDR setting, 0/1 = independence/dependence of coef. (default = 1)
;    [-n value] number of scalexy (default = 3)
;    [-N value] number of scalez (default = 5)
;    [-F value] first detection scalexy (default = 1)
;    [-f value] first detection scalez (default = 1)
;    [-K] ignore the last approximation band (used with iteration)
;    [-p] detect only the positive coefficients (used with iteration)
;    [-I value] iteration modes for M = 0,1
;                      I = 0 : Direct iterative mode
;                      I = 1 : L1-regularized iterative mode (default)
;    [-i value] number of iteration (default = 10)
;    [-Q value] write SNR file for every band
;    [-v] verbose mode
;
; OUTPUTS:
;           Result: result of the Poisson denoising
;
; EXTERNAL CALLS:
;           msvst_2d1d (C++ program)
;
; EXAMPLE:
;           filter an image with all default options 
;                msvst_2d1d, Imag, Result
;
; HISTORY:
;	Written: Jean-Luc Starck 2013.
;	November, 2013 File creation
;-

PRO msvst_2d1d, imag, result, Model=Model, OPT=OPT
COMMON MR1ENV
 
if N_PARAMS() LT 2 then begin 
        spawn, 'msvst_2d1d'
        print, 'CALL SEQUENCE: msvst_2d1d, Data, Data_Out, Model=Model, OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

NameImag = gettmpfilename() 
NameResult = gettmpfilename() 
writefits, NameImag, imag
NameModel = gettmpfilename() 
if keyword_set(Model) then writefits, NameModel, Model

com = BIN_ISAPCXX + '/msvst_2d1d  -T ' + OPT + ' ' + NameImag + ' ' + NameResult
if keyword_set(Model) then  com = com +  ' ' + NameModel


; com = 'msvst_iwt2d_coupled ' + OPT + ' ' + NameImag + ' ' + NameResult

spawn, com
Result = readfits(NameResult, /silent)

delete, NameImag
delete, NameResult
if keyword_set(Model) then delete, NameModel

DONE:
end
