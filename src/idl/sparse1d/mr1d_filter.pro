;+
; NAME: 
;       MR1D_FILTER
;
; PURPOSE: 
;       Filter a 1D signal by the multiresolution support method. If the input is a 2D array Data[*,0:N], then 
;       apply independantly N times the signal Data[*,i].
;       
; CALLING:
;
;       MR1D_FILTER, Signal, Result, Opt
;
; INPUTS:
;       Signal -- IDL 1D array: Signal to filter
;
; KEYWORDS:
;
;      Opt: string which contains the differents options. Options are:
;
;         [-f type_of_filtering]
;              1: Multiresolution Hard K-Sigma Thresholding 
;              2: Multiresolution Soft K-Sigma Thresholding 
;              3: Iterative Multiresolution Thresholding 
;              4: Universal Hard Thesholding 
;              5: Universal Soft Thesholding 
;              6: SURE Hard Thesholding 
;              7: SURE Soft Thesholding 
;              8: MULTI-SURE Hard Thesholding 
;              9: MULTI-SURE Soft Thesholding 
;              10: Median Absolute Deviation (MAD) Hard Thesholding 
;              11: Median Absolute Deviation (MAD) Soft Thesholding 
;              default is Multiresolution Hard K-Sigma Thresholding.
;
;         [-t type_of_multiresolution_transform]
;               1: linear wavelet transform: a trous algorithm 
;               2: b1spline wavelet transform: a trous algorithm 
;               3: b3spline wavelet transform: a trous algorithm 
;               4: Derivative of a b3spline: a trous algorithm 
;               5: undecimated Haar wavelet transform: a trous algorithm 
;               6: morphological median transform 
;               7: Undecimated (bi-) orthogonal wavelet transform
;               8: pyramidal linear wavelet transform 
;               9: pyramidal b3spline wavelet transform 
;               10: pyramidal median transform 
;
;          [-T type_of_filters]
;              1: Antonini 7/9 filters 
;              2: Daubechies filter 4 
;              3: Biorthogonal 2/6 Haar filters 
;              4: Biorthogonal 2/10 Haar filters 
;              5: Odegard 7/9 filters 
;              6: User's filters 
;              default is Antonini 7/9 filters
;
;    	    [-m type_of_noise]
;              1: Gaussian Noise 
;              2: Poisson Noise 
;              3: Poisson Noise + Gaussian Noise 
;              4: Multiplicative Noise 
;              5: Non uniform additive noise 
;              6: Non uniform multiplicative noise 
;              7: Undefined uniform Noise 
;              8: Undefined Noise 
;              9: Stationary correlated noise 
;              10: Poisson noise with few events 
;      	       Default is Gaussian noise.
;
;           [-g sigma]
;                Gaussian noise
;                  sigma = noise standard deviation 
;                by default, the noise is gaussian, and the standard 
;                devaition is automatically estimated. 
;
;           [-c gain,sigma,mean]
;                case of a CCD: noise = Poisson noise + read-out noise
;                  gain = CCD gain 
;                  sigma = standard deviation of the read-out noise
;                  mean = mean of the read-out noise
;                if this option is set, 
;                           Noise = Poisson + Gaussian read-out Noise
;                it is generally the case with the CCD.
;                Attention, these parameters must be separated by a comma 
;                without space. example: -c 0.133,1.733,0.
;                If mean, or sigma and mean are omitted, default values are 0.
;                gain can not be omitted. 
;
;           [-E Epsilon]
;                Epsilon = precision for computing thresholds
;                          (only used in case of poisson noise with 
;                           few events). Default is 1.00e-03 
;
;           [-n number_of_scales]
;                number of scales used in the multiresolution transform
;                default is 4.
;                default is 6 in case of poisson noise with few events 
;
;           [-s NSigma]
;                Thresolding at NSigma * SigmaNoise at each scale
;                default is 3
;
;           [-i number_of_iterations]
;              Maximum number of iterations
;              default is 10
;
;           [-e Epsilon]
;                Convergence parameter
;                default is 0.0001
;
;           [-K]
;             Suppress the last scale. Default is no.
;
;           [-p]
;             Detect only positive structure. Default is no.
;
;           [-k]
;             Suppress isolated pixels in the support. Default is no.
;
;           [-v]
;             Verbose. Default is no.
;
; OUTPUTS:
;           Result: result of the filtering
;
; EXTERNAL CALLS:
;           mr1d_filter (C++ program)
;
; EXAMPLE:
;           filter a signal with all default options 
;                MR1D_FILTER, Signal, Result
;
;          same example, but impose the number of scales  to be 3, and 
;          a thresolding at 5 sigma
;                MR1D_FILTER, Signal, Result, OPT='-n 3 -s 5'
;
;
; HISTORY:
;	Written: Jean-Luc Starck 1994.
;	December, 1995 File creation
;       January, 1998 more noise options
;-


PRO mr1d_filter, signal, result, OPT=OPT, rmsnoise=rmsnoise 

COMMON MR1ENV
COMMON C_PLANCK

if N_PARAMS() LT 2 then begin 
        spawn, 'mr_filter'
        print, 'CALL SEQUENCE: mr1d_filter Signal Signal_Out OPT=Opt'
        goto, DONE
        end

vs = size(signal)
Ndim = vs[0]
if vs[0] EQ 2 then OptFil = ' -M ' else OptFil = ' '

if keyword_set(Opt) then OptFil = OptFil + Opt   

NameSig = gettmpfilename()   
NameResult = gettmpfilename()  

writefits, NameSig, signal

if keyword_set(rmsnoise) then begin
   NameRMS = gettmpfilename()  
   OptFil = OptFil + ' -m5 -R '+ NameRMS
   writefits, NameRMS, rmsnoise
end   

com = BIN_ISAPCXX + '/mr1d_filter ' + OptFil + ' ' + NameSig + ' ' + NameResult
; print, com
spawn, com
Result = readfits(NameResult, /silent)

delete, NameSig
delete, NameResult
if keyword_set(rmsnoise) then delete, NameRMS
DONE:
end
