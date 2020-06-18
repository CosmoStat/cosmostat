PRO MW1D_FILTER, Signal, Result, OPT=OPT
;+ 
; NAME: 
;        MW1D_FILTER
;
; PURPOSE: 
;       Filter a 1D signal by the multiscale entropy.
;
; CALLING SEQUENCE:
; 
;       MW1D_FILTER, Signal, Result, Opt=Opt
;
; INPUTS: 
;   Signal -- 1D IDL array: signal to filter 
;    
; KEYED INPUTS: 
;      Opt: string which contains the differents options. Options are:
;
;          [-t type_of_multiresolution_transform]
;              1: linear wavelet transform: a trous algorithm 
;              2: b1spline wavelet transform: a trous algorithm 
;              3: b3spline wavelet transform: a trous algorithm 
;              4: Derivative of a b3spline: a trous algorithm 
;              5: undecimated Haar wavelet transform: a trous algorithm 
;              6: morphological median transform 
;              7: Undecimated (bi-) orthogonal wavelet transform
;              8: pyramidal linear wavelet transform 
;              9: pyramidal b3spline wavelet transform 
;              10: pyramidal median transform 
;              Default is b3spline wavelet transform: a trous algorithm 
;
;          [-m type_of_noise]
;              1: Gaussian noise 
;              2: Poisson noise 
;              3: Poisson noise + Gaussian noise 
;              4: Multiplicative noise 
;              5: Non-stationary additive noise 
;              6: Non-stationary multiplicative noise 
;              7: Undefined stationary noise 
;              8: Undefined noise 
;             Default is Gaussian noise.
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
;           [-s NSigma]
;               Thresolding at NSigma * SigmaNoise at each scale
;               default is 3
;
;           [-i number_of_iterations]
;               Maximum number of iterations
;               default is 10.
;
;           [-e epsilon]
;               Convergence parameter
;               default is 0.000100.
;
;           [-n number_of_scales]
;               number of scales used in the multiresolution transform.
;
;           [-G RegulParam]
;               Regularization parameter 
;               default is 1
;
;           [-D]
;               Alpha is modified using the data SNR.
;               default is no.
;
;           [-w FilterCoefFileName]
;               Write to the disk the filtered wavelet coefficient.
;
;           [-v]
;                Verbose. Default is no.
;
;
; OUTPUTS: 
;   Result -- 1D IDL array: filtered signal
;
;
; MODIFICATION HISTORY: 
;	Written: Jean-Luc Starck 1999.
;       October, 1999 File creation
;-
 
 ;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 2 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', $ 
    'MW1D_FILTER, Signal, Result, Opt=Opt'
   GOTO, DONE
 ENDIF
 
if not keyword_set(Opt) then Opt = ' '  

NameSignal = 'xx_signal.fits'
NameResult = 'xx_result.fits'

writefits, NameSignal, Signal

com = 'mw1d_filter ' + OPT + ' ' + NameSignal + ' ' + NameResult
spawn, com
Result = readfits(NameResult, /silent)

delete, NameSignal
delete, NameResult

 DONE:
  RETURN
 
 END
