;+
; NAME: 
;       MW_FILTER
;
; PURPOSE: 
;       filter  an image by the multiscale entropy.
;       
; CALLING:
;
;       MW_FILTER, Imag, Result, Opt=Opt
;
; INPUTS:
;       Imag: image to filter
;
; KEYWORDS:
;      Opt: string which contains the differents options. Options are:
;
;
;          [-t type_of_multiresolution_transform]
;                  1: linear wavelet transform: a trous algorithm 
;                  2: bspline wavelet transform: a trous algorithm 
;                  3: wavelet transform in Fourier space 
;                  4: morphological median transform 
;                  5: morphological minmax transform 
;                  6: pyramidal linear wavelet transform 
;                  7: pyramidal bspline wavelet transform 
;                  8: pyramidal wavelet transform in Fourier space: 
;                     wavelet =  between two resolutions 
;                  9: pyramidal wavelet transform in Fourier space: 
;                     wavelet = difference  between the square 
;                                                of two resolutions
;                 10: pyramidal median transform 
;                 11: pyramidal laplacian 
;                 12: morphological pyramidal minmax transform 
;                 13: decomposition on scaling function 
;                 14: Mallat's wavelet transform 
;                 15: Feauveau's wavelet transform 
;                 16: Feauveau's wavelet transform without undersampling 
;                 17: G transform (morphological min-max algorithm 
;                 18: Haar's wavelet transform 
;                 19: half-pyramidal transform 
;                 20: mixed Half-pyramidal WT and Median method (WT-HPMT) 
;                 21: undecimated diadic wavelet transform (two bands per scale) 
;                 22: mixed WT and PMT method (WT-PMT)
;                 23: undecimated Haar transform: a trous algorithm (one band per scale) 
;                 24: undecimated mallat transform (three bands per scale)
;                 default is bspline wavelet transform: a trous algorithm
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
;           [-m]
;              1: Gaussian Noise 
;              2: Poisson Noise 
;              3: Poisson Noise + Gaussian Noise 
;              4: Multiplicative Noise 
;              5: Non uniform additive noise 
;              6: Non uniform multiplicative noise 
;              7: Undefined uniform Noise 
;              8: Undefined Noise 
;              9: Stationary correlated noise 
;             default is Gaussian noise
;
;           [-n number_of_scales]
;               number of scales used in the multiresolution transform
;               default is 4.
;               default is 6 in case of poisson noise with few events 
;
;           [-s NSigma]
;               Thresolding at NSigma * SigmaNoise at each scale
;               default is 3
;
;          [-S SizeBlock]
;             Size of the  blocks used for local variance estimation.
;             default is 7 
;
;           [-N NiterSigmaClip]
;               iteration number used for local variance estimation.
;               default is 1 
;
;           [-R RMS_Map_File_Name]
;               RMS Map (only used with -m 5 and -m 9 options). 
;
;           [-G RegulParam]
;               Regularization parameter 
;               default is 1
;
;           [-C CvgParam]
;               Convergence parameter. 
;               default is 0.01.
;
;           [-T Type_of_Regularization]
;               1: Use a fixed user Alpha value.
;               2: Estimate the optimal Alpha.
;               3: Estimate one  Alpha value per band.
;               default is 1.
;
;           [-D]
;               Alpha is modified using the data SNR.
;              default is no.
;
;           [-i MaxIter]
;              Maximum number of iterations
;              default is 20.
;
;           [-P]
;              Apply the positivity constraint.
; 
;           [-v]
;             Verbose. Default is no.
;
; OUTPUTS:
;           Result: result of the deconvolution
;
; EXTERNAL CALLS:
;           mw_filter (C++ program)
;
; EXAMPLE:
;           filter an image with all default options 
;                MW_FILTER, Imag, Result
;
;
; HISTORY:
;	Written: Jean-Luc Starck 1999.
;       October, 1999 File creation
;-

PRO mw_filter, imag, result, OPT=OPT

if N_PARAMS() LT 2 then begin 
        spawn, 'mw_filter'
        print, 'CALL SEQUENCE: mw_filter Imag Imag_Out OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

Nl = (size(imag))(2)
Nc = (size(imag))(1)

NameImag = 'xx_imag.fits'
NameResult = 'xx_result.fits'


writefits, NameImag, imag

com = 'mw_filter ' + OPT + ' ' + NameImag + ' ' + NameResult
spawn, com
Result = readfits(NameResult, /silent)

delete, NameImag
delete, NameResult
DONE:
end
