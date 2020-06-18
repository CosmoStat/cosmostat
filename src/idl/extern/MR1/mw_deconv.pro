;+
; NAME: 
;       MW_DECONV
;
; PURPOSE: 
;       Deconvolve of an image by the multiscale entropy.
;
; CALLING:
;
;       MW_DECONV, Imag, Psf, Result, Opt=Opt
;
; INPUTS:
;       Imag: image to deconvolve
;       Psf:  Point Spread Function
;
; KEYWORDS:
;      Opt: string which contains the differents options. Options are:
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
;                 11: morphological pyramidal minmax transform 
;                 12: pyramidal laplacian 
;                 13: decomposition on scaling function 
;                 14: Mallat's wavelet transform (7/9 filters) 
;                 15: Feauveau's wavelet transform 
;                 16: Feauveau's wavelet transform without undersampling 
;                 17: G transform (morphological min-max algorithm) 
;                 18: Haar's wavelet transform 
;                 19: half-pyramidal transform 
;                 20: mixed Half-pyramidal WT and Median method (WT-HPMT) 
;                 21: undecimated diadic wavelet transform (two bands per scale) 
;                 22: mixed WT and PMT method (WT-PMT) 
;                 23: undecimated Haar transform: a trous algorithm (one band per scale) 
;                 24: undecimated mallat transform (three bands per scale)
;                default is bspline wavelet transform: a trous algorithm
;
;           [-H]
;               1: H = Wavelet coef. Energy (N1-MSE).
;               2: H = Noise Information (N2-MSE, for Gaussian noise only).
;               Default is 1.
;               For Gaussian noise, default 2.
;
;           [-m type_of_noise]
;                1: Gaussian noise 
;                2: Poisson noise 
;                3: Poisson noise + Gaussian noise 
;                4: Multiplicative noise 
;                5: Non-stationary additive noise 
;                6: Non-stationary multiplicative noise 
;                7: Undefined stationary noise 
;                8: Undefined noise 
;                9: Stationary correlated noise 
;               default is Gaussian noise
;
;           [-g sigma]
;                Gaussian noise
;                  sigma = noise standard deviation 
;                by default, the noise is gaussian, and the standard 
;                deviation is automatically estimated. 
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
;           [-n number_of_scales]
;                number of scales used in the multiresolution transform
;                default is 4
;
;           [-s NSigma]
;                Thresolding at NSigma * SigmaNoise at each scale
;                default is 3
;
;           [-i number_of_iterations]
;              Maximum number of iterations
;              default is 500
;
;           [-e Epsilon]
;                Convergence parameter
;                default is 0.001
;
;           [-K]
;               Suppress the last scale. Default is no. 
;
;           [-R RMS_Map_File_Name]
;               RMS Map (only used with -m 5 and -m 9 options).
;
;           [-P]
;               Suppress the positivity constraint. 
;
;           [-f ICF_Fwhm]
;               Fwhm = full width at half maximum
;               only used if type_of_deconvolution in 3 or 5.
;
;           [-I ICF_FileName]
;               Fwhm = Full Width at Half Maximum.
;
;           [-F First_Guess]
;               Input solution file.
;
;           [-W DataWeightingType]
;               1: no weighting  
;               2: soft weighting 
;               3: hard weighting 
;               default is 3.
;
;           [-A RegulProtectType]
;               0: no regularization (all protected) 
;               1: no protection from regularization 
;               2: soft protection 
;               3: hard protection 
;               4: soft + hard protection 
;               default is 3.
;
;           [-G RegulParam]
;               Regularization parameter 
;                default is 1
;                with hard weighting (-W 3), default is 0.1
;	     
;           [-M NSigmaObj]
;                NSigma level for object multiresolution determination
;                default is 3
;                with hard weighting (-W 3), default is 1
;
;           [-C ConvergParam]
;                Convergence parameter. 
;                default is 1.
;
;           [-r residual_file_name]
;                if this option is set, the residual is written to 
;                the file of name residual_file_name. By default, the
;                residual is not written.
;
;           [-S]
;              do not shift automatically the maximum  
;              of the PSF the center.
;
;           [-v]
;              Verbose. Default is no.
; 
;
; OUTPUTS:
;           Result: result of the deconvolution
;
; EXTERNAL CALLS:
;           mw_deconv (C++ program)
;
; EXAMPLE:
;           deconvolve an image with all default options 
;                MW_DECONV, Imag, Psf, Result
;
;          same example, but impose the number of iterations to be 30
;                MW_DECONV, Imag, Psf, Result, OPT='-i 30 -e 0'
;
;
; HISTORY:
;	Written: Jean-Luc Starck 1999.
;	October, 1999 File creation
;-

PRO mw_deconv, imag, Psf, result, OPT=OPT

if N_PARAMS() LT 3 then begin 
        spawn, 'mw_deconv'
        print, 'CALL SEQUENCE: mw_deconv Imag Psf OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

Nl = (size(imag))(2)
Nc = (size(imag))(1)

NameImag = 'xx_imag.fits'
NamePsf = 'xx_psf.fits'
NameResult = 'xx_result.fits'

writefits, NameImag, imag
writefits, NamePsf, Psf

com = 'mw_deconv ' + OPT + ' ' + NameImag + ' ' + NamePsf + ' ' + NameResult
spawn, com
Result = readfits(NameResult, /silent)

delete, NameImag
delete, NameResult
delete, NamePsf
DONE:
end
