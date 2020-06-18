;+
; NAME: 
;       MW_STAT
;
; PURPOSE: 
;       Computes the standard deviation, the skewness,the curtosis and
;       the multiscale entropy at each resolution level of the curvelet 
;       transform of the 2D input data.
;       
; CALLING:
;
;       MW_STAT, Image, Result, Opt=Opt
;
; INPUTS:
;       Image -- IDL 2D array: Image to analyze
;
; KEYWORDS:
;
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
;                 11: pyramidal laplacian 
;                 12: morphological pyramidal minmax transform 
;                 13: decomposition on scaling function 
;                 14: Mallat's wavelet transform 
;                 15: Feauveau's wavelet transform 
;                 16: Feauveau's wavelet transform without undersampling 
;                 17: G transform (morphological min-max algorithm 
;                 18: Haar's wavelet transform 
;                 19: Half-pyramidal transform 
;                 20: Mixed Half-pyramidal WT and Median method (WT-HPMT) 
;                 21: diadic wavelet transform 
;                 22: Mixed WT and PMT method (WT-PMT) 
;                 23: undecimated Haar transform: a trous algorithm 
;                 24: undecimated mallat transform (three bands per scale)
;                 default is 2
;
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
;
;           [-n number_of_scales]
;                number of scales used in the multiresolution transform
;                default is 4.
;                default is 6 in case of poisson noise with few events 
;
;           [-N NiterSigmaClip]
;             iteration number used for local variance estimation.
;             default is 1 
;
;           [-R RMS_Map_File_Name]
;               RMS Map (only used with -m 5 and -m 9 options).
;
;           [-S SizeBlock]
;             Size of the  blocks used for local variance estimation.
;             default is 7 
;
;           [-v]
;             Verbose. Default is no.
;
; OUTPUTS:
;           Result: IDL array = cotnains the results
;                     Result(i,0) = standard deviation of the band i
;                     Result(i,1) = skewness of the band i
;                     Result(i,2) = curtosis of the band i
;                     Result(i,3) = multiscale entropy of the band i
;
; EXTERNAL CALLS:
;           mw_stat (C++ program)
;
; EXAMPLE:
;           caculate the statistics of an image 
;                MW_STAT, Image, Result
;
;
; HISTORY:
;	Written: Jean-Luc Starck 2002.
;	January, 2002 File creation
;-

PRO mw_stat, image, result, OPT=OPT, NameRes=NameRes, nodel=nodel

if N_PARAMS() LT 2 then begin 
        spawn, 'mw_stat'
        print, 'CALL SEQUENCE: mw_stat, image, result, OPT=Opt'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' '  

NameSig = 'xx_imag.fits'
if not keyword_set(NameRes) then NameRes = 'xx_result.fits' 
NameResult = NameRes 

writefits, NameSig, Image

com = 'mw_stat ' + OPT + ' ' + NameSig + ' ' + NameResult
spawn, com
Result = readfits(NameResult, /silent)

if not keyword_set(nodel) then begin
delete, NameSig
delete, NameResult
end

DONE:
end
