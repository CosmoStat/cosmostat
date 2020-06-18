;+
; NAME:
;        NFIRSTWT1D
;
; PURPOSE:
;	Computes a multiresolution transform of a signal. Set to zero
;       the smallest coefficients and reconstruct the signal.
;       The number of coefficient we want to preserve from the thresholding
;       can be given either by the keyword NBR or using the keyword  NFPER.
;       NFPER is the percentage of coefficients we want to keep.
;       Default transform is the orthogonal wavelet transform using five scales,
;       and a L2 normalization. The default value is 10%. 
;
; CALLING:
;
;      NFIRSTWT1D, Imag, Rec, Opt=Opt, nbr=nbr, nfper=nfper
;       
; INPUTS:
;     Imag -- 1D IDL array: image we want transform
;    
;
; OUTPUTS:
;     Rec -- 1D: reconstructed image 
;
; KEYWORDS:
;      nbr -- int:  Number of coefficients we want to keep
;
;      nfper -- int: Percentage of coefficients we want to keep.
;                    Default is 10%.
;
;      Opt -- string: string which contains the differents options. Options are:
;
;        where options = 
;
;          [-t type_of_multiresolution_transform]
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
;               11: Morlet's wavelet transform 
;               12: mexican hat wavelet transform 
;               13: french hat wavelet transform 
;               14: Gaussian Derivative wavelet transform 
;               15: bi-orthogonam transform   
;               16: bi-orthogonam transform via lifting sheme (CDF filters) 
;               17: Wavelet packets (CDF filters) 
;               18: Wavelet packets from lifting sheme 
;               19: Wavelet packets using the a-trous algorithm 
;               Default is 3.
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
;          [-l type_of_lifting_transform]
;              1: Lifting scheme: CDF WT 
;              2: Lifting scheme: median prediction 
;              3: Lifting scheme: integer Haar WT 
;              4: Lifting scheme: integer CDF WT 
;              5: Lifting scheme: integer (4,2) interpolating transform 
;              6: Lifting scheme: 7/9 WT 
;              7: Lifting scheme: integer 7/9 WT 
;             default is Lifting scheme: integer Haar WT
;   
;           [-n number_of_scales]
;    	        number of scales used in the multiresolution transform 
;    	        Default value is automatically calculated.
;
;          [-l type_of_lifting_transform]
;              1: Lifting scheme: CDF WT 
;              2: Lifting scheme: median prediction 
;              3: Lifting scheme: integer Haar WT 
;              4: Lifting scheme: integer CDF WT 
;              5: Lifting scheme: integer (4,2) interpolating transform 
;              6: Lifting scheme: 7/9 WT 
;              7: Lifting scheme: integer 7/9 WT 
;             default is Lifting scheme: integer Haar WT
;
;           [-L]
;                Use a L2 normalization. Default is L1
;
; EXTERNAL CALLS:
;       mr1d_trans (C++ program)
;       mr_recons (C++ program)
;
; EXAMPLE:
;
;       Compute the multiresolution of the signal I with default options
;       (i.e. OWT with 5 scales, 10% of the coefficients are not thresholded). 
;               NFIRSTWT1D, I, Rec
;
;       Ditto, but keep only 5%
;              NFIRSTWT1D, I, Rec, nfper=5.
;
;       Ditto, but use an undecimated WT instead of the OWT
;              NFIRSTWT1D, I, Rec, nfper=5., OPT='-t 7 -n 5 -L'
;
;-

pro nfirstwt1d, ima, rec, opt=opt, nbr=nbr, nfper=nfper

if not keyword_set(opt) then opt='-t15 -L -n5 -v'

mr1d_trans, ima, trans, opt=opt
wt = trans.coef
n = long(n_elements(wt))
if not keyword_set(nbr) then nbr = -1 $
else nfper =  float(nbr)  / float(n) * 100.

if not keyword_set(nfper) then nfper = -1
if nbr lt 0 and nfper lt 0 then nfper = 10.

nf = long(n * float(nfper) / 100.)
print, " WT, OPT = ", opt, "Total nbr of coeff = ", n
print, " Keep percentage P = ", nfper, " ==> KEEP NFIRST COEF: NF = ", nf

t = fltarr(n)
t(*) = wt
t1 = abs(t)
a = sort(t1)

t(a(0:n-1-nf)) = 0
wt(*)= t(*)
trans.coef = wt
mr1d_recons, trans, rec

plot, rec
end
