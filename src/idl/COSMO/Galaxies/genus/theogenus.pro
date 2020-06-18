;+
; NAME:
;        THEOGENUS
;
; PURPOSE:
;	Computes the theoretical Genus curve
;
; CALLING:
;
;      THEOGENUS, Nu, G, sig=sig, NuIndex=NuIndex, NbrPix=NbrPix, MaxNu=MaxNu
;       
;
; INPUTS:
;    
; OUTPUTS:
;     Nu -- 1D array: X-axis coordinate of the genus
;
;     G -- 1D array: genus
;
; KEYWORDS:
;
;     NuIndex -- float: Power spectrum index of the Gaussian Random Field 
;                       Default is -1.
;
;     sig   -- float: Standard deviation of the Gaussian used for the convolution with
;           the Gaussian Random Field. Default is 1.
;  
;     NbrPix -- int: Number of pixels used in the genus calculation.
;           Default is 512.
;
;     MaxNu -- float: genus is calculated for nu values between [-MaxNu,MaxNu].
;           Default is 3.
;
; HISTORY:
;-

pro theogenus, nu, g, sig=sig, NuIndex=NuIndex, NbrPix=NbrPix, MaxNu=MaxNu

if not keyword_set(NbrPix) then NbrPix=512

if not keyword_set(sig) then sig=2.
if not keyword_set(NuIndex) then NuIndex=-1.
if not keyword_set(MaxNu) then MaxNu=3.

N = NbrPix
nu = 2.*MaxNu*(findgen(N) / float(N-1) - 0.5)
Cst = 1./sqrt(2*!PI)
L=1./Sig*sqrt((NuIndex+3.)/6.)
Vlambda = L /sqrt(2*!PI)
V2 = Vlambda*Vlambda
Nu2 = nu^2
g = (Vlambda*V2*Cst*(Nu2-1)*exp(-Nu2/2))

end
