;+
; NAME: 
;				CST_G2EB
;
; PURPOSE: 
;				Calculate weak lensing E and B mode from a shear data set, assuming constraint on the borders.
;
; CALLING:
;				EB =  cst_g2eb(Gamma, psize=psize)
;
; INPUT: 
;				Gamma ---- IDL array [*,*,0:1] : shear data set
;
; OUTPUT:
;                 EB    ---- IDL array [*,*,0:1] :  E and B modes. E mode = EB[*,*,0] and B mode = EB[*,*,1] 
;
; INPUT KEYWORD:			
;               Niter --  int : number of iterations. Default is 20.
;               BZero  --  saclar : if set, then a B zero mode contraint is applied
;               BBorderZero --  saclar : if set, then a B zero mode contraint is applied on the border
;               sparse --  saclar : if set, then a wavelet sparsity is applied
;               soft --  saclar : : if set, a soft thresholding is applied instead of the hard thresholding. Valid only if sparse is set.
;               psize --- float : pixel size 
;               TrueSolEB -- IDL [*,*,0:1] = true solution in case of simu. Then print the reconstruction error at each iteration.
;
; HISTORY:
;				Written by J.-L. Starck, Oct 2013
;-
;-------------------------------------------------------------------------------

function cst_g2eb, GammaShear, BZero=BZero, BBorderZero=BBorderZero, sparse=sparse, soft=soft, Niter=Niter, TrueSolEB=TrueSolEB, psize = psize, nscale=nscale
; BZero : if set, then a B zero mode contraint is applied
; BBorderZero : if set, then a B zero mode contraint is applied on the border
; sparse: if set, then a wavelet sparsity is applied
; soft: if set, a soft thresholding is applied instead of the hard thresholding. Valid only if sparse is set.
; Niter: number of iterations
; TrueSolEB -- IDL [*,*,2] = true solution in case of simu. Then print the reconstruction error at each iteration.

;print,niter
;read,idum
if not keyword_set(Niter) then Niter=20
 if not keyword_set(nscale) then nscale=5

g1 = GammaShear
vs = size(g1)
N = vs[1]
N1 = 2*N
Dep = N/2
g = dblarr(N1,N1,2)
g [Dep: Dep + N -1 , Dep: Dep + N -1 , *, *] = g1
M = g * 0
M[Dep: Dep + N -1 , Dep: Dep + N -1 , *, *] = 1.
 
;Niter=20
 RecEB = g2eb(G,psize=psize)
gcf1 = star2d(RecEB(*,*,0), nscale=nscale)
gcf2 = star2d(RecEB(*,*,1), nscale=nscale)
MaxL = 0
Lambda = 0
if keyword_set(sparse) then begin
	dirac = g[*,*,0]
	dirac [*] = 0
	vs = size(dirac)
	Nx = vs[1]
	Ny = vs[2]
	dirac[Nx/2, Ny/2] = 1
	w  = star2d(dirac, Nscale=Nscale, Gen2=Gen2)
	TabNorm = dblarr(Nscale)
	for j=0, Nscale-1 do TabNorm[j] = sqrt( total(w[*,*,j]^2) )
	MaxL = max( [max(gcf1[*,*,0:nscale-2]), max(gcf2[*,*,0:nscale-2])])
    Lambda = MaxL
end   
StepL = Lambda / float(niter-1)
LastThreshold = 0
FirstThreshold = MaxL
DeltaThreshold = FirstThreshold - LastThreshold

RecG = g
RecEB = g2eb(RecG,psize=psize)

for i=0,Niter-1 do begin
    Lambda =  LastThreshold  + DeltaThreshold  * (1.-erf(2.8*double(i) / double(niter)))

   if keyword_set(BZero) then RecEB[*,*,1] = 0
   if keyword_set(BBorderZero) then RecEB[*,*,1] =  RecEB[*,*,1]  * M[*,*,1]
   if keyword_set(sparse) then wtsparsecstr, RecEB, lambda, soft=soft,TabNorm=TabNorm     ; l1cstr, RecEB, lambda   ; sparsecstr, RecEB, Lambda
   ; load, RecEB[*,*,1]
  ;  load, RecEB[*,*,0]
   RecG = eb2g(RecEB,psize=psize)
   Resi =  M*(G - RecG)
   RecG = RecG  + Resi
   RecEB = RecEB  + g2eb(Resi,psize=psize)
   if keyword_set(sparse) then print,  i+1, ", Lambda = ", Lambda, ", SigmaResi = ", sigma(Resi) $
   else print,  i+1, ", SigmaResi = ", sigma(Resi) 
   if keyword_set(TrueSolEB) then  print, '          SigmaErrK = ',  sigma( RecEB[Dep: Dep + N -1 , Dep: Dep + N -1 ] - TrueSolEB[*,*,0] ) / sigma(TrueSolEB[*,*,0]) * 100. 
   Lambda = Lambda - StepL
end

RecEB = RecEB[Dep: Dep + N -1 , Dep: Dep + N -1, * ] 
return, RecEB
end



