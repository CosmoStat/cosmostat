;+
; NAME: 
;				CST_EB2G
;
; PURPOSE: 
;				Calculate weak lensing shear data  from E and B mode., assuming constraints on the border or on the B mode..
;
; CALLING:
;				G =  cst_eb2g(EBmode, psize=psize)
;
; INPUT: 
;				EB ---- IDL array [*,*,0:1] :  E and B modes. E mode = EB[*,*,0] and B mode = EB[*,*,1] 
;
; OUTPUT:
;               G  ---- IDL array [*,*,0:1] :   shear. gamma 1  = G[*,*,0] and gamma 2 = G[*,*,1] 
;
; INPUT KEYWORD:	
;               Niter --  int : number of iterations. Default is 20.
;               BZero  --  saclar : if set, then a B zero mode contraint is applied
;               BBorderZero --  saclar : if set, then a B zero mode contraint is applied on the border
;               sparse --  saclar : if set, then a wavelet sparsity is applied
;               soft --  saclar : : if set, a soft thresholding is applied instead of the hard thresholding. Valid only if sparse is set.
;               psize --- float : pixel size 
;               TrueSolGamma -- IDL [*,*,0:1] = true solution in case of simu. Then print the reconstruction error at each iteration.
;
; HISTORY:
;				Written by J.-L. Starck, Oct 2013
;-
;-------------------------------------------------------------------------------

function cst_eb2g, EBMode, BZero=BZero, BBorderZero=BBorderZero, sparse=sparse, soft=soft, Niter=Niter, TrueSolGamma=TrueSolGamma, psize = psize, nscale=nscale
; BZero : if set, then a B zero mode contraint is applied
; BBorderZero : if set, then a B zero mode contraint is applied on the border
; sparse: if set, then a wavelet sparsity is applied
; soft: if set, a soft thresholding is applied instead of the hard thresholding. Valid only if sparse is set.
; Niter: number of iterations
; TrueSolGamma -- IDL [*,*,2] = true solution in case of simu. Then print the reconstruction error at each iteration.

;print,niter
;read,idum
if not keyword_set(Niter) then Niter=20
 if not keyword_set(nscale) then nscale=5

g1 = EBMode
vs = size(g1)
N = vs[1]
N1 = 2*N
Dep = N/2
g = dblarr(N1,N1,2)
g [Dep: Dep + N -1 , Dep: Dep + N -1 , *, *] = g1
M = g * 0
M[Dep: Dep + N -1 , Dep: Dep + N -1 , *, *] = 1.
 
 nscale=5
;Niter=20
 RecGamma = eb2g(G,psize=psize)
gcf1 = star2d(RecGamma(*,*,0), nscale=nscale)
gcf2 = star2d(RecGamma(*,*,1), nscale=nscale)
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
RecGamma = eb2g(RecG,psize=psize)

for i=0,Niter-1 do begin
    Lambda =  LastThreshold  + DeltaThreshold  * (1.-erf(2.8* double(i) / double(niter)))

   if keyword_set(BZero) then RecGamma[*,*,1] = 0
   if keyword_set(BBorderZero) then RecGamma[*,*,1] =  RecGamma[*,*,1]  * M[*,*,1]
   if keyword_set(sparse) then wtsparsecstr, RecGamma, lambda, soft=soft,TabNorm=TabNorm     ; l1cstr, RecGamma, lambda   ; sparsecstr, RecGamma, Lambda
   ; load, RecGamma[*,*,1]
  ;  load, RecGamma[*,*,0]
   RecG = g2eb(RecGamma,psize=psize)
   Resi =  M*(G - RecG)
   RecG = RecG  + Resi
   RecGamma = RecGamma  + eb2g(Resi,psize=psize)
   if keyword_set(sparse) then print,  i+1, ", Lambda = ", Lambda, ", SigmaResi = ", sigma(Resi) $
   else print,  i+1, ", SigmaResi = ", sigma(Resi) 
   if keyword_set(TrueSolGamma) then  print, '          SigmaErrK = ',  sigma( RecGamma[Dep: Dep + N -1 , Dep: Dep + N -1 ] - TrueSolGamma[*,*,0] ) / sigma(TrueSolGamma[*,*,0]) * 100. 
   Lambda = Lambda - StepL
end

RecGamma = RecGamma[Dep: Dep + N -1 , Dep: Dep + N -1, * ] 
return, RecGamma
end
