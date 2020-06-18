;+
; NAME: 
;				VF_WT_THRESHOLDING
;
; PURPOSE: 
;				Appy a wavelet thresholding to a vector field. The thresholding is applied independantly on the two components.
;
; CALLING:
;				vf_wt_thresholding, VectorField, lambda, nscale=nscale, soft=soft, TabNorm=TabNorm
;
; INPUT/OUTPUT: 
;				VectorField ---- IDL array [*,*,0:1] : shear or EB data set
;               TabNorm ---- IDL array [0:nscale-1] : l2 normalization table for the wavelet coefficients
;
; INPUTS
;                 Lambda    ---- float: threshold level 
;
; INPUT KEYWORD:			
;              nscale --- int: number of scales used in the wavelet thresholding 
;				soft -- scalar : if set, a soft thresholding is applied instead of the hard thresholding
;
; HISTORY:
;				Written by J.-L. Starck, Oct 2013
;-
;-------------------------------------------------------------------------------

pro vf_wt_thresholding, r, lambda, nscale=nscale, soft=soft,TabNorm=TabNorm
vs = size(u)
n = vs(1)
u1 = star2d(r(*,*,0), nscale=nscale, /gen)
u2 = star2d(r(*,*,1), nscale=nscale, /gen)
 
 if not keyword_set(TabNorm) then begin
	dirac = double(r[*,*,0])
	dirac [*] = 0
	vs = size(dirac)
	Nx = vs[1]
	Ny = vs[2]
	dirac[Nx/2, Ny/2] = 1
	w  = star2d(dirac, Nscale=Nscale,  /gen)
	TabNorm = dblarr(Nscale)
	for j=0, Nscale-1 do TabNorm[j] = sqrt( total(w[*,*,j]^2) )
	MaxL = max( [max(gcf1[*,*,0:nscale-2]), max(gcf2[*,*,0:nscale-2])])
    Lambda = MaxL
end   

 for j=0, nscale-2 do begin
    z1 = u1[*,*,j]
	; softthreshold, z1, lambda
	z1 = mrs_absthreshold(z1, lambda*TabNorm[j],soft=soft)
	u1[*,*,j] = z1
	z2  = u2[*,*,j]
	; softthreshold, z2,lambda
	z2 = mrs_absthreshold(z2, lambda*TabNorm[j], soft=soft)
	u2[*,*,j] = z2
end
;u1 =  cfthres(u1, lambda)
;u2 =  cfthres(u2, lambda)

; r(*,*,0) = float ( dfti(u1) )
; r(*,*,1) = float ( dfti(u2) )

r(*,*,0) = istar2d(u1, /gen)
r(*,*,1) = istar2d(u2, /gen)
end


