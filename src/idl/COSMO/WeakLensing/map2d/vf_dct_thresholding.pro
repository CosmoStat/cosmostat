;+
; NAME: 
;				VF_DCT_THRESHOLDING
;
; PURPOSE: 
;				Appy a DCT thresholding to a vector field. The thresholding is applied independantly on the two components.
;
; CALLING:
;				vf_dct_thresholding, VectorField, lambda,  
;
; INPUT/OUTPUT: 
;				VectorField ---- IDL array [*,*,0:1] : shear or EB data set
;
; INPUTS
;                 Lambda    ---- float: threshold level 
;
; HISTORY:
;				Written by J.-L. Starck, Oct 2013
;-
;-------------------------------------------------------------------------------


pro vf_dct_thresholding, r, lambda
vs = size(u)
n = vs(1)

u1 = dct(r(*,*,0))
u2 = dct(r(*,*,1))
 
softthreshold, u1, lambda
softthreshold, u2, lambda

;u1 =  cfthres(u1, lambda)
;u2 =  cfthres(u2, lambda)

; r(*,*,0) = float ( dfti(u1) )
; r(*,*,1) = float ( dfti(u2) )

r(*,*,0) = float ( dct(u1,/inverse) )
r(*,*,1) = float ( dct(u2, /inverse) )
end

