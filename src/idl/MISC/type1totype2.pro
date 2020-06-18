;-------------------------------------------------------------------------------
;+
; NAME:
;	type1totype2
;	
; PURPOSE:
;	changes a Param structure from type 1 to type 2
;
;   type 1 = white noise
;   type 2 = colored noise
;
; EXPLANATION:
;   simply modifies the fields Param.type and Param.noise
;   if input is already of type 2 , Param_out = Param_in
;
; CALLING SEQUENCE:
;	type1totype2, Param_in, Param_out
;
; INPUTS:
;	Param_in : structured model parameters, field Param_in.type = 1 and field Param_in.noise is an m*1 matrix
;				where m is the number of channels
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	Param_out : structured model parameters, field Param_out.type = 2 and field Param_out.noise is an m*q matrix
;				where m is the number of channels and q is the number of covariances to be jointly matched.
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;	type1totype2, Param_in, Param_out 
;
; MODIFICATION HISTORY:
;	
;-------------------------------------------------------------------------------   

pro type1totype2, Param_in, Param_out

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 2 $
	then begin 
        print, 'CALLING SEQUENCE:  type1totype2, Param_in, Param_out'
        goto, DONE
	end
	
	if Param_in.type GT 2 $
	then begin
		print, 'Incorrect parameter type.'
		goto, DONE
	endif

	;------------------------------------------------------------------------

	
	;------------------------------------------------------------------------
	;reading the inpt structure
	if Param_in.type eq 2 then Param_out = Param_in $
	else begin
		q = (size(Param_in.source[0,*]))(2) 
		m = (size(Param_in.mixmat[*,0]))(1)
		n = (size(Param_in.mixmat[0,*]))(2)

		Param_out = {mixmat:dblarr(m, n) , source:dblarr(n, Q) , noise:dblarr(m,Q) , type:2}

		Param_out.mixmat = Param_in.mixmat 
		Param_out.source  = Param_in.source  
		Param_out.noise = (Param_in.noise) # (dblarr(1,q)+1.)
		Param_out.type = 2

	endelse

DONE:

return
end	
