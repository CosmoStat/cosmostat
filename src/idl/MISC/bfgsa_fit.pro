;-------------------------------------------------------------------------------
;+
; NAME:
;   bfgsa_fit
;	
; PURPOSE:
;	uses bfgs algorithm to optimise mismatch with respect to all the 
;   the entries of the mixing matrix A only.
;
;   this is mostly intended to speed up convergence in
;   the colored noise case. Indeed, some noise variance 
;   parameters can become very small and the em_fit procedure
;   is no longer fast enough.
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;  bfgsa_fit, Stats, Param_in, Param_out
; 
; INPUTS:
;
;   Stats : structure, groups data statistics
;			required fields are Stats.covmat and Stats.weight
;	Param_in : structure, groups model parameter values
;
; COMMENT:
;		future version will have a Mask as optional input to specify over which 
;		subset of entries of the mixing matrix optimization is to be conducted
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;		Param_out : structure, groups model parameter values after bfgsa step
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;		kron, smica_mismatch
;
; EXAMPLE:
;  bfgsa_fit, Stats, Param_in, Param_out 
;
; MODIFICATION HISTORY:
;
;		Written: Yassir Moudden & Jean-Francois Cardoso 2005.
;-------------------------------------------------------------------------------   

pro bfgsa_fit, Stats, Param_in, Param_out


	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 3 $
	then begin 
        print, 'CALLING SEQUENCE:  bfgsa_fit, Stats, Param_in, Param_out '
        goto, DONE
	end
	
	if not ( Param_in.type EQ 1 or Param_in.type EQ 2 )$
	then begin
		print, 'Parameter type not supported.'
		goto, DONE
	end
	
	;------------------------------------------------------------------------
	
	;-----------------------------------------------------------------------------
	; Defining a few constants
	max_nb_iter   = 100.  
	epsilon  = 1.0e-10 
	;stops when the gradient's Fisher-norm is smaller 
	
	; Backtracking  parameters
	max_backtrack = 20		;Maximum number of backtrack steps in each line search
	C1           = 0.01		;First Wolfe  condition
	RHO          = 0.66		;backtracking factor in line search
	;-----------------------------------------------------------------------------



	;-----------------------------------------------------------------------------
	;reading input structures Stats and Param_in
	A = Param_in.mixmat 
	Rp = Param_in.source
	Rn = Param_in.noise 
	w = Stats.weight

	m = ( size(A))(1) 
	n = ( size(A))(2)
	q = ( size(Rp))(2) 
	
	;extends the white noise case to be the same as the colored noise case
	if ( Param_in.type EQ 1)  then Rn = Rn # (dblarr(1, q) + 1.) 

	
	;-----------------------------------------------------------------------------
	;defining a few temporary matrices, variables
	Al = dblarr(m,n)
	grad = dblarr(m,n)  
	oldgrad = dblarr(m,n) 
	sdir = dblarr(m,n)   ;direction of line search
	dX = dblarr(m*n,1)  
	dG = dblarr(m*n,1) 
	iHdG = dblarr(m*n,1)  
	diH = dblarr(m*n)   

	val    = 0. 
	oldval = 0. 
	gain = 1. 
	
	Param_out = Param_in  
	
	;-----------------------------------------------------------------------------
	;initializing the hessian and gradient
	smica_mismatch, Stats, Param_in, val, grad, hess
	initiH = invert(hess, invalid, /double)
	if invalid then begin
		print, 'Invalid matrix inversion.'
		goto, Done
	endif
	iH = initiH
	
	;-----------------------------------------------------------------------------
	for num_iter=1, max_nb_iter do begin
		
		sdir(*)   = -iH#grad(*)      
		;Search along the rectfied gradient
		
		dirder    = transpose( grad(*)) # sdir(*) 
		;The derivative in the search direction. 
			
		if (-dirder  LT epsilon) or (gain LT 1.0e-6 ) then break 

		if (dirder>0) then begin 
			print, 'Resetting the inverse Hessian.'		; why not a new hessian??
			iH      = initiH 
			sdir(*) = -iH#grad(*)
			dirder  = transpose(grad(*)) # sdir(*)  
		endif

		oldval  = val 
		oldgrad = grad 
    
		; line search, backtracking
		lambda = 1. 
		num_bt    = max_backtrack 
		Al     = A + lambda * sdir  
		Param_out.mixmat = Al
		smica_mismatch, Stats, Param_out, val1, grad1

		while ( (num_bt-1) AND  (val1 GT (val + C1 * lambda * dirder)) ) do begin
			lambda = RHO * lambda 
			Al     = A + lambda * sdir  
			Param_out.mixmat = Al
			smica_mismatch, Stats, Param_out, val1, grad1
		endwhile


	
		;-----------------------------------------------------------------------------
		;updating the mixing matrix A with the line search result		
		if (num_bt EQ 0) then begin 
			A = A
			Param_out.mixmat= A 
			print, "Problem in line search: too much backtracking."
			print, "Doing without this optimisation step."
		endif else begin
			A = Al
			Param_out.mixmat = A
		endelse
		
		smica_mismatch, Stats, Param_out, val, grad
		gain = oldval-val		
		;-----------------------------------------------------------------------------
		;updating the hessian using the bfgs rule	
		;why not use the approximate hessian computed with smica_mismatch??	
		;there should be an option to use either possibility
		
		dX   = lambda * sdir(*)  
		dG   = grad(*) - oldgrad(*) 
		idGdX  = (1./(transpose(dG)#dX))(0)		
		iHdG   = iH#dG 
		dGiHdG = (0.5*(1.+idGdX*(transpose(dG)#iHdG))   )(0)
		diH    = (idGdX*dX) # ( transpose( ( dGiHdG*dX  - iHdG) ))  
		iH     = iH + ( diH + transpose(diH) )  
    
	endfor

	
DONE:

return
end



