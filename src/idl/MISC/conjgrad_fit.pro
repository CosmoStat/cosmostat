;-------------------------------------------------------------------------------
;+
; NAME:
;   conjgrad_fit
;	
; PURPOSE:
;	use the optimisation library in IDL to perform conjugate gradient optimization of the smica mismatch
;over all parameters.
;
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;  conjgrad_fit, Stats, Param_in, Param_out
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
;		positivity : if the key word is set, the parameters are reprojected before output to make 
;					the source and noise varaince parameters positive.
;					(maybe one might want to try  a projection after each gradient step)
; OUTPUTS:
;		Param_out : structure, groups model parameter values after conjugate gradient steps
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;		kron, smica_mismatch,minF_conj_grad, conjgrad_smica_mismatch,Struct2Vec, Vec2Struct
;
; EXAMPLE:
;  conjgrad_fit, Stats, Param_in, Param_out 
;
; MODIFICATION HISTORY:
;       Written: Yassir Moudden 2005.
;       September, 2005 File creation
;
;-------------------------------------------------------------------------------   

pro conjgrad_fit, Stats, Param_in, Param_out, positivity = positivity

COMMON SMICAenv

;initializing global variables in the SMICAenv. NB these variables are first defined in mrs.pro
Statistics = Stats
A = Param_in.mixmat 
Rp = Param_in.source
number_c = ( size(A))(1) 
number_s = ( size(A))(2)
number_b = ( size(Rp))(2) 
	

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 3 $
	then begin 
        print, 'CALLING SEQUENCE:  conjgrad_fit, Stats, Param_in, Param_out '
        goto, DONE
	end
	
	if not (  Param_in.type EQ 2 )$
	then begin
		print, 'Only parameter type 2 , id est free noise, is supported.'
		goto, DONE
	end
	
	
	;------------------------------------------------------------------------

	;convert the parameter structure to a vector
	;this is our initial guess for the optimization procedure
		
		vec_0 = Struct2Vec( Param_in ) 
		misma_0 = conjgrad_smica_mismatch( vec_0 )
		vec_1 = vec_0
		misma_1 = misma_0
	
		minF_conj_grad, vec_1, misma_1, conv_factor, FUNC_NAME="conjgrad_smica_mismatch", /INITIALIZE
	
		;imposing positivity after each gradient step
		if keyword_set(positivity) then begin
			param_out = Vec2Struct(vec_1)
			Rp = param_out.source
			Rn = param_out.noise
			for num_q = 0, number_b -1 do begin
				Rp_tmp = Rp(*,num_q) 
				index_pos = where(Rp_tmp gt 0, count_neg)
				index_neg = where(Rp_tmp le 0, count_pos)
				if (count_neg ne 0) and (count_pos ne 0) then Rp_tmp(index_neg) = 0.1 * min(Rp_tmp(index_pos) )			
				Rp(*,num_q) = Rp_tmp																					
				Rn_tmp = Rn
				index_pos = where(Rn_tmp gt 0, count_neg)
				index_neg = where(Rn_tmp le 0, count_pos)
				if (count_neg ne 0) and (count_pos ne 0) then Rn_tmp(index_neg) = 0.1 * min(Rn_tmp(index_pos) )				
				Rn = Rn_tmp																									
			endfor 
			param_out.source = Rp
			param_out.noise = Rn
			vec_1 = Struct2Vec( param_out ) 
			misma_1 = conjgrad_smica_mismatch( vec_1 )
		endif

		if (finite(misma_1) ne 1) OR ( misma_1 gt misma_0)$
		then begin 
			vec_1 = vec_0
			misma_1 = misma_0
			goto, DONE
		endif else begin
			misma_0 = misma_1
			vec_0 = vec_1
		endelse

		num_iter = 1.
		max_number_iter = 2
		while (num_iter LT max_number_iter) do begin
			
			vec_1 = vec_0
			misma_1 = misma_0
			minF_conj_grad, vec_1, misma_1, conv_factor, FUNC_NAME="conjgrad_smica_mismatch"			
			
			;imposing positivity after each gradient step
			if keyword_set(positivity) then begin
				param_out = Vec2Struct(vec_1)
				Rp = param_out.source
				Rn = param_out.noise
				for num_q = 0, number_b -1 do begin
					Rp_tmp = Rp(*,num_q) 
					index_pos = where(Rp_tmp gt 0, count_neg)
					index_neg = where(Rp_tmp le 0, count_pos)
					if (count_neg ne 0) and (count_pos ne 0) then Rp_tmp(index_neg) = 0.1 * min(Rp_tmp(index_pos) )			
					Rp(*,num_q) = Rp_tmp																					
					Rn_tmp = Rn
					index_pos = where(Rn_tmp gt 0, count_neg)
					index_neg = where(Rn_tmp le 0, count_pos)
					if (count_neg ne 0) and (count_pos ne 0) then Rn_tmp(index_neg) = 0.1 * min(Rn_tmp(index_pos) )				
					Rn = Rn_tmp																									
				endfor 
				param_out.source = Rp
				param_out.noise = Rn
				vec_1 = Struct2Vec( param_out ) 
				misma_1 = conjgrad_smica_mismatch( vec_1 )
			endif
			

			if (finite(misma_1) ne 1) OR ( misma_1 gt misma_0)$
			then begin 
				vec_1 = vec_0
				misma_1 = misma_0
				goto, DONE
			endif else begin
				misma_0 = misma_1
				vec_0 = vec_1
			endelse
			
			num_iter = 1. + num_iter
						
		endwhile

DONE:
	vec_0 = vec_1
	misma_0 = misma_1
	param_out = Vec2Struct(vec_0)
	
return
end


;-------------------------------------------------------------------------------   
pro smica_mismatch_bis, Param, mismatch, mismatch_grad, mismatch_hess

COMMON SMICAenv
 q = number_b
 n = number_s
 m = number_c
 

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 2 $
	then begin 
        print, 'CALLING SEQUENCE: smica_mismatch_bis, Param, mismatch, mismatch_grad, mismatch_hess'
        goto, DONE
	endif
	
	;checking that the noise type is supported
	if not ( Param.type EQ 2 )$
	then begin
		print, "Noise parameter type not supported."
		goto, DONE
	endif

	
	;checking whether the mismatch gradient or and hessian wrt to A should be computed
	if (N_PARAMS(0) EQ 3) or (N_PARAMS(0) EQ 4)  then Compgrad = 1 else Compgrad = 0
	if (N_PARAMS(0) EQ 4) then Comphess = 1 else Comphess = 0
	;------------------------------------------------------------------------
	
	
	;------------------------------------------------------------------------
	;reading the Statistics and Param structures
	A = Param.mixmat     
	Rp = Param.source
	Rn = Param.noise
	w = Statistics.weight
	;allocating a few temporary matrices
	R    = dblarr(m,m) 
	dR   = dblarr(m,m) 
	hR   = dblarr(m,m) 
	iR   = dblarr(m,m) 
	mm   = dblarr(1,q) 
	ARp  = dblarr(m, n)
	iRARp = dblarr(m,n)
	Im   = identity(m) 
	
	
	if (Compgrad EQ 1) then mismatch_grad = {mixmat:dblarr(m, n) , source:dblarr(n, q) , noise:dblarr(m,q)}
	if (Comphess EQ 1) then begin
		mismatch_hess = {mixmat:dblarr(m*n, m*n) , PN:dblarr((n+m)*q, (n+m)*q ) }
		BB = [[A], [Im]]
		HPN = dblarr((n+m)*q, (n+m)*q)
		h1 = dblarr(m*n, m*n)
		h2 = dblarr(m*n, m*n)
	endif

	;------------------------------------------------------------------------
	for num_q=0, q-1 do begin
		hR = Statistics.covmat(*,*,num_q)
		ARp = A # diag_matrix(Rp(*,num_q)) 
		R  = ARp # transpose(A) + diag_matrix(Rn(*,num_q)) 
			
		iR = invert(R, invalid,/double)
		if invalid $
		then begin
			print, 'Invalid matrix inversion while computing mismatch.'
			goto, DONE
		endif	
			
		dR = hR # iR 
		;mm(num_q) = 0.5 * ( trace(dR, /double) - alog(determ(dR, /double)) - m ) 
		
		val_hR = EIGENQL( 0.5*(hR + transpose(hR) ) , /DOUBLE, EIGENVECTORS=vec_hR)
		val_R =  EIGENQL( 0.5*(R +transpose(R) ) , /DOUBLE, EIGENVECTORS=vec_R)
		

		
		determ_dR = abs( product(val_hR, /double) / product( val_R, /double) )
		trace_dR = trace( hR # vec_R # diag_matrix(1./val_R) #transpose( vec_R )  , /double )

		sign_val = total( sign(val_R) )+ total( sign(val_hR) )
		if ( sign_val eq 2.*double(m) ) $
		then mm(num_q)= 0.5 * ( trace_dR - alog(abs(determ_dR))- double(m) ) $
		else mm(num_q) = alog(-1.)

		if (Compgrad EQ 1) $
		then begin
			dRbis = iR # (dR-Im) 
			mismatch_grad.mixmat = mismatch_grad.mixmat - w(num_q) * ( dRbis #ARp)			
			mismatch_grad.source(*, num_q) = - w(num_q)*0.5 * diag_matrix( transpose(A )#dRbis#A )
			mismatch_grad.noise(*, num_q) = - w(num_q)*0.5 * diag_matrix( dRbis )
		endif
		
		if (Comphess EQ 1) $		
		then begin
			;hessian with respect to the mixing matrix		
			iRARp = iR # ARp
			kron, transpose(iRARp)#ARp, iR, kron1
			kron, transpose(iRARp), iRARp, kron2
			h1 = h1 + w(num_q) * kron1
			h2 = h2 + w(num_q) * kron2
			
			;hessian with respect to source and noise
			HPN((n+m)*num_q:(n+m)*(num_q+1) -1 , (n+m)*num_q:(n+m)*(num_q+1) -1) = 0.5 * w(num_q) * (  (transpose(BB)#iR# BB)^2     )
		endif		

	endfor
	
	mismatch = total( w * mm )
		
	
	if (Comphess EQ 1) $
	then begin
		Idx_A = findgen(n,m )		
		Idx_A = transpose(Idx_A)
		Idx_A = Idx_A(*) 
		mismatch_hess.mixmat = h1 + h2(*, Idx_A) 

	
		Idx_PN = transpose( findgen(n+m, q ) )		
		Idx_PN = Idx_PN(*) 
		HPN = HPN(Idx_PN,*)
		HPN = HPN(*, Idx_PN)
		mismatch_hess.PN = HPN
	endif		

DONE:	
	
end	

;-------------------------------------------------------------------------------   
function conjgrad_smica_mismatch, vec_param, gradient

COMMON SMICAenv
q = number_b
n = number_s
m = number_c
 
	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if (N_PARAMS() LT 1 )  or (N_PARAMS() gT 2) $
	then begin 
        print, 'CALLING SEQUENCE:  conjgrad_smica_mismatch, vec_param, gradient'
        goto, DONE
	end
	
	if N_PARAMS() eq 2 then compute_gradient = 1 else compute_gradient = 0
	
	
	;------------------------------------------------------------------------


	;compute the mismatch and the derivatives
	param = Vec2Struct(vec_param) 
	if compute_gradient then begin 
	smica_mismatch_bis, param, mismatch, gradient_mismatch, hessian_mismatch
	endif else begin													
	smica_mismatch_bis, param, mismatch				
	endelse															
	
	
		if compute_gradient then begin
			gradient = Struct2Vec(gradient_mismatch)
			iHA = invert(hessian_mismatch.mixmat, /double, invalid_A)
			iHPN = invert(hessian_mismatch.PN, /double, invalid_PN)
			if invalid_A or invalid_PN then begin
					print, 'invalid hessian inverse'
					print, 'gradient on output is not corrected by the hessian'
					goto, DONE
			endif else begin
				temp = gradient
				temp(0:n*m-1) = iHa#gradient(0:n*m-1)
				temp(n*m:n*m+(n+m)*q -1) = iHPN#gradient(n*m:n*m+(n+m)*q -1)
				gradient = temp
			endelse	
		endif 
	

DONE:	
	
return, mismatch	
end	

;-------------------------------------------------------------------------------   
function Struct2Vec, param 

COMMON SMICAenv
 q = number_b
 n = number_s
 m = number_c
 


theta = dblarr(m*n + n*q + m*q) 
theta(0:m*n-1) = param.mixmat(*)
temp_source = transpose( param.source )
theta(m*n:n*q+m*n-1) = temp_source(*)
temp_noise = transpose(param.noise)
theta(n*q+m*n:n*q+m*n+ m*q-1) = temp_noise(*)

return, theta
end

;-------------------------------------------------------------------------------   

function Vec2Struct, theta 

COMMON SMICAenv
 q = number_b
 n = number_s
 m = number_c
 

param = {mixmat:dblarr(m, n) , source:dblarr(n, q) , noise:dblarr(m,q) , type:2 }

param.mixmat(*) = theta(0:m*n-1)

temp_source = dblarr(q,n)
temp_source(*) = theta(m*n:n*q+m*n-1)
param.source = transpose( temp_source ) 

temp_noise = dblarr(q,m)
temp_noise(*) =theta(n*q+m*n:n*q+m*n+ m*q-1)  
param.noise  =  transpose( temp_noise )

return, param
end

;-------------------------------------------------------------------------------   

function GradFisherNorm, mismatch_grad, mismatch_hess 

COMMON SMICAenv
 q = number_b
 n = number_s
 m = number_c
 
 	
	gradient = Struct2Vec(mismatch_grad)
	iHA = invert(mismatch_hess.mixmat, /double, invalid_A)
	iHPN = invert(mismatch_hess.PN, /double, invalid_PN)
	if invalid_A or invalid_PN then begin
			print, 'invalid hessian inverse'
			gfn = alog(-1)
	endif else begin
		gfn = 0.
		gfn = gfn + transpose(gradient(0:n*m-1))#iHa#gradient(0:n*m-1)
		gfn = gfn + transpose(gradient(n*m:n*m+(n+m)*q -1))#iHPN#gradient(n*m:n*m+(n+m)*q -1)
	endelse

DONE:
return, gfn
end
	

