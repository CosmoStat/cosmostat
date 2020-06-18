;-------------------------------------------------------------------------------
;+
; NAME:
;   qn_fit
;	
; PURPOSE:
;	uses quadratic approximation and fixed point to optimize mismatch wrt to noise spectra
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;  qn_fit, Stats, Param_in, Param_out, nb_steps = nb_steps, epsilon = epsilon
; 
; INPUTS:
;	Stats		: structure grouping covariance matrices computed from the data
;					required fields are Stats.covmat and Stats.weight
;   Param_in	: initial guess of the parameter values, structure
;
; OPTIONAL INPUTS:
;
;   nb_steps	: maximum number of iterations
;
; KEYWORD PARAMETERS:
;  
;
; OUTPUTS:
;	
;   Param_out  : structure, estimated values of the model parameters after nb_steps EM steps 
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;		
;
; EXAMPLE:
;  qn_fit, Stats, Param_in, Param_out
;
; MODIFICATION HISTORY:
; Yassir Moudden & Jean-Francois Cardoso, 2005
;
;----------------------------------------------------------------------

pro qn_fit, Stats, Param_in, Param_out, nb_steps = nb_steps, epsilon = epsilon, positivity = positivity

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 3 $
	then begin 
        print, 'CALLING SEQUENCE: qn_fit, Stats, Param_in, Param_out, nb_steps = nb_steps, epsilon = epsilon'
        goto, DONE
	end
	;------------------------------------------------------------------------


	;------------------------------------------------------------------------
	;set default values of optional input parameters
	
	if not keyword_set(nb_steps) then nb_steps = 3			;this a fixed point algorithm. nb_iter is a maximum number of iterations
	if not keyword_set(epsilon) then epsilon =0.000001			;relative threshold on the decrease of the mismatch



	;------------------------------------------------------------------------
	;reading structures Param_in and Stats	
	A  = Param_in.mixmat 
	Rp = Param_in.source 
	Rn = Param_in.noise 
	w  = Stats.weight 
	m = (size(A))(1)	;number of channels
	n = (size(A))(2)	;number of sources
	q = (size(Rp))(2)   ;number of data covariance matrices to be jointly matched

	
	;------------------------------------------------------------------------
	;optimization

	case Param_in.type of			1: begin				;in the white noise case
	
			
			;initialisation
			Param_out = {mixmat:dblarr(m, n) , source:dblarr(n, q) , noise:dblarr(m,1) , type:1 }
			Param_out = Param_in 
			iR   = dblarr(m,m)
			AiRA = dblarr(n,n) 
      	
			smica_mismatch, Stats, Param_in, mismatch_old
			mismatch_new = mismatch_old
			Rn_new = Rn
			num_step = 1
			temp_Rn = dblarr(m ,q)					;temporary variable
			
			;iterations
			while (num_step LE nb_steps) do begin
			
					mismatch_old = mismatch_new
					Rn = Rn_new

					for num_q=0, Q-1 do begin				;BEGIN loop over the covariance matrices
						
						hR = Stats.covmat(*,*,num_q)
						ARp = A # diag_matrix(Rp(*,num_q)) 
						R  = ARp # transpose(A) + diag_matrix(Rn) 
						R = 0.5 * (R + transpose(R))
						
						iR = invert(R, invalid,/double)
						if invalid $
						then begin
							print, 'Invalid matrix inversion in qp_fit.pro'
							goto, DONE
						endif	
											
						temp_Rn(*,num_q) = invert( iR * iR, /double) # diag_matrix(iR # (hR - A#diag_matrix(Rp(*,num_q)) # transpose(A)  ) # iR) 
						
					endfor								;END loop over the covariance matrices
					
					
					Rn_new = temp_Rn#transpose(w)		;averaging over the bands			
					
					;impose positivity
					if keyword_set(positivity) then begin
						Rn_new_tmp = Rn_new(*,num_q) 
						index_pos = where(Rn_new_tmp gt 0, count_neg)
						index_neg = where(Rn_new_tmp le 0, count_pos)
						if (count_neg ne 0) and (count_pos ne 0) then Rn_new_tmp(index_neg) = 0.1 * min(Rn_new_tmp(index_pos) )  
						Rn_new(*,num_q) = Rn_new_tmp
					endif
					
					Param_out.noise = Rn_new
					smica_mismatch, Stats, Param_out, mismatch_new
						
					;updating only if mismatch is finite and has decreased
					if (finite(mismatch_new) ne 1 ) OR (mismatch_new gt mismatch_old)$
					then begin
						mismatch_new = mismatch_old
						Param_out.noise = Rn
						Rn_new = Rn
					endif 

					if (mismatch_new gt (1.-epsilon)*mismatch_old) $
					then  num_step = nb_steps +1 $
					else num_step = num_step+1
			endwhile
				
			end		;of the white noise case
	


			2: begin  ;the free noise case				
		
			;initialisation
			Param_out = {mixmat:dblarr(m, n) , source:dblarr(n, q) , noise:dblarr(m,1) , type:2 }
			Param_out = Param_in 
			iR   = dblarr(m,m)
      	
			smica_mismatch, Stats, Param_in, mismatch_old
			mismatch_new = mismatch_old
			Rn_new = Rn
			num_step = 1
			
			;iterations
			while (num_step LE nb_steps) do begin
			
					for num_q=0, Q-1 do begin				;BEGIN loop over the covariance matrices

						mismatch_old = mismatch_new
						Rn = Rn_new
					
						hR = Stats.covmat(*,*,num_q)
						ARp = A # diag_matrix(Rp(*,num_q)) 
						R  = ARp # transpose(A) + diag_matrix(Rn(*,num_q)) 
						R = 0.5 * (R + transpose(R))
						
						iR = invert(R, invalid,/double)
						if invalid $
						then begin
							print, 'Invalid matrix inversion in qp_fit.pro'
							goto, DONE
						endif	
											
						Rn_new(*,num_q) = invert( iR * iR, /double) # diag_matrix(iR # (hR - A#diag_matrix(Rp(*,num_q)) # transpose(A)  ) # iR) 
	
						;impose positivity
						if keyword_set(positivity) then begin
							Rn_new_tmp = reform( Rn_new(*,num_q) )
							index_pos = where(Rn_new_tmp gt 0, count_pos)
							index_neg = where(Rn_new_tmp le 0, count_neg)
							if (count_neg ne 0) and (count_pos ne 0) then Rn_new_tmp(index_neg) = 0.1 * min(Rn_new_tmp(index_pos) )  
							Rn_new(*,num_q) = Rn_new_tmp
						endif
					
						Param_out.noise = Rn_new
						smica_mismatch, Stats, Param_out, mismatch_new
						
						;updating only if mismatch is finite and has decreased
						if (finite(mismatch_new) ne 1 ) OR (mismatch_new gt  mismatch_old)$
						then begin
							mismatch_new = mismatch_old
							Param_out.noise = Rn
							Rn_new = Rn
						endif

					endfor								;END loop over the covariance matrices
	
					if (mismatch_new gt (1.-epsilon)*mismatch_old) $
					then  num_step = nb_steps +1 $
					else num_step = num_step+1
					
				endwhile
				
				end		;of the free noise case
				
	endcase



DONE:
return
end























