;-------------------------------------------------------------------------------
;+
; NAME:
;	em_fit
;	
; PURPOSE:
;	Estimate the full set of parameters of the covariance matching model using 
;   the expectation maximization algorithm. This procedure does not allow for 
;   constraints to be set on the model parameters.
;
;
; EXPLANATION:
;   
;
; CALLING SEQUENCE:
;	em_fit, Stats, Param_in, Param_out, nb_steps=nb_steps, RESCALE = rescale
;
; INPUTS:
;	Stats		: structure grouping covariance matrices computed from the data
;					required fields are Stats.covmat and Stats.weight
;   Param_in	: initial guess of the parameter values, structure
;
; OPTIONAL INPUTS:
;   nb_steps	: number of EM steps
;
;
; KEYWORD PARAMETERS:
;   RESCALE		: boolean,  to rescale the columns of the mixing matrix at each step
;				and the corresponding lines of the componenet matrix
;
; OUTPUTS:
;	Param_out  : structure, estimated values of the model parameters after nb_steps EM steps 
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;   rescale_ap
;
; EXAMPLE:
;	em_fit, Stats, Param_in, Param_out
;
; MODIFICATION HISTORY:
;		Written: Yassir Moudden & Jean-Francois Cardoso 2005.
;
;-------------------------------------------------------------------------------   


pro em_fit, Stats, Param_in, Param_out, nb_steps=nb_steps, RESCALE = rescale, positivity = positivity

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 3 $
	then begin 
        print, 'CALLING SEQUENCE: em_fit, Stats, Param_in, Param_out, nb_steps=nb_steps,  RESCALE = rescale'
        goto, DONE
	end
	;------------------------------------------------------------------------

	;------------------------------------------------------------------------
	;set default values of optional input parameters
	if not keyword_set(nb_steps) then nb_steps = 100			



	;------------------------------------------------------------------------
	;reading structures Param_in and Stats	
	A  = Param_in.mixmat 
	Rp = Param_in.source 
	Rn = Param_in.noise 
	w  = Stats.weight 
	m = (size(A))(1)	;number of channels
	n = (size(A))(2)	;number of sources
	q = (size(Rp))(2)   ;number of data covariance matrices to be jointly matched


	;defining a few temporary matrices
	Rxs  = dblarr(m,n)
	Rss  = dblarr(n,n)
	Rxsq = dblarr(m,n)
	Rssq = dblarr(n,n)    
	Cq    = dblarr(n,n)
	Wq    = dblarr(n,m)
	iNA   = dblarr(m,n) 
	Pr    = dblarr(n,1) 
	Pr2   = dblarr(n,1) 
	Rxx  = dblarr(m,m) 
	
	for num_q=0, q-1 do Rxx  = Rxx + w(num_q) * Stats.covmat[*,*,num_q]
	
	;------------------------------------------------------------------------	
	;iterating em steps depending on noise model ie either white or colored
	case Param_in.type of
		;white noise
		1: begin
			
			;starting value of the mismatch
			smica_mismatch, Stats, Param_in, mismatch_old
			mismatch_new = mismatch_old
			Param_out = {mixmat:dblarr(m, n) , source:dblarr(n, q) , noise:dblarr(m,1) , type:1 }
			Param_out = Param_in 
			Param_old = {mixmat:dblarr(m, n) , source:dblarr(n, q) , noise:dblarr(m,1) , type:1 }
			Param_old = Param_in 
			
			for num_step = 1l, nb_steps do begin
				
				Param_old = Param_out
				mismatch_old = mismatch_new
				A  = Param_old.mixmat 
				Rp = Param_old.source 
				Rn = Param_old.noise 
			
				Rss  = 0.*dblarr(n,n)
				Rxs  = 0.*dblarr(m,n)
				iNA  = diag_matrix(1. / Rn) # A
			
				for num_q= 0, q-1 do begin
					Pr = sqrt(Rp(*,num_q))
					Pr2 = Pr # transpose(Pr) 
					Cq = Pr2 * invert( identity(n) + Pr2 * (transpose(A) # iNA) , /double) 
					Wq = iNA # Cq 
					Rxsq  = Stats.covmat(*,*,num_q)# Wq 
					Rssq = Cq + transpose(Wq) # Rxsq
					Rxs = Rxs + w(num_q) * Rxsq 
					Rss = Rss + w(num_q) * Rssq 
					Rp(*,num_q) = diag_matrix(Rssq) 
					
					;impose positivity
					if keyword_set(positivity) then begin
						Rp_tmp = Rp(*,num_q) 
						index_pos = where(Rp_tmp gt 0, count_neg)
						index_neg = where(Rp_tmp le 0, count_pos)
						if (count_neg ne 0) and (count_pos ne 0) then Rp_tmp(index_neg) = 0.1 * min(Rp_tmp(index_pos) )			
						Rp(*,num_q) = Rp_tmp																					
					endif
					
				endfor

				A = Rxs # invert(Rss, /double) 
				Rn = diag_matrix( Rxx - A # transpose(Rxs) ) 
				
				;impose positivity
					if keyword_set(positivity) then begin
						Rn_tmp = Rn
						index_pos = where(Rn_tmp gt 0, count_neg)
						index_neg = where(Rn_tmp le 0, count_pos)
						if (count_neg ne 0) and (count_pos ne 0) then Rn_tmp(index_neg) = 0.1 * min(Rn_tmp(index_pos) )				
						Rn = Rn_tmp																									
					endif
					
				if keyword_set(RESCALE) then rescale_ap, A, Rp
				
				;compute new mismatch
				Param_out.mixmat = A 
				Param_out.source = Rp 
				Param_out.noise = Rn 
				smica_mismatch, Stats, Param_out, mismatch_new


				;updating only if mismatch is finite and has decreased
				if (finite(mismatch_new) ne 1 ) OR (mismatch_new gt mismatch_old)$
				then begin
					mismatch_new = mismatch_old
					Param_out  = Param_old
				endif 

			endfor
		  end 

		;free noise
		2:begin
			;starting value of the mismatch
			smica_mismatch, Stats, Param_in, mismatch_old
			mismatch_new = mismatch_old
			
			Param_out = {mixmat:dblarr(m, n) , source:dblarr(n, q) , noise:dblarr(m,1) , type:2 }
			Param_out = Param_in 
			Param_old = {mixmat:dblarr(m, n) , source:dblarr(n, q) , noise:dblarr(m,1) , type:2 }
			Param_old = Param_in 
			
			iN = dblarr(m,m)
			ImDR = dblarr(m,m)
			Big = dblarr(n,n,m)
			BBig = dblarr(n,n,m)
			ACq = dblarr(m,n)
		
			for num_step = 1, nb_steps do begin
				Param_old = Param_out
				mismatch_old = mismatch_new
				A  = Param_old.mixmat 
				Rp = Param_old.source 
				Rn = Param_old.noise 
			
				Rxs = 0.*dblarr(m,n)
				Big = 0.*dblarr(n,n,m)

				for num_q = 0, q-1 do begin
					;EM step with respect to A and RP
					iN = diag_matrix(1. / Rn(*, num_q))
					iNA = iN # A
					Pr = sqrt( Rp(*, num_q) )
					Pr2 = Pr # transpose(Pr)
					Cq = Pr2 * invert(identity(n) + Pr2*(transpose(A)#iNA) , /double)
					Wq = iNA # Cq
					Rxsq  = Stats.covmat[*,*,num_q] # Wq 
					Rssq  = Cq + transpose(Wq) # Rxsq
					Rxs  = Rxs + w(num_q) * (iN # Rxsq) 
					;Rss  = Rss + w(num_q) * Rssq 
		
					for num_channel = 0, m-1 do BBig(*,*,num_channel) = ( w(num_q) * iN(num_channel , num_channel) ) *Rssq  
					Big = Big + BBig
					Rp(*,num_q) = diag_matrix(Rssq)		;new matrix Rp
					
					;impose positivity
					if keyword_set(positivity) then begin
						Rp_tmp = Rp(*,num_q) 
						index_pos = where(Rp_tmp gt 0, count_neg)
						index_neg = where(Rp_tmp le 0, count_pos)
						if (count_neg ne 0) and (count_pos ne 0) then Rp_tmp(index_neg) = 0.1 * min(Rp_tmp(index_pos) )			
						Rp(*,num_q) = Rp_tmp																						
					endif
					
				endfor
		
				for num_channel = 0, m-1 do begin
					iBig = invert( Big(*,*,num_channel), invalid,  /double )
					if invalid $ 
					then begin
						print, 'Invalid matrix inversion in EM_fit.pro'
						goto, DONE
					endif	
					A(num_channel, *) = iBig # transpose( Rxs(num_channel, * ) ) 
				endfor
				;new matrix A
		
				if keyword_set(RESCALE) then rescale_ap, A, Rp


				;compute new mismatch
				Param_out.mixmat = A 
				Param_out.source = Rp 
				Param_out.noise = Rn 
				smica_mismatch, Stats, Param_out, mismatch_new
				
				;updating only if mismatch is finite and has decreased
				if (finite(mismatch_new) ne 1 ) OR (mismatch_new gt mismatch_old)$
				then begin
					mismatch_new = mismatch_old
					Param_out  = Param_old
				endif 

				Param_old = Param_out
				mismatch_old = mismatch_new
				A  = Param_old.mixmat 
				Rp = Param_old.source 
				Rn = Param_old.noise 

				;EM step with respect to Rn and Rp with the new value of A 
				for num_q = 0, q-1 do begin
					iNA    = diag_matrix( 1. / Rn(*,num_q) ) # A
					;iN = diag_matrix(1. / Rn(*, num_q))
					;iNA    = iN # A
					Pr     = sqrt(Rp(*,num_q))
					Pr2    = Pr# transpose(Pr) 
					Cq     = Pr2 * invert( identity(n) + Pr2 * ( transpose(A)#iNA), /double) 
					Wq     = iNA # Cq
					Rxsq   = Stats.covmat[*,*,num_q] # Wq 
					Rssq   = Cq + transpose(Wq) # Rxsq
					Rxs    = Rxs + w(num_q) * (iN # Rxsq) 
					;Rss    = Rss + w(num_q) * Rssq 
					ACq    = A # Cq 
					ImDR   = identity(m) - ACq # transpose( iNA)		;slightly redundant
					
					Rn(*,num_q) = diag_matrix( ImDR # Stats.covmat[*,*,num_q] # transpose(ImDR) + ACq # transpose(A) )  ;new matrix Rn
					Rp(*,num_q) = diag_matrix(Rssq)		;new matrix Rp
					
					;impose positivity
					if keyword_set(positivity) then begin
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
					endif
																												
				endfor
				
				;compute new mismatch
				Param_out.mixmat = A 
				Param_out.source = Rp 
				Param_out.noise = Rn 
				smica_mismatch, Stats, Param_out, mismatch_new

				;updating only if mismatch is finite and has decreased
				if (finite(mismatch_new) ne 1 ) OR (mismatch_new gt mismatch_old)$
				then begin
					mismatch_new = mismatch_old
					Param_out  = Param_old
				endif 
			endfor
		end
		
		else:   begin
					print, 'Noise type not supported. Param_in.type should either be 1:white noise or  2:colored noise.'
					goto, DONE
				end
	endcase   
   
 

DONE:

return
end
