;-------------------------------------------------------------------------------
;+
; NAME:
;	rescale_Param
;	
; PURPOSE:
;	Rescales the mixing matrix A = Param.mixmat and the source variance profiles 
;   Rp = Param.source so that the l2 norms of the columns A(*,i) of the mixing 
;   matrix are all equal to one, while leaving unchanged
;   the matrix product A#Rp#transpose(A)
;
; EXPLANATION:   
;
; CALLING SEQUENCE:
;	rescale_Param, Param
;
; INPUTS:
;	Param : structure grouping the covariance matching model parameters
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   The output is simply the input modified
;   
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;	rescale_Param, Param
;
; MODIFICATION HISTORY:
;
;-------------------------------------------------------------------------------   

pro rescale_Param, Param

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 1 $
	then begin 
        print, 'CALLING SEQUENCE: rescale_Param, Param'
        goto, DONE
	end
	;------------------------------------------------------------------------



	;------------------------------------------------------------------------
	A = Param.mixmat  
	Rp = Param.Source 
	l2normA =total( A*A, 1, /double) 
	Rp = diag_matrix(l2normA ) # Rp 
    A = A # diag_matrix(1./ sqrt(l2normA) ) 
	Param.mixmat = A  
	Param.Source = Rp


DONE:  
  
return 
end



