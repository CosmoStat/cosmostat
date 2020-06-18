;-------------------------------------------------------------------------------
;+
; NAME:
;	rescale_ap
;	
; PURPOSE:
;	Rescales the mixing matrix A and the source variance profiles Rp so that the 
;   l2 norms of the columns A(*,i) of the mixing matrix are all equal to one, while 
;   leaving unchanged the matrix product A#Rp#transpose(A)
;
; EXPLANATION:   
;
; CALLING SEQUENCE:
;	rescale_ap, A, Rp
;
; INPUTS:
;	A : mixing matrix
;   Rp: source covariance profiles
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;   The outputs are simply the inputs modified
;   
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;	rescale_ap, A, Rp
;
; MODIFICATION HISTORY:
;
;-------------------------------------------------------------------------------   

pro rescale_ap, A, Rp

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 2 $
	then begin 
        print, 'CALLING SEQUENCE: rescale_ap, A, Rp'
        goto, DONE
	end
	;------------------------------------------------------------------------



	;------------------------------------------------------------------------
	l2normA =total( A*A, 1, /double) 
	Rp = diag_matrix(l2normA ) # Rp 
    A = A # diag_matrix(1./ sqrt(l2normA) ) 

DONE:  
  
return 
end
