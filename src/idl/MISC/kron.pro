;-------------------------------------------------------------------------------
;+
; NAME:
;	kron
;	
; PURPOSE:
;	procedure to compute the kronecker product of matrices A and B
;
; EXPLANATION:
;
; CALLING SEQUENCE:
;	kron, A, B, C
;
; INPUTS:
;	A , B : two matrices i.e. 2D arrays
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	C : matrix equal to the kronecker product of A and B
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;
; EXAMPLE:
;	kron, A, B, C 
;
; MODIFICATION HISTORY:
;   File written january 2005, Yassir Moudden
;-------------------------------------------------------------------------------   
pro kron, A, B, C

	;------------------------------------------------------------------------
	;checking for proper number of i/o parameters
	if N_PARAMS() LT 3 $
	then begin 
        print, 'CALLING SEQUENCE:  kron, A, B, C'
        goto, DONE
	end
	;------------------------------------------------------------------------


	;------------------------------------------------------------------------
	na = (size(A))(1)
	nb = (size(B))(1)   
	ma = (size(A))(2)
	mb = (size(B))(2)


	C = dblarr(na*nb, ma*mb)

	for ia = 0, na-1 do begin
		for ja = 0,ma-1 do begin
					C( ia*nb : (ia+1)*nb -1  ,  ja*mb : (ja+1)*mb - 1   ) = A(ia,ja) * B  
		endfor
	endfor

DONE:


return
end

